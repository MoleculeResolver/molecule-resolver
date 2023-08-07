from datetime import datetime
import sqlite3
from typing import Optional, Sequence
import threading
import uuid

from moleculeresolver.molecule import Molecule

class SqliteMoleculeCache:
    def __init__(self, db_path: str = ':memory:', expiration_datetime: Optional[datetime] = None):
        self.db_path = db_path
        self.expiration_datetime = expiration_datetime
        self._connections = {}
        self._main_thread_id = threading.get_ident()

    def __enter__(self):
        self._create_tables()
        self.delete_expired()
        return self

    def __exit__(self, *args, **kwargs):
        # close all other connections except the main one
        for thread_id, thread_connection in self._connections.items():
            if thread_id != self._main_thread_id:
                if thread_connection:
                    thread_connection.close()
                    self._connections[thread_id] = None

        # close the connection from the main thread
        this_thread_id = threading.get_ident()
        if this_thread_id  == self._main_thread_id:
            main_thread_connection = self._connections[self._main_thread_id]
            main_thread_connection.execute('PRAGMA analysis_limit=8192')
            main_thread_connection.execute('PRAGMA optimize')
            main_thread_connection.close()
            self._connections.clear()

    def get_connection(self):
        thread_id = threading.get_ident()
        if thread_id not in self._connections:
            self._connections[thread_id] = sqlite3.connect(self.db_path, check_same_thread=False)
            self._connections[thread_id].execute("PRAGMA foreign_keys = 1")
            self._connections[thread_id].execute('PRAGMA journal_mode=WAL')
            self._connections[thread_id].execute('PRAGMA synchronous=NORMAL')

        return self._connections[thread_id]

    def _create_tables(self):
        this_thread_connection = self.get_connection()
        with this_thread_connection:
            this_thread_connection.execute('''
                CREATE TABLE IF NOT EXISTS molecules (
                    id INTEGER PRIMARY KEY,
                    service TEXT NOT NULL,
                    identifier_mode TEXT NOT NULL,
                    identifier TEXT NOT NULL,
                    SMILES TEXT,
                    additional_information TEXT,
                    datetime_added DATETIME DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            this_thread_connection.execute('''
                CREATE TABLE IF NOT EXISTS synonyms (
                    id INTEGER PRIMARY KEY,
                    molecule_id INTEGER NOT NULL,
                    synonym_index INTEGER NOT NULL,
                    synonym TEXT NOT NULL,
                    CONSTRAINT fk_molecules_synonyms
                        FOREIGN KEY (molecule_id)
                        REFERENCES molecules(id)
                        ON DELETE CASCADE
                )
            ''')
            this_thread_connection.execute('''
                CREATE TABLE IF NOT EXISTS cas_numbers (
                    id INTEGER PRIMARY KEY,
                    molecule_id INTEGER NOT NULL,
                    cas_number_index INTEGER NOT NULL,
                    cas_number TEXT NOT NULL,
                    CONSTRAINT fk_molecules_cas_numbers
                        FOREIGN KEY (molecule_id)
                        REFERENCES molecules(id)
                        ON DELETE CASCADE
                )
            ''')
            this_thread_connection.execute('''
                CREATE INDEX IF NOT EXISTS idx_molecules_service_identifier_mode_identifier
                ON molecules(service, identifier_mode, identifier)
            ''')
            this_thread_connection.execute('''
                CREATE INDEX IF NOT EXISTS idx_service_synonyms ON synonyms (synonym)
            ''')
            this_thread_connection.execute('''
                CREATE INDEX IF NOT EXISTS idx_service_cas_number ON cas_numbers (cas_number)
            ''')

    def save(self, service: str | Sequence[str], identifier_mode: str | Sequence[str], identifier: str | Sequence[str], molecules: Molecule | Sequence[Molecule]) -> None:

        if isinstance(molecules, Molecule) or molecules is None:
            molecules = [molecules]

        for molecule in molecules:
            if molecule:
                if any(['|' in synonym for synonym in molecule.synonyms]):
                    raise ValueError('molecule names i.e. synonyms must not contain pipe symbols: "|"')

        if isinstance(service, str):
            service = [service] * len(molecules)
            identifier_mode = [identifier_mode] * len(molecules)
            identifier = [identifier] * len(molecules)
            
        this_thread_connection = self.get_connection()
        with this_thread_connection:
            # unfortunately it seems, that python sqlite3 does not support executemany while returning
            # the inserted rows. And even if it would be supported, the order of returned ids is not 
            # guaranteed to be the same order of insertion. Therefore we have to do it one by one.
            # https://discuss.python.org/t/sqlite3-executemany-with-returning-clauses/26291
            # It could be circumvented by constructing the insert statement manually, running with execute
            # and then matching the returned ids to the inserted data. Idk what is faster though.
            molecule_ids = []
            for p, m, i, c in zip(service, identifier_mode, identifier, molecules, strict=True):
                
                if c is None:
                    this_data = (p, m, i, None, None)
                else:
                    
                    if molecule.service != p or molecule.mode != m or molecule.identifier != i:
                        raise ValueError('The molecule properties do not match the input values to the save function.')

                    this_data = (
                        molecule.service.strip(),
                        molecule.mode.strip(),
                        i,
                        molecule.SMILES.strip(),
                        str(molecule.additional_information).strip() if molecule.additional_information else None
                    )

                cursor = this_thread_connection.execute('''
                    INSERT INTO molecules (service, identifier_mode, identifier, SMILES, additional_information)
                    VALUES (?, ?, ?, ?, ?)
                    ''', this_data)
                molecule_ids.append(cursor.lastrowid)

            name_rows_to_insert = []
            cas_number_rows_to_insert = []
            for molecule_id, molecule in zip(molecule_ids, molecules, strict=True):
                if molecule:
                    if molecule.synonyms:
                        this_molecule_synonyms = [(molecule_id, synonym_index, synonym.strip()) for synonym_index, synonym in enumerate(molecule.synonyms)]
                        name_rows_to_insert.extend(this_molecule_synonyms)

                    if molecule.CAS:
                        this_molecule_cas_numbers = [(molecule_id, cas_number_index, cas_number.strip()) for cas_number_index, cas_number in enumerate(molecule.CAS)]
                        cas_number_rows_to_insert.extend(this_molecule_cas_numbers)

            this_thread_connection.executemany('''
                INSERT INTO synonyms (molecule_id, synonym_index, synonym)
                VALUES (?, ?, ?)
            ''', name_rows_to_insert)
            
            this_thread_connection.executemany('''
                INSERT INTO cas_numbers (molecule_id, cas_number_index, cas_number)
                VALUES (?, ?, ?)
            ''', cas_number_rows_to_insert)

    def _search(self, service: str | Sequence[str], identifier_mode: str | Sequence[str], identifier: str | Sequence[str], only_check_for_existance:bool = False) -> Optional[list[Molecule]] | list[Optional[list[Molecule]]] | bool | list[bool]:
        
        if not isinstance(identifier, str):
            search_mode = 'multiple'
            if not (isinstance(identifier_mode, str) and isinstance(service, str)):
                if isinstance(service, str) or \
                    isinstance(identifier_mode, str) or \
                    len(service) != len(identifier_mode) or \
                    len(identifier_mode) != len(identifier):
                    raise ValueError('When searching for multiple molecules, service, mode and identifier, all must be provided as str or same sized lists.')
        else:
            search_mode = 'single'
            if not (isinstance(identifier_mode, str) and isinstance(identifier, str)):
                raise ValueError('When searching for a single molecule, service, mode and identifier, all must be provided as string.')
            
        def rows_to_molecules(service_, identifier_mode_, identifier_, rows):
            molecules = []
            for row in rows:
                SMILES, additional_information, temp_synonyms, temp_cas_numbers = row
                if SMILES:
                    synonyms = []
                    cas_numbers = []

                    # this workaround is needed as GROUP_CONCAT does not preserve order of the values
                    if temp_synonyms:
                        temp_synonyms = {k:v for k,v in (kv.split('|') for kv in temp_synonyms.split('||'))}
                        synonyms = [temp_synonyms[k] for k in sorted(temp_synonyms.keys())]
                    if temp_cas_numbers:
                        temp_cas_numbers = {k:v for k,v in (kv.split('|') for kv in temp_cas_numbers.split('||'))}
                        cas_numbers = [temp_cas_numbers[k] for k in sorted(temp_cas_numbers.keys())]
                    
                    molecules.append(Molecule(SMILES, synonyms, cas_numbers, additional_information, identifier_mode_, service_, 1, identifier_))
            return molecules
            
        this_thread_connection = self.get_connection()
        with this_thread_connection:

            if search_mode == 'single':
                identifier_clause = 'identifier = ?'
                identifier_mode_clause = 'identifier_mode = ? AND'
                values = (service, identifier_mode, identifier)

                if identifier_mode == 'name':
                    identifier_clause = '(identifier = ? COLLATE NOCASE OR synonyms.synonym = ? COLLATE NOCASE)'
                    identifier_mode_clause = ''
                    values = (service, identifier, identifier)
                if identifier_mode == 'cas':
                    identifier_clause = 'cas_numbers.cas_number = ?'
                    identifier_mode_clause = ''
                    values = (service, identifier)

                sql = f'''
                SELECT molecules.id,
                    SMILES,
                    additional_information,
                    GROUP_CONCAT(synonym_index || '|' || synonym, '||'),
                    GROUP_CONCAT( cas_number_index || '|' || cas_number, '||')
                FROM molecules
                    LEFT JOIN synonyms ON molecules.id = synonyms.molecule_id
                    LEFT JOIN cas_numbers ON molecules.id = cas_numbers.molecule_id
                WHERE service = ? AND  {identifier_mode_clause} {identifier_clause}
                GROUP BY molecules.id
                '''
                cursor = this_thread_connection.execute(sql, values)

                molecule_rows = [row[1:] for row in cursor if row[0]]

                if only_check_for_existance:
                    return len(molecule_rows) != 0

                if not molecule_rows:
                    return None
                
                return rows_to_molecules(service, identifier_mode, identifier, molecule_rows)
            
            else:
                this_transaction_unique_temp_table_name = f'tmp_{uuid.uuid4().hex}'

                this_thread_connection.execute(f'''
                    CREATE TABLE {this_transaction_unique_temp_table_name} (
                        search_index INTEGER NOT NULL,
                        service TEXT NOT NULL,
                        identifier_mode TEXT NOT NULL,
                        identifier TEXT NOT NULL
                    )
                ''')

                this_thread_connection.executemany(f'''
                    INSERT INTO {this_transaction_unique_temp_table_name} (search_index, service, identifier_mode, identifier)
                    VALUES (?, ?, ?, ?)
                ''', list(zip(range(len(service)), service, identifier_mode, identifier, strict=True)))

                if only_check_for_existance:
                    optional_columns = ''
                else:
                    optional_columns = ''',
                            SMILES,
                            additional_information,
                            GROUP_CONCAT(synonym_index || '|' || synonym, '||'),
                            GROUP_CONCAT(cas_number_index || '|' || cas_number, '||')
                    '''
                cursor = this_thread_connection.execute(f'''
                    SELECT search_index,
                        molecules.id{optional_columns}
                    FROM molecules
                        LEFT JOIN synonyms ON molecules.id = synonyms.molecule_id
                        LEFT JOIN cas_numbers ON molecules.id = cas_numbers.molecule_id
                        RIGHT JOIN {this_transaction_unique_temp_table_name} AS temp_table ON 
                            molecules.service = temp_table.service AND
                            (
                                (temp_table.identifier_mode = 'cas' and cas_numbers.cas_number  = temp_table.identifier) OR
                                (temp_table.identifier_mode = 'name' and (synonyms.synonym = temp_table.identifier COLLATE NOCASE or molecules.identifier = temp_table.identifier COLLATE NOCASE)) OR
                                (temp_table.identifier_mode != 'cas' and temp_table.identifier_mode != 'name' and molecules.identifier = temp_table.identifier)
                            )
                    GROUP BY search_index, molecules.id
                    ORDER BY search_index
                ''')
                results = [None] * len(service)
                rows = cursor.fetchall()
                if only_check_for_existance:
                    for row in rows:
                        search_index, molecule_id = row
                        results[search_index] = molecule_id is not None
                else:
                    rows_by_search_index = {}
                    for row in rows:
                        search_index, molecule_id, SMILES, additional_information, temp_synonyms, temp_cas_numbers = row

                        entry_found =  molecule_id is not None
                        if entry_found:
                            if search_index not in rows_by_search_index:
                                rows_by_search_index[search_index] = []

                            if SMILES:
                                rows_by_search_index[search_index].append((SMILES, additional_information, temp_synonyms, temp_cas_numbers))

                    for search_index, rows in rows_by_search_index.items():
                        results[search_index] = rows_to_molecules(service[search_index], identifier_mode[search_index], identifier[search_index], rows)

                return results

    def exists(self, service: str | Sequence[str], identifier_mode: str | Sequence[str], identifier: str | Sequence[str]) -> bool | list[bool]:
        return self._search(service, identifier_mode, identifier, only_check_for_existance=True)

    def search(self, service: str | Sequence[str], identifier_mode: str | Sequence[str], identifier: str | Sequence[str]) -> Optional[list[Molecule]] | list[Optional[list[Molecule]]]:
        return self._search(service, identifier_mode, identifier)

    def delete_expired(self) -> None:
        if self.expiration_datetime:
            this_thread_connection = self.get_connection()
            with this_thread_connection:
                if self.expiration_datetime:
                    this_thread_connection.execute('''
                        DELETE FROM molecules
                        WHERE datetime_added < ?
                    ''', (self.expiration_datetime,)) 

    def count(self, service: str = None):
        this_thread_connection = self.get_connection()
        with this_thread_connection:
            if service:
                cursor = this_thread_connection.execute('''
                    SELECT COUNT(*)
                    FROM molecules
                    WHERE service = ?
                ''', (service,))
            else:
                cursor = this_thread_connection.execute('''
                    SELECT COUNT(*)
                    FROM molecules
                ''')

            return cursor.fetchone()[0]
