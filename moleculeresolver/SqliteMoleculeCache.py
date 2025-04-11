from datetime import datetime
import os
import sqlite3
from typing import Optional, Union
import threading
import uuid

from moleculeresolver.molecule import Molecule


class SqliteMoleculeCache:
    """
    A class for caching molecule information using SQLite.

    This class provides methods to initialize, manage, and query a SQLite database
    for storing molecule information. It supports multi-threading and implements
    context management for proper resource handling.

    Attributes:
        db_path (str): Path to the SQLite database file. Defaults to ":memory:".

        expiration_datetime (Optional[datetime]): Expiration date for cached entries.

        _connections (dict): Thread-specific database connections.

        _main_thread_id (int): ID of the main thread.
    """

    def __init__(
        self, db_path: Optional[str] = ":memory:", expiration_datetime: Optional[datetime] = None
    ):
        """
        Initialize a new SqliteMoleculeCache instance.

        Args:
            db_path (Optional[str]): Path to the SQLite database file. Defaults to ":memory:".

            expiration_datetime (Optional[datetime]): Expiration date for cached entries.
        """
        self.db_path = db_path
        self.expiration_datetime = expiration_datetime
        self._connections = {}
        self._main_thread_id = threading.get_ident()

    def __enter__(self) -> "SqliteMoleculeCache":
        """
        Enter the runtime context related to this object.

        Creates tables and deletes expired entries.

        Returns:
            SqliteMoleculeCache: The instance of the class.
        """
        self._create_tables()
        self.delete_expired()
        return self

    def close_child_connections(self) -> None:
        """
        Close all child thread database connections.
        """
        for thread_id, thread_connection in self._connections.items():
            if thread_id != self._main_thread_id:
                if thread_connection:
                    thread_connection.close()
                    self._connections[thread_id] = None

    def __exit__(self, exception_type, exception_value, exception_traceback) -> None:
        """
        Exit the runtime context and close all database connections.

        Closes all child thread connections and optimizes the main thread's connection before closing.
        """
        self.close_child_connections()

        # Close the connection from the main thread
        this_thread_id = threading.get_ident()
        if this_thread_id == self._main_thread_id:
            if self._main_thread_id in self._connections:
                main_thread_connection = self._connections[self._main_thread_id]
                main_thread_connection.execute("PRAGMA analysis_limit=8192")
                main_thread_connection.execute("PRAGMA optimize")
                main_thread_connection.close()
                self._connections.clear()

    def get_connection(self) -> sqlite3.Connection:
        """
        Get or create a thread-specific database connection.

        Returns:
            sqlite3.Connection: A SQLite database connection for the current thread.
        """
        thread_id = threading.get_ident()
        if thread_id not in self._connections:
            self._connections[thread_id] = sqlite3.connect(
                self.db_path, check_same_thread=False
            )
            self._connections[thread_id].execute("PRAGMA foreign_keys = 1")
            self._connections[thread_id].execute("PRAGMA journal_mode=WAL")
            self._connections[thread_id].execute("PRAGMA synchronous=NORMAL")
            self._connections[thread_id].execute("PRAGMA temp_store=MEMORY")

        return self._connections[thread_id]

    def _create_tables(self) -> None:
        """
        Create the necessary tables in the SQLite database if they don't exist.
        """
        this_thread_connection = self.get_connection()
        with this_thread_connection:
            this_thread_connection.execute(
                """
                CREATE TABLE IF NOT EXISTS molecules (
                    id INTEGER PRIMARY KEY,
                    service TEXT NOT NULL,
                    identifier_mode TEXT NOT NULL,
                    identifier TEXT NOT NULL,
                    SMILES TEXT,
                    additional_information TEXT,
                    datetime_added DATETIME DEFAULT CURRENT_TIMESTAMP
                )
            """
            )
            this_thread_connection.execute(
                """
                CREATE TABLE IF NOT EXISTS synonyms (
                    id INTEGER PRIMARY KEY,
                    molecule_id INTEGER NOT NULL,
                    synonym_index INTEGER NOT NULL,
                    synonym TEXT NOT NULL COLLATE NOCASE,
                    CONSTRAINT fk_molecules_synonyms
                        FOREIGN KEY (molecule_id)
                        REFERENCES molecules(id)
                        ON DELETE CASCADE
                )
            """
            )
            this_thread_connection.execute(
                """
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
            """
            )
            this_thread_connection.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_molecules_service_identifier_mode_identifier
                ON molecules(service, identifier_mode, identifier)
            """
            )
            this_thread_connection.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_covering_synonyms ON synonyms (molecule_id, synonym COLLATE NOCASE, synonym_index)
            """
            )
            this_thread_connection.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_covering_cas_number ON cas_numbers (molecule_id, cas_number, cas_number_index)
            """
            )

    def save(
        self,
        service: Union[str, list[str]],
        identifier_mode: Union[str, list[str]],
        identifier: Union[str, list[str]],
        molecules: Union[Molecule, list[Molecule]],
    ) -> None:
        """
        Save molecule information to the database.

        Saves one or multiple Molecule objects to the database, along with their associated service, identifier_mode, and identifier.

        Args:
            service (Union[str, list[str]]): The service(s) associated with the molecule(s).

            identifier_mode (Union[str, list[str]]): The identifier mode(s) for the molecule(s).

            identifier (Union[str, list[str]]): The identifier(s) for the molecule(s).

            molecules (Union[Molecule, list[Molecule]]): The molecule(s) to be saved.

        Raises:
            ValueError: If a molecule's synonyms contain a pipe symbol or if molecule properties don't match the input values.
        """
        if isinstance(molecules, Molecule) or molecules is None:
            molecules = [molecules]

        for molecule in molecules:
            if molecule:
                if any(["|" in synonym for synonym in molecule.synonyms]):
                    raise ValueError(
                        'molecule names i.e. synonyms must not contain pipe symbols: "|"'
                    )

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
            for s, m, i, molecule in zip(
                service, identifier_mode, identifier, molecules
            ):
                if molecule is None:
                    this_data = (s, m, i, None, None)
                else:
                    if (
                        molecule.service != s
                        or molecule.mode != m
                        or molecule.identifier != i
                    ):
                        raise ValueError(
                            "The molecule properties do not match the input values to the save function."
                        )

                    this_data = (
                        molecule.service.strip(),
                        molecule.mode.strip(),
                        i,
                        molecule.SMILES.strip(),
                        (
                            str(molecule.additional_information).strip()
                            if molecule.additional_information
                            else None
                        ),
                    )

                cursor = this_thread_connection.execute(
                    """
                    INSERT INTO molecules (service, identifier_mode, identifier, SMILES, additional_information)
                    VALUES (?, ?, ?, ?, ?)
                    """,
                    this_data,
                )
                molecule_ids.append(cursor.lastrowid)

            name_rows_to_insert = []
            cas_number_rows_to_insert = []
            for molecule_id, molecule in zip(molecule_ids, molecules):
                if molecule:
                    if molecule.synonyms:
                        this_molecule_synonyms = [
                            (molecule_id, synonym_index, synonym.strip())
                            for synonym_index, synonym in enumerate(molecule.synonyms)
                        ]
                        name_rows_to_insert.extend(this_molecule_synonyms)

                    if molecule.CAS:
                        this_molecule_cas_numbers = [
                            (molecule_id, cas_number_index, cas_number.strip())
                            for cas_number_index, cas_number in enumerate(molecule.CAS)
                        ]
                        cas_number_rows_to_insert.extend(this_molecule_cas_numbers)

            this_thread_connection.executemany(
                """
                INSERT INTO synonyms (molecule_id, synonym_index, synonym)
                VALUES (?, ?, ?)
            """,
                name_rows_to_insert,
            )

            this_thread_connection.executemany(
                """
                INSERT INTO cas_numbers (molecule_id, cas_number_index, cas_number)
                VALUES (?, ?, ?)
            """,
                cas_number_rows_to_insert,
            )

    def _search(
        self,
        service: Union[str, list[str]],
        identifier_mode: Union[str, list[str]],
        identifier: Union[str, list[str]],
        only_check_for_existence: Optional[bool] = False,
    ) -> Union[
        Optional[list[Molecule]], list[Optional[list[Molecule]]], bool, list[bool]
    ]:
        """
        Search for molecules in the database based on the provided criteria.

        Supports single and multiple molecule searches. It can either return the full molecule information or just check for existence.

        Args:
            service (Union[str, list[str]]): The service(s) to search in.

            identifier_mode (Union[str, list[str]]): The mode(s) of identification (e.g., 'name', 'cas').

            identifier (Union[str, list[str]]): The identifier(s) to search for.

            only_check_for_existence (Optional[bool]): If True, only check if the molecule exists. Defaults to False.

        Returns:

            Union[Optional[list[Molecule]], list[Optional[list[Molecule]]], bool, list[bool]]:
            - If searching for a single molecule:
                - If only_check_for_existence is False: returns Optional[list[Molecule]]
                - If only_check_for_existence is True: returns bool
            - If searching for multiple molecules:
                - If only_check_for_existence is False: returns list[Optional[list[Molecule]]]
                - If only_check_for_existence is True: returns list[bool]

        Raises:
            ValueError: If the input parameters are inconsistent or invalid for multiple searches.
        """
        if not isinstance(identifier, str):
            search_mode = "multiple"
            if not (isinstance(identifier_mode, str) and isinstance(service, str)):
                if (
                    isinstance(service, str)
                    or isinstance(identifier_mode, str)
                    or len(service) != len(identifier_mode)
                    or len(identifier_mode) != len(identifier)
                ):
                    raise ValueError(
                        "When searching for multiple molecules, service, mode and identifier all must be provided as str or same sized lists."
                    )
        else:
            search_mode = "single"
            if not (isinstance(identifier_mode, str) and isinstance(service, str)):
                raise ValueError(
                    "When searching for a single molecule, service, mode and identifier all must be provided as strings."
                )

        def rows_to_molecules(service_, identifier_mode_, identifier_, rows):
            molecules = []
            for row in rows:
                SMILES, additional_information, temp_synonyms, temp_cas_numbers = row
                if SMILES:
                    synonyms = []
                    cas_numbers = []

                    # Workaround as GROUP_CONCAT does not preserve order of the values
                    if temp_synonyms:
                        temp_synonyms = {
                            int(k): v
                            for k, v in (
                                kv.split("|") for kv in temp_synonyms.split("||")
                            )
                        }
                        synonyms = [
                            temp_synonyms[k] for k in sorted(temp_synonyms.keys())
                        ]
                    if temp_cas_numbers:
                        temp_cas_numbers = {
                            int(k): v
                            for k, v in (
                                kv.split("|") for kv in temp_cas_numbers.split("||")
                            )
                        }
                        cas_numbers = [
                            temp_cas_numbers[k] for k in sorted(temp_cas_numbers.keys())
                        ]

                    molecules.append(
                        Molecule(
                            SMILES,
                            synonyms,
                            cas_numbers,
                            additional_information if additional_information else "",
                            identifier_mode_,
                            service_,
                            1,
                            identifier_,
                        )
                    )
            return molecules

        this_thread_connection = self.get_connection()
        with this_thread_connection:
            if search_mode == "single":
                identifier_clause = "identifier = ?"
                identifier_mode_clause = "identifier_mode = ? AND"
                values = (service, identifier_mode, identifier)

                if identifier_mode == "name":
                    identifier_clause = "identifier = ? COLLATE NOCASE"
                if identifier_mode == "cas":
                    identifier_clause = (
                        "identifier = ? "  # "cas_numbers.cas_number = ?"
                    )

                sql = f"""
                SELECT molecules.id,
                    SMILES,
                    additional_information,
                    GROUP_CONCAT(synonym_index || '|' || synonym, '||'),
                    GROUP_CONCAT(cas_number_index || '|' || cas_number, '||')
                FROM molecules
                    LEFT JOIN synonyms ON molecules.id = synonyms.molecule_id
                    LEFT JOIN cas_numbers ON molecules.id = cas_numbers.molecule_id
                WHERE service = ? AND  {identifier_mode_clause} {identifier_clause}
                GROUP BY molecules.id
                """
                cursor = this_thread_connection.execute(sql, values)

                molecule_rows = [row[1:] for row in cursor if row[0]]

                if only_check_for_existence:
                    return len(molecule_rows) != 0

                if not molecule_rows:
                    return None

                return rows_to_molecules(
                    service, identifier_mode, identifier, molecule_rows
                )

            else:
                this_transaction_unique_temp_table_name = f"tmp_{uuid.uuid4().hex}"

                this_thread_connection.execute(
                    f"""
                    CREATE TEMPORARY TABLE {this_transaction_unique_temp_table_name} (
                        search_index INTEGER NOT NULL,
                        service TEXT NOT NULL,
                        identifier TEXT NOT NULL
                    )
                """
                )

                this_thread_connection.executemany(
                    f"""
                    INSERT INTO {this_transaction_unique_temp_table_name} (search_index, service, identifier)
                    VALUES (?, ?, ?)
                """,
                    list(
                        zip(
                            range(len(service)),
                            service,
                            identifier,
                        )
                    ),
                )

                if only_check_for_existence:
                    optional_columns = ""
                else:
                    optional_columns = """,
                            SMILES,
                            additional_information,
                            GROUP_CONCAT(synonym_index || '|' || synonym, '||'),
                            GROUP_CONCAT(cas_number_index || '|' || cas_number, '||')
                    """

                # Distinction makes queries run much faster
                all_one_service = len(set(service)) == 1
                molecule_join_on_service = "t.service"
                if all_one_service:
                    molecule_join_on_service = f"'{service[0]}'"

                all_one_identifier_mode = len(set(identifier_mode)) == 1
                if not all_one_identifier_mode:
                    raise ValueError(
                        "This class expects all identifier modes to be the same."
                    )

                collation = ""
                if identifier_mode[0] == "name":
                    collation = "COLLATE NOCASE"

                cursor = this_thread_connection.execute(
                    f"""
                    SELECT search_index,
                        m.id{optional_columns}
                        FROM {this_transaction_unique_temp_table_name} AS t
                        INNER JOIN molecules AS m
                            ON m.identifier_mode = '{identifier_mode[0]}'
                            AND m.service = {molecule_join_on_service}
                        LEFT JOIN synonyms AS s
                            ON m.id = s.molecule_id
                        LEFT JOIN cas_numbers AS c
                            ON m.id = c.molecule_id
                    WHERE m.identifier = t.identifier {collation}
                    GROUP BY search_index, m.id
                """
                )
                # TODO: search also the synonyms and cas_numbers tables
                results = [None] * len(service)
                rows = cursor.fetchall()
                if only_check_for_existence:
                    for row in rows:
                        search_index, molecule_id = row
                        results[search_index] = molecule_id is not None
                else:
                    rows_by_search_index = {}
                    for row in rows:
                        (
                            search_index,
                            molecule_id,
                            SMILES,
                            additional_information,
                            temp_synonyms,
                            temp_cas_numbers,
                        ) = row

                        entry_found = molecule_id is not None
                        if entry_found:
                            if search_index not in rows_by_search_index:
                                rows_by_search_index[search_index] = []

                            if SMILES:
                                rows_by_search_index[search_index].append(
                                    (
                                        SMILES,
                                        additional_information,
                                        temp_synonyms,
                                        temp_cas_numbers,
                                    )
                                )

                    for search_index, rows in rows_by_search_index.items():
                        results[search_index] = rows_to_molecules(
                            service[search_index],
                            identifier_mode[search_index],
                            identifier[search_index],
                            rows,
                        )

                return results

    def exists(
        self,
        service: Union[str, list[str]],
        identifier_mode: Union[str, list[str]],
        identifier: Union[str, list[str]],
    ) -> Union[bool, list[bool]]:
        """
        Check if molecule(s) exist in the database based on the provided criteria.

        Supports both single and multiple molecule existence checks.

        Args:
            service (Union[str, list[str]]): The service(s) to search in.
            Can be a single string or a sequence of strings for multiple checks.

            identifier_mode (Union[str, list[str]]): The mode(s) of identification (e.g., 'name', 'cas').
            Can be a single string or a sequence of strings for multiple checks.

            identifier (Union[str, list[str]]): The identifier(s) to search for.
            Can be a single string or a sequence of strings for multiple checks.

        Returns:

            Union[bool, list[bool]]:
                
                - For a single check: A boolean indicating whether the molecule exists.
                - For multiple checks: A list of booleans, each indicating whether the corresponding molecule exists.

        Note:
            This method uses the internal _search method with the 'only_check_for_existence' flag set to True.
        """
        return self._search(
            service, identifier_mode, identifier, only_check_for_existence=True
        )

    def search(
        self,
        service: Union[str, list[str]],
        identifier_mode: Union[str, list[str]],
        identifier: Union[str, list[str]],
    ) -> Union[Optional[list[Molecule]], list[Optional[list[Molecule]]]]:
        """
        Search for molecules based on the given parameters.

        Searches for molecules using the specified service, identifier mode, and identifier.
        Supports both single and multiple searches.

        Args:
            service (Union[str, list[str]]): The service(s) to use for the search.
            Can be a single string or a sequence of strings.
            
            identifier_mode (Union[str, list[str]]): The identifier mode(s) to use.
            Can be a single string or a sequence of strings.

            identifier (Union[str, list[str]]): The identifier(s) to search for.
            Can be a single string or a sequence of strings.

        Returns:

            Union[Optional[list[Molecule]], list[Optional[list[Molecule]]]]:

                - If a single search is performed, returns either None or a list of Molecule objects.
                - If multiple searches are performed, returns a list of results, where each result
                  is either None or a list of Molecule objects.

        Note:
            This method internally calls the _search method to perform the actual search operation.
        """
        return self._search(service, identifier_mode, identifier)

    def delete_expired(self) -> None:
        """
        Delete expired molecules from the cache.

        Removes all molecules from the database that were added before the expiration datetime, if set.

        Note:
            This method only performs the deletion if 'self.expiration_datetime' is set.
        """
        if self.expiration_datetime:
            this_thread_connection = self.get_connection()
            with this_thread_connection:
                this_thread_connection.execute(
                    """
                    DELETE FROM molecules
                    WHERE datetime_added < ?
                """,
                    (self.expiration_datetime,),
                )

    def delete_by_service(self, service: str, mode: Optional[str] = '%') -> None:
        """
        Delete all molecules associated with a specific service from the cache.

        Args:
            service (str): The name of the service whose molecules should be deleted.
        """
        this_thread_connection = self.get_connection()
        with this_thread_connection:
            sql = """
                DELETE FROM molecules
                WHERE service = ? AND identifier_mode LIKE ?
            """
            this_thread_connection.execute(
                sql,
                (service, mode),
            )

    def recreate_all_tables(self) -> None:
        """
        Recreate all tables in the database.

        Closes any existing connections, deletes the database files,
        and then recreates the tables. Use with caution, as it will
        result in data loss.

        Raises:
            RuntimeError: If called in a multi-threaded environment (more than one connection).
        """
        if len(self._connections) > 1:
            raise RuntimeError(
                "Cannot delete cache files in a multi-threaded environment."
            )
        else:
            if len(self._connections) == 1:
                this_thread_connection = self.get_connection()
                this_thread_connection.close()
                self._connections.clear()

            files = [self.db_path, f"{self.db_path}-shm", f"{self.db_path}-wal"]
            for file in files:
                if os.path.exists(file):
                    os.remove(file)

        self._create_tables()

    def count(self, service: Optional[str] = None) -> int:
        """
        Count the number of molecules in the database, optionally filtered by service.

        Args:
            service (Optional[str]): The service to filter by. If None, counts all molecules.

        Returns:
            int: The number of molecules matching the criteria.
        """
        this_thread_connection = self.get_connection()
        with this_thread_connection:
            if service:
                cursor = this_thread_connection.execute(
                    """
                    SELECT COUNT(*)
                    FROM molecules
                    WHERE service = ?
                """,
                    (service,),
                )
            else:
                cursor = this_thread_connection.execute(
                    """
                    SELECT COUNT(*)
                    FROM molecules
                """
                )

            return cursor.fetchone()[0]
