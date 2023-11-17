import collections
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing, contextmanager
import copy
from datetime import datetime
from functools import cache
import gzip
import html
import json
import os
from PIL import ImageFont, ImageDraw
import platform
import requests
import subprocess
import tempfile
import time
from typing import Any, Optional, Sequence, Tuple, Union
import traceback
import unicodedata
import urllib
import uuid
import warnings
import ssl

import openpyxl
from prompt_toolkit import PromptSession
from prompt_toolkit.shortcuts import yes_no_dialog
import regex
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.MolStandardize import canonicalize_tautomer_smiles
from rdkit.Chem import RegistrationHash
from rdkit.Chem.rdchem import ResonanceMolSupplierCallback
from tqdm import tqdm
import urllib3
import xmltodict
from moleculeresolver.rdkitmods import disabling_rdkit_logger
from moleculeresolver.molecule import Molecule
from moleculeresolver.SqliteMoleculeCache import SqliteMoleculeCache


class EmptyResonanceMolSupplierCallback(ResonanceMolSupplierCallback):
    """Workaround for https://github.com/rdkit/rdkit/issues/6704"""

    def __call__(self):
        pass


class CustomHttpAdapter(requests.adapters.HTTPAdapter):
    """Workaround for SSL error from https://stackoverflow.com/questions/71603314/ssl-error-unsafe-legacy-renegotiation-disabled/71646353#71646353"""

    def __init__(self, ssl_context=None, **kwargs):
        self.ssl_context = ssl_context
        super().__init__(**kwargs)

    def init_poolmanager(self, connections, maxsize, block=False):
        self.poolmanager = urllib3.poolmanager.PoolManager(
            num_pools=connections,
            maxsize=maxsize,
            block=block,
            ssl_context=self.ssl_context,
        )


class MoleculeResolver:
    _supported_modes_by_services = {
        "cas_registry": ["name", "smiles", "inchi", "cas"],
        "chebi": ["name", "cas", "formula", "smiles", "inchi", "inchikey"],
        "chemeo": ["name", "smiles", "inchi", "inchikey"],
        "cir": ["formula", "name", "cas", "smiles", "inchi", "inchikey"],
        "comptox": ["name", "cas", "inchikey"],
        "cts": ["name", "cas", "smiles"],
        "nist": ["formula", "name", "cas", "smiles"],
        "opsin": ["name"],
        "pubchem": ["name", "cas", "smiles", "formula", "inchi", "inchikey", "cid"],
        "srs": ["name", "cas"],
    }
    _available_services = sorted(list(_supported_modes_by_services.keys()))
    _supported_modes = []
    for modes in _supported_modes_by_services.values():
        _supported_modes.extend(modes)
    _supported_modes = sorted(list(set(_supported_modes)))

    CAS_regex_with_groups = regex.compile(r"^(\d{2,7})-(\d{2})-(\d)$")
    CAS_regex = r"(\d{2,7}-\d{2}-\d)"
    empirical_formula_regex_compiled = regex.compile(
        r"([A-IK-Z][a-ik-z]*)([0-9]+(?:[.][0-9]+)?)?"
    )
    formula_bracket_group_regex_compiled = regex.compile(
        r"(\((?:[^()]|(?R))*\))(\d+(?:\.\d+)?)"
    )
    non_generic_SMILES_regex_compiled = regex.compile(
        r"^[a-ik-zA-IK-Z0-9%=#$@+\-\[\]\(\)\\\/\:\.]*$"
    )  # non-generic SMILES, no wildcards or unspecified bonds
    InChI_regex_compiled = regex.compile(
        r"^InChI=\dS?\/[0-9a-ik-zA-IK-Z]+\/[0-9a-ik-zA-IK-Z+\-\(\)\\\/,\?]*$"
    )
    InChICode_regex_compiled = regex.compile(r"^[A-Z]{14}\-[A-Z]{8}[SN][A-Z]\-[A-Z]$")
    chemeo_API_token_regex_compiled = regex.compile(r"[a-zA-Z0-9_]+")
    comptox_API_token_regex_compiled = regex.compile(r"[a-z0-9\-]+")
    html_tag_regex_compiled = regex.compile(r"<.*?>")

    @staticmethod
    def chunker(seq: list, size: int) -> set:
        return (seq[pos : pos + size] for pos in range(0, len(seq), size))

    def __init__(
        self,
        available_service_API_keys: Optional[dict[str, Optional[str]]] = None,
        molecule_cache_db_path: Optional[str] = None,
        molecule_cache_expiration: Optional[datetime] = None,
        differentiate_isotopes: bool = False,
    ):
        if not available_service_API_keys:
            available_service_API_keys = {}

        if "chemeo" not in available_service_API_keys:
            available_service_API_keys["chemeo"] = None

        if "comptox" not in available_service_API_keys:
            available_service_API_keys["comptox"] = None

        if not molecule_cache_db_path:
            module_path = os.path.dirname(__file__)
            molecule_cache_db_path = os.path.join(module_path, "molecule_cache.db")

        self.molecule_cache_db_path = molecule_cache_db_path
        self.molecule_cache_expiration = molecule_cache_expiration

        self.available_service_API_keys = available_service_API_keys
        self._differentiate_isotopes = differentiate_isotopes
        self._available_services_with_batch_capabilities = ["srs", "comptox", "pubchem"]
        self._message_slugs_shown = []
        self._session = None
        self._session_CompTox = None
        self._java_path = MoleculeResolver.get_java_path()
        if self._java_path:
            self._available_services_with_batch_capabilities.insert(0, "opsin")
        self._OPSIN_tempfolder = None

        self._init_session()

    def __enter__(self):
        self._disabling_rdkit_logger = disabling_rdkit_logger()
        self._disabling_rdkit_logger.__enter__()
        self.molecule_cache = SqliteMoleculeCache(
            self.molecule_cache_db_path, self.molecule_cache_expiration
        )
        self.molecule_cache.__enter__()
        if "opsin" in self._available_services_with_batch_capabilities:
            self._OPSIN_tempfolder = tempfile.TemporaryDirectory(
                prefix="OPSIN_tempfolder_"
            )
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):

        error_ocurred = (exc_type is not None or exc_value is not None or exc_traceback is not None)

        self._disabling_rdkit_logger.__exit__()
        self.molecule_cache.__exit__()
        self._disabling_rdkit_logger.__exit__()
        if self._OPSIN_tempfolder and not error_ocurred:
            self._OPSIN_tempfolder.cleanup()

    @contextmanager
    def query_molecule_cache(
        self, service: str, identifier_mode: str, identifier: str
    ) -> Tuple[bool, list[Molecule]]:
        molecules = self.molecule_cache.search(service, identifier_mode, identifier)
        entry_available = molecules is not None

        if entry_available:
            yield entry_available, molecules
        else:
            molecules = []
            yield entry_available, molecules
            for molecule in molecules:
                molecule.identifier = identifier

            should_save = True
            if service == "cts" and "CTS_is_down" in self._message_slugs_shown:
                should_save = False
            elif service == "cir" and "CIR_is_down" in self._message_slugs_shown:
                should_save = False

            if should_save:
                self.molecule_cache.save(
                    service,
                    identifier_mode,
                    identifier,
                    molecules if molecules else None,
                )

        return

    @contextmanager
    def query_molecule_cache_batchmode(
        self,
        service: str,
        identifier_mode: str,
        identifiers: list[str],
        save_not_found: bool = True,
    ) -> Tuple[list[str], list[int], Optional[list[list[Molecule], None]]]:
        results = self.molecule_cache.search(
            [service] * len(identifiers),
            [identifier_mode] * len(identifiers),
            identifiers,
        )

        identifiers_to_search = []
        indices_of_identifiers_to_search = []
        for (molecule_index, molecule), identifier in zip(
            enumerate(results), identifiers
        ):
            if molecule is None:
                identifiers_to_search.append(identifier)
                indices_of_identifiers_to_search.append(molecule_index)

        yield identifiers_to_search, indices_of_identifiers_to_search, results
        identifiers_to_save = []
        molecules_to_save = []
        for molecule_index, identifier in zip(
            indices_of_identifiers_to_search, identifiers_to_search
        ):
            molecules = results[molecule_index]
            if molecules:
                for molecule in molecules:
                    identifiers_to_save.append(identifier)
                    molecules_to_save.append(molecule)
            else:
                if save_not_found:
                    identifiers_to_save.append(identifier)
                    molecules_to_save.append(None)

        self.molecule_cache.save(
            [service] * len(identifiers_to_save),
            [identifier_mode] * len(identifiers_to_save),
            identifiers_to_save,
            molecules_to_save,
        )
        return

    def _init_session(
        self, pool_connections: Optional[int] = None, pool_maxsize: Optional[int] = None
    ) -> None:
        if self._session is not None:
            return

        if pool_connections is None:
            # reserve two connections for each service in case of different hosts are used.
            pool_connections = len(MoleculeResolver._available_services) * 2

        if pool_maxsize is None:
            # default from documentation
            pool_maxsize = 10
        else:
            pool_maxsize = max(pool_maxsize, 10)

        self._session = requests.Session()
        self._session.mount(
            "https://",
            requests.adapters.HTTPAdapter(
                pool_connections=pool_connections,
                pool_maxsize=pool_maxsize,
                max_retries=2,
            ),
        )

        self._session_CompTox = requests.Session()
        ctx = ssl.create_default_context(ssl.Purpose.SERVER_AUTH)
        ctx.options |= 0x4
        self._session_CompTox.mount(
            "https://",
            CustomHttpAdapter(
                ctx,
                pool_connections=pool_connections,
                pool_maxsize=pool_maxsize,
                max_retries=2,
            ),
        )

    def _resilient_request(
        self,
        url: str,
        kwargs: Optional[dict[str, Any]] = None,
        request_type: str = "get",
        accepted_status_codes: list[int] = [200],
        rejected_status_codes: list[int] = [404],
        max_retries: int = 10,
        sleep_time: Union[int, float] = 2,
        allow_redirects: bool = False,
        json: Any = None,
        return_response: bool = False,
    ) -> Optional[requests.Response]:
        if not kwargs:
            kwargs = {}

        if request_type not in ["get", "post"]:
            raise ValueError("The request_type must be either 'get' or 'post'.")

        if "timeout" not in kwargs:
            kwargs["timeout"] = 5

        headers = {}
        user_agent_is_set = False
        if "headers" in kwargs:
            headers = kwargs["headers"]
            user_agent_is_set = "user-agent" in [key.lower() for key in headers.keys()]

        if user_agent_is_set is False:
            headers[
                "user-agent"
            ] = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/116.0.0.0 Safari/537.36"

        kwargs["headers"] = headers

        session = (
            self._session if "comptox" not in url.lower() else self._session_CompTox
        )
        n_try = 0
        while True:
            n_try += 1
            try:
                response = None

                if request_type == "get":
                    response = session.get(
                        url, **kwargs, json=json, allow_redirects=allow_redirects
                    )
                if request_type == "post":
                    response = session.post(
                        url, **kwargs, json=json, allow_redirects=allow_redirects
                    )

                if response.status_code in accepted_status_codes:
                    if return_response:
                        return response
                    response_text = response.content.decode("utf8", errors="ignore")
                    response_text = response_text.replace("\u200b", "")
                    # find encoding errors
                    if response_text != response.text:
                        if response_text.count(r"\u") + response_text.count(r"\x") > 0:
                            raise UnicodeEncodeError(
                                "Wrong charachter encoding was used."
                            )
                    return response_text
                elif response.status_code in rejected_status_codes:
                    return None
                else:
                    raise requests.exceptions.HTTPError("Wrong status_code.")

            except requests.exceptions.RequestException as error:
                time.sleep(sleep_time)

                if n_try > max_retries:
                    if isinstance(error, requests.exceptions.ConnectionError):
                        raise error

                    if response is not None:
                        print(
                            f"ERROR ocurred during request to:\n{url}\n{kwargs}\n HTTP request status code: {response.status_code}\n"
                        )
                    else:
                        print(
                            f"ERROR ocurred during request to:\n{url}\n{kwargs}\n",
                            traceback.format_exc(),
                        )

                    return None

    @staticmethod
    @cache
    def standardize_SMILES(
        SMILES: str,
        standardize: bool,
        /,
        remove_atom_mapping_number: bool = True,
        disconnect_metals: bool = False,
        normalize: bool = True,
        reionize: bool = True,
        uncharge: bool = False,
        stereo: bool = True,
        canonicalize_tautomer: bool = False,
        isomeric_SMILES: bool = True,
    ) -> str:
        mol = MoleculeResolver.get_from_SMILES(SMILES)
        if mol is None:
            return None
        return (
            Chem.MolToSmiles(
                MoleculeResolver.standardize_molecule(
                    mol,
                    remove_atom_mapping_number,
                    disconnect_metals,
                    normalize,
                    reionize,
                    uncharge,
                    stereo,
                    canonicalize_tautomer,
                ),
                isomericSmiles=isomeric_SMILES,
            )
            if standardize
            else SMILES
        )

    @staticmethod
    @cache
    def try_disconnect_more_metals(SMILES, standardize):
        # this does something very similar like the metal disconnector from rdkit
        # but includes more metals and behaves correctly on Hg: https://github.com/rdkit/rdkit/discussions/6729
        # if the problem for Hg is a bug, this should be replaced by the metal disconnector from rdkit
        metals = "[#3,#11,#19,#37,#55,#87,#4,#12,#20,#38,#56,#88,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#13,#31,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#49,#50,#72,#73,#74,#75,#76,#77,#78,#79,#80,#81,#82,#83]"
        SMARTS = f"{metals}~[B,C,#14,P,#33,#51,S,#34,#52,F,Cl,Br,I,#85]"
        mol = Chem.MolFromSmiles(SMILES)
        if not mol:
            return SMILES
        patt = Chem.MolFromSmarts(SMARTS)
        hit_ats = list(mol.GetSubstructMatches(patt))
        if not hit_ats:
            return SMILES

        property_name = "molecule_resolver_charge"
        cation_charges = {}
        bonds_to_be_broken = []
        for cation_id, anion_id in hit_ats:
            cation = mol.GetAtomWithIdx(cation_id)
            if cation_id not in cation_charges:
                cation_charge = cation.GetTotalDegree()
                cation_charges[cation_id] = cation_charge
                cation.SetIntProp(property_name, cation_charge)
            else:
                cation_charge = cation_charges[cation_id]

            bond_to_be_broken = mol.GetBondBetweenAtoms(cation_id, anion_id)
            bonds_to_be_broken.append(bond_to_be_broken.GetIdx())

            anion_charge = -1 * int(bond_to_be_broken.GetBondType())
            anion = mol.GetAtomWithIdx(anion_id)
            anion.SetIntProp(property_name, anion_charge)
            cation_charges[cation_id] += anion_charge

        if not all(v == 0 for v in cation_charges.values()):
            return SMILES  # charge missmatch

        fragment_mols = rdmolops.GetMolFrags(
            rdmolops.FragmentOnBonds(mol, bonds_to_be_broken), asMols=True
        )

        ion_mols = []
        for m in fragment_mols:
            rwmol = Chem.RWMol(m)
            ai_dummy_atoms = []
            for a in rwmol.GetAtoms():
                if a.HasProp(property_name):
                    a.SetFormalCharge(a.GetIntProp(property_name))
                if a.GetSymbol() == "*":
                    ai_dummy_atoms.append(a.GetIdx())

            for ai in sorted(ai_dummy_atoms, reverse=True):
                rwmol.RemoveAtom(ai)

            ion_mols.append(rwmol.GetMol())

        return MoleculeResolver.standardize_SMILES(
            ".".join([Chem.MolToSmiles(m) for m in ion_mols]), standardize
        )

    @staticmethod
    def standardize_molecule(
        mol: Chem.rdchem.Mol,
        /,
        remove_atom_mapping_number: bool = True,
        disconnect_metals: bool = False,
        normalize: bool = True,
        reionize: bool = True,
        uncharge: bool = False,
        stereo: bool = True,
        canonicalize_tautomer: bool = False,
    ) -> Optional[Chem.rdchem.Mol]:
        if mol is None:
            return None

        mol_fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

        if len(mol_fragments) > 1:
            all_mol_fragments = []
            for mol_fragment in mol_fragments:
                standardized_mol_fragment = MoleculeResolver.standardize_molecule(
                    mol_fragment,
                    remove_atom_mapping_number=remove_atom_mapping_number,
                    disconnect_metals=disconnect_metals,
                    normalize=normalize,
                    reionize=reionize,
                    uncharge=uncharge,
                    stereo=stereo,
                    canonicalize_tautomer=canonicalize_tautomer,
                )
                if standardized_mol_fragment is None:
                    return None
                all_mol_fragments.append(standardized_mol_fragment)

            standardize_mol = all_mol_fragments[0]
            for mol_fragment in all_mol_fragments[1:]:
                standardize_mol = Chem.CombineMols(standardize_mol, mol_fragment)

            return standardize_mol

        # based on datamol version on https://github.com/datamol-org/datamol
        mol = copy.deepcopy(mol)

        if disconnect_metals:
            md = rdMolStandardize.MetalDisconnector()
            mol = md.Disconnect(mol)

        if normalize:
            # for some molecules this gives weird results. e.g. DMSO
            if Chem.MolToSmiles(mol) not in ["CS(C)=O"]:
                mol = rdMolStandardize.Normalize(mol)

        if reionize:
            reionizer = rdMolStandardize.Reionizer()
            try:
                mol = reionizer.reionize(mol)
            except:
                warnings.warn(
                    "Reionization step could not be performed. It will be skipped."
                )

        if uncharge:
            uncharger = rdMolStandardize.Uncharger()
            mol = uncharger.uncharge(mol)

        if stereo:
            Chem.AssignStereochemistry(mol, force=False, cleanIt=True)

        if canonicalize_tautomer:
            mol = MoleculeResolver.get_from_SMILES(
                canonicalize_tautomer_smiles(Chem.MolToSmiles(mol))
            )

        if remove_atom_mapping_number:
            [a.SetAtomMapNum(0) for a in mol.GetAtoms()]

        return mol

    @staticmethod
    def molecule_contains_isotope(mol):
        return any([atom.GetIsotope() != 0 for atom in mol.GetAtoms()])

    @staticmethod
    def _check_and_flatten_identifiers_and_modes(
        identifiers: Union[str, list[str], list[list[str]]],
        modes: Union[str, list[str]],
    ) -> tuple[list[str], list[str], list[str], list[str], str]:
        if isinstance(identifiers, str):
            identifiers = [identifiers]
        if isinstance(modes, str):
            modes = [modes]

        flattened_identifiers = []
        flattened_modes = []
        for identifier, mode in zip(identifiers, modes):
            if all([not isinstance(identifier, t) for t in [int, str, list, tuple]]):
                raise TypeError("An identifier can only be an int, str, list or tuple.")

            temp_identifiers = identifier
            if not isinstance(identifier, tuple) and not isinstance(identifier, list):
                temp_identifiers = [identifier]

            for i in temp_identifiers:
                flattened_identifiers.append(str(i).strip())
                flattened_modes.append(mode.strip().lower())

        synonyms = []
        CAS = set()
        given_SMILES = []
        for identifier, mode in zip(flattened_identifiers, flattened_modes):
            if mode == "name":
                synonyms.append(identifier)
            elif mode == "cas":
                if not MoleculeResolver.is_valid_CAS(identifier):
                    raise ValueError("You provided an invalid CAS.")
                CAS.add(identifier)
            elif mode == "inchi":
                if not MoleculeResolver.is_valid_InChI(identifier):
                    raise ValueError("You provided an invalid InChI.")
                given_SMILES.append(MoleculeResolver.InChI_to_SMILES(identifier))
            elif mode == "smiles":
                if not MoleculeResolver.is_valid_SMILES(identifier):
                    raise ValueError("You provided an invalid SMILES.")
                given_SMILES.append(identifier)

        if not len(given_SMILES) < 2:
            raise ValueError(
                "Only one valid structure is accepted to search for (this includes SMILES and InChI)."
            )
        if len(given_SMILES) > 0:
            given_SMILES = given_SMILES[0]
        else:
            given_SMILES = None

        MoleculeResolver._check_parameters(
            modes=flattened_modes,
            identifiers=flattened_identifiers,
            context="find_single",
        )

        return flattened_identifiers, flattened_modes, synonyms, CAS, given_SMILES

    @staticmethod
    def _is_list_of_list_of_str(value: list[list[str]]) -> bool:
        for items in value:
            if not isinstance(items, list):
                return False
            if not all([isinstance(item, str) for item in items]):
                return False
        return True

    @staticmethod
    def _check_parameters(
        *,
        modes=None,
        services=None,
        identifiers=None,
        required_formulas=None,
        required_charges=None,
        required_structure_types=None,
        context="get_molecule",
    ):
        if modes is not None:
            if context == "get_molecule" or context == "get_molecules_batch":
                if not isinstance(modes, str):
                    raise TypeError("The mode parameter can only be a string.")
                if modes not in MoleculeResolver._supported_modes:
                    raise ValueError("The mode parameter can only be a supported mode.")
                if (
                    services is not None
                    and modes
                    not in MoleculeResolver._supported_modes_by_services[services]
                ):
                    raise ValueError(
                        f"The chosen mode ({modes}) is not compatible with the service {services}."
                    )

            if context == "batch" or context == "find_single":
                if not isinstance(modes, list):
                    raise TypeError(
                        "The modes parameter can only be a list of strings."
                    )
                if not all([isinstance(temp_mode, str) for temp_mode in modes]):
                    if not MoleculeResolver._is_list_of_list_of_str(modes):
                        raise TypeError(
                            "The modes parameter can only be a list of strings."
                        )
                if (
                    identifiers is not None
                    and isinstance(identifiers, list)
                    and not len(modes) == len(identifiers)
                ):  # you can give several identifiers e.g. [['sodium chloride', 'NaCl']] with modes = ['name']
                    raise ValueError(
                        "The modes and identifiers parameters must be of the same length."
                    )

            if context == "find_multiple" and not all(
                [len(mode) > 0 for mode in modes]
            ):
                raise ValueError(
                    "The list of modes cannot include an empty string or an empty list."
                )

        if identifiers is not None:
            if not identifiers:
                raise ValueError("The identifiers parameter cannot be empty.")

        if services is not None:
            if context == "get_molecule" or context == "get_molecules_batch":
                if not isinstance(services, str):
                    raise TypeError("The service parameter can only be a string.")
                if services not in MoleculeResolver._available_services:
                    raise ValueError(f"The service {services} is not supported.")

            if context == "batch" or context == "find_single":
                if not isinstance(services, list):
                    raise TypeError(
                        "The services parameter can only be a list of strings."
                    )
                if not all(
                    [isinstance(temp_service, str) for temp_service in services]
                ):
                    raise TypeError(
                        "The services parameter can only be a list of strings."
                    )
                if not all(
                    [
                        temp_service in MoleculeResolver._available_services
                        for temp_service in services
                    ]
                ):
                    raise ValueError(
                        "The list of services can only include supported services."
                    )

        if required_formulas is not None:
            if (
                context == "find_multiple"
                and identifiers is not None
                and not len(identifiers) == len(required_formulas)
            ):
                raise ValueError(
                    "The identifiers and required_formulas parameters must be of the same length."
                )

            if context == "get_molecule" and not isinstance(required_formulas, str):
                raise TypeError("required_formula must be a string.")

        if required_charges is not None:
            if context == "get_molecule":
                if not isinstance(required_charges, str) and not isinstance(
                    required_charges, int
                ):
                    raise TypeError(
                        "required_charge can be zero, not_zero, positive, negative or any integer."
                    )
                if isinstance(required_charges, str):
                    if required_charges not in [
                        "zero",
                        "non_zero",
                        "positive",
                        "negative",
                    ]:
                        raise ValueError(
                            "required_charge can be 'zero', 'non_zero', 'positive', 'negative' or any integer."
                        )

            if context == "find_multiple":
                if not isinstance(required_charges, list):
                    raise TypeError("The required_charges must be a list.")
                if identifiers is not None and not len(identifiers) == len(
                    required_charges
                ):
                    raise ValueError(
                        "The parameters identifiers and required_charges must be of the same length."
                    )

        if required_structure_types is not None:
            if (
                context == "find_multiple"
                and identifiers is not None
                and not len(identifiers) == len(required_structure_types)
            ):
                raise ValueError(
                    "The parameters identifiers and required_structure_types must be of the same length."
                )

            if context == "get_molecule":
                if not isinstance(required_structure_types, str):
                    raise TypeError("required_structure_type must be a string.")
                if required_structure_types not in [
                    "mixture_neutrals",
                    "mixture_ions",
                    "neutral",
                    "salt",
                    "ion",
                    "mixture_neutrals_salts",
                    "mixture_neutrals_ions",
                ]:
                    raise ValueError(
                        "required_structure_type must be one of the following: 'mixture_neutrals','mixture_ions', 'neutral', 'salt', 'ion','mixture_neutrals_salts', 'mixture_neutrals_ions'"
                    )

    @staticmethod
    def take_most_common(
        container: list[Any], number_to_take: Optional[int] = None
    ) -> list[Any]:
        if not number_to_take:
            number_to_take = len(container)

        if len(container) < 2:
            return container

        original_container = container

        is_str = isinstance(container[0], str)
        if is_str:
            original_container = [item.strip() for item in container]
            container = [item.strip().lower() for item in container]

        index_container = []
        items_seen = []
        for item_index, item in enumerate(container):
            if item in items_seen:
                continue
            n = container.count(item)
            for i in range(n):
                index_container.append(item_index)
            items_seen.append(item)

        counter = collections.Counter(index_container)
        sorted_item_indices_to_take = [
            item_index for item_index, _ in counter.most_common(number_to_take)
        ]

        final_container = []
        for item_index in sorted_item_indices_to_take:
            final_container.append(original_container[item_index])

        return final_container

    @staticmethod
    def filter_and_sort_synonyms(
        synonyms: list[str], number_of_synonyms_to_take: int = 5, strict: bool = False
    ) -> list[str]:
        # heuristics used here are mainly to remove synonyms from pubchem
        synonyms = [synonym.strip() for synonym in synonyms if synonym is not None]
        if len(synonyms) == 0:
            return []

        temp = []
        for synonym in synonyms:
            temp.extend(synonym.split("|"))

        synonyms_taken = []
        for synonym in temp:
            if synonym is None:
                continue

            synonym = synonym.strip()
            if len(synonym) == 0:
                continue

            if synonym is not None:
                # typical tokens to be rejected
                tokens = [
                    "SCHEMBL",
                    "AKOS",
                    "MFCD",
                    "AKOS",
                    "CYPHOS",
                    "DTXSID",
                    "UNII",
                    "NSC",
                    "DSSTox",
                    "FLUKA",
                    "ACMC",
                    "%",
                    "(at)",
                    "SIGMA",
                    "ALDRICH",
                    "RIEDEL",
                    "SIAL",
                ]
                n_tokens_found = sum([synonym.lower().count(t.lower()) for t in tokens])
                if n_tokens_found > 0:
                    continue

                if not synonym[-1].isalpha() and synonym[-1] != "-":
                    continue

                # if completely uppercase convert to lower
                if synonym.isupper():
                    synonym = synonym.lower()

                ratio_uppercase_letters = sum(1 for c in synonym if c.isupper()) / len(
                    synonym
                )

                if ratio_uppercase_letters < 0.5:
                    match = regex.match(
                        MoleculeResolver.CAS_regex_with_groups, synonym
                    )  # filter CAS numbers
                    # filter names with 3 consecutive numbers or containing some keywords
                    match2 = regex.search(
                        r"\d{3}|\bin\b|\bgrade\b|\bsolution\b|\bstandard\b|\bstabilized\b|\bdimer\b|\bfor\b|\btablet\b|\btotal\b|\bcode\b|\bNo\.\b|\banhydrous\b",
                        synonym,
                        regex.IGNORECASE,
                    )
                    if not match and not match2:
                        synonym = synonym.replace("<em>", "").replace(
                            "</em>", ""
                        )  # replace(';', ', ')

                        temp_synonym = synonym.split(", ")
                        if len(temp_synonym) > 1:
                            n_parts_ending_in_hyphen = 0
                            for i, part in enumerate(temp_synonym):
                                if i == 0:
                                    continue
                                if part.endswith("-"):
                                    n_parts_ending_in_hyphen += 1
                            if n_parts_ending_in_hyphen == len(temp_synonym) - 1:
                                temp_synonym.reverse()
                                synonym = "".join(temp_synonym)

                        synonyms_taken.append(synonym)

        final_synonyms_taken = MoleculeResolver.take_most_common(
            synonyms_taken, number_of_synonyms_to_take
        )

        if len(final_synonyms_taken) == 0:
            if not strict:
                if len(synonyms) == 0:
                    synonyms.append("Noname")
                final_synonyms_taken = [synonyms[0]]

        return final_synonyms_taken

    @staticmethod
    def expand_name_heuristically(
        name: str,
        prefixes_to_delete: Optional[Sequence[str]] = None,
        suffixes_to_use_as_prefix: Optional[Sequence[str]] = None,
        suffixes_to_delete: Optional[Sequence[str]] = None,
        parts_to_delete: Optional[Sequence[str]] = None,
        maps_to_replace: Optional[dict[str, str]] = None,
    ) -> str:
        name = name.strip()
        original_name = name
        names = [original_name]
        name = html.unescape(name)

        name_without_html_tags = MoleculeResolver.html_tag_regex_compiled.sub("", name)
        if name != name_without_html_tags:
            if "()" not in name_without_html_tags:
                names.append(name_without_html_tags)

        # add
        # acid, alskjdalskjd acid
        # amine,N,N-  amine
        # delete empty spaces in some places 2, 3-dimethil
        # α -> alpha
        # β -> beta   (dl)

        name_parts = name.split(", ")
        if len(name_parts) > 1:
            n_parts_ending_in_hyphen = 0
            for i, part in enumerate(name_parts):
                if i == 0:
                    continue
                if part.endswith("-"):
                    n_parts_ending_in_hyphen += 1
            if n_parts_ending_in_hyphen == len(name_parts) - 1:
                name_parts.reverse()
                names.append("".join(name_parts))

        if prefixes_to_delete is None:
            prefixes_to_delete = [
                "±",
                "\+\,-",
                "\+/-",
                "\+-",
                "\.\+\,-\.",
                "\.\+/-\.",
                "\.\+-\.",
                "(dl)",
                "DL",
                "RS",
                "\(<\+\->\)\-",
                "\(2H\d+\)",
                "\[2H\d+\]",
            ]

        if suffixes_to_delete is None:
            suffixes_to_delete = [
                "mixed isomers",
                "isomers",
                "tautomers",
                "dl and meso",
                "±",
                "\+\,-",
                "\+/-",
                "\+-",
                "\.\+\,-\.",
                "\.\+/-\.",
                "\.\+-\.",
                "cis and trans",
                "-d2",
            ]

        if suffixes_to_use_as_prefix is None:
            suffixes_to_use_as_prefix = ["R", "S"]

        if parts_to_delete is None:
            parts_to_delete = ["±"]

        if maps_to_replace is None:
            maps_to_replace = [
                [r"([a-zA-Z]{2,})-([a-zA-Z]{2,})", r"\1\2"],
                [r"trans(-\d+-[a-z]{3,}ene)\b", r"(E)\1"],
                [r"trans(-\d+-[a-z]{3,}ene)\b", r"E\1"],
                [r"cis(-\d+-[a-z]{3,}ene)\b", r"(Z)\1"],
                [r"cis(-\d+-[a-z]{3,}ene)\b", r"Z\1"],
                [r"(.*)(\((?:Z|E)\)-)(.*)", r"\2\1\3"],
                [r"«|»", r""],
                [r"- -", "-"],
            ]

        new_name = name
        for prefix in prefixes_to_delete:
            m = regex.match(rf"\s*\(?\s*{prefix}\s*\)?\s*-?(.*)$", new_name)
            if m:
                new_name = m.group(1)

        for suffix in suffixes_to_delete:
            m = regex.match(rf"(.*)[;,]+\s*\(?\s*{suffix}\s*\)?\s*-?\s*$", new_name)
            if m:
                new_name = m.group(1)

        if new_name not in names:
            names.append(new_name)

        for suffix in suffixes_to_use_as_prefix:
            m = regex.match(rf"(.*)[;,]+\s+\(?\s*{suffix}\s*\)?\s*-?\s*$", new_name)
            if m:
                new_name = f"({suffix})-{m.group(1)}"
                if new_name not in names:
                    names.append(new_name)

        for part in parts_to_delete:
            for temporary_name in names:
                new_name = regex.sub(
                    rf"([;, ]*\(?\s*{part}\s*\)?\s*-?)", "", temporary_name
                )
                if new_name not in names:
                    names.append(new_name)

        for temporary_name in names:
            for pattern, replacement in maps_to_replace:
                ms = regex.findall(pattern, temporary_name)
                if ms:
                    if all(
                        [
                            v not in ms[0]
                            for v in [
                                "trans",
                                "cis",
                                "alpha",
                                "beta",
                                "tert",
                                "sec",
                                "gamma",
                                "erythro",
                                "threo",
                            ]
                        ]
                    ):
                        new_name = regex.sub(pattern, replacement, temporary_name)
                        if new_name not in names:
                            names.append(new_name)

        if names[0] != original_name:
            return names.insert(0, original_name)
        return names

    @staticmethod
    def combine_molecules(
        molecules: list[Molecule], return_None_when_different_SMILES: bool = False
    ) -> Molecule:
        if not molecules:
            return None

        if len(molecules) == 1:
            return molecules[0]

        all_SMILES_are_equal = all(
            [
                MoleculeResolver.are_equal(
                    MoleculeResolver.get_from_SMILES(cmp.SMILES),
                    MoleculeResolver.get_from_SMILES(molecules[0].SMILES),
                )
                for cmp in molecules
            ]
        )
        all_identifiers_are_equal = all(
            [cmp.identifier == molecules[0].identifier for cmp in molecules]
        )
        all_services_are_equal = all(
            [cmp.service == molecules[0].service for cmp in molecules]
        )

        if not all_SMILES_are_equal or (
            all_services_are_equal and not all_identifiers_are_equal
        ):
            # if not all_SMILES_are_equal and all_services_are_equal and all_identifiers_are_equal:
            #     if return_None_when_different_SMILES:
            #         return None # service does not deliver one single component back for the same identifier
            raise ValueError(
                "Combining molecules is only allowed if the structures are the same and found by the same identifier in the case of the same service."
            )

        all_molecules_used_same_service = all(
            [cmp.service == molecules[0].service for cmp in molecules]
        )

        merged_synonyms = []
        merged_CAS = []
        merged_modes = []
        merged_services = []
        merged_additional_information = []

        for molecule in molecules:
            merged_additional_information.append(molecule.additional_information)
            merged_CAS.extend(molecule.CAS)
            merged_synonyms.extend(molecule.synonyms)
            merged_modes.append(molecule.mode)
            merged_services.append(molecule.service)

        merged_synonyms = MoleculeResolver.filter_and_sort_synonyms(merged_synonyms)

        # prefer CAS from official registry
        if "cas_registry" in merged_services:
            index_cas_registry = merged_services.index("cas_registry")
            merged_CAS = molecules[index_cas_registry].CAS
        else:
            merged_CAS = MoleculeResolver.filter_and_sort_CAS(merged_CAS)

        if all_molecules_used_same_service:
            merged_additional_information = [
                v for v in merged_additional_information if v is not None
            ]
            if merged_additional_information:
                merged_additional_information = str(
                    sorted(merged_additional_information)
                )
            else:
                merged_additional_information = ""
            merged_modes = str(sorted(merged_modes))
            merged_services = merged_services[0]
            number_of_crosschecks = 1
        else:
            merged_modes = "; ".join(merged_modes)
            merged_additional_information = "; ".join(merged_additional_information)
            number_of_crosschecks = len(merged_services)
            merged_services = "; ".join(merged_services)

        return Molecule(
            molecules[0].SMILES,
            merged_synonyms,
            merged_CAS,
            merged_additional_information,
            merged_modes,
            merged_services,
            number_of_crosschecks,
            molecules[0].identifier,
        )

    @staticmethod
    def group_molecules_by_structure(
        molecules: list[Molecule], group_also_by_services: bool = True
    ) -> dict[str, list[Molecule]]:
        if len(molecules) == 1:
            return {molecules[0].SMILES: molecules}

        # if all molecules given are from a primary service, combine them for each service if needed
        all_molecules_from_primary_source = all(
            [
                cmp.service in MoleculeResolver._available_services
                and cmp.number_of_crosschecks == 1
                for cmp in molecules
            ]
        )
        if all_molecules_from_primary_source:
            grouped_molecules = {}
            grouped_molecules_SMILES = {}
            grouped_molecules_molecules = {}
            all_keys = []
            for cmp in molecules:
                cmp_mol = MoleculeResolver.get_from_SMILES(cmp.SMILES)
                equal_SMILES = None
                old_key = None
                for key, grouped_molecule_mol in grouped_molecules_molecules.items():
                    if MoleculeResolver.are_equal(cmp_mol, grouped_molecule_mol):
                        if equal_SMILES is None:
                            equal_SMILES = grouped_molecules_SMILES[key]
                        equal_SMILES = sorted(
                            [equal_SMILES, grouped_molecules_SMILES[key]]
                        )[0]
                        old_key = key
                        break

                if equal_SMILES is None:
                    equal_SMILES = cmp.SMILES

                if group_also_by_services:
                    key = f"{cmp.service}_{equal_SMILES}"
                else:
                    key = equal_SMILES

                if old_key:
                    if key != old_key:
                        grouped_molecules[key] = grouped_molecules[old_key]
                        del grouped_molecules[old_key]

                if key not in grouped_molecules:
                    grouped_molecules[key] = []
                    grouped_molecules_SMILES[key] = equal_SMILES
                    grouped_molecules_molecules[key] = MoleculeResolver.get_from_SMILES(
                        equal_SMILES
                    )

                grouped_molecules[key].append(cmp)
                all_keys.append(key)

            # if all where searched with same service, mode and are a primary source
            # take the most common structure
            if len(grouped_molecules) > 1:
                all_molecules_used_same_service_and_mode = all(
                    [
                        cmp.mode == molecules[0].mode
                        and cmp.service == molecules[0].service
                        for cmp in molecules
                    ]
                )
                if all_molecules_used_same_service_and_mode:
                    smiles_count = collections.Counter(all_keys).most_common()
                    if smiles_count[0][1] > smiles_count[1][1]:
                        key = smiles_count[0][0]
                        return {grouped_molecules_SMILES[key]: grouped_molecules[key]}

            grouped_molecules = {
                grouped_molecules_SMILES[key]: grouped_molecules[key]
                for key in all_keys
            }
        else:
            raise ValueError(
                "This function only takes molecules that were not previously combined from different sources."
            )

        return grouped_molecules

    @staticmethod
    def filter_molecules(
        molecules: list[Molecule],
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = False,
    ) -> list[Molecule]:
        filtered_molecules = []

        if not isinstance(molecules, list):
            molecules = [molecules]
        if len(molecules) == 0:
            return None

        for cmp in molecules:
            if cmp is None:
                continue
            if cmp.SMILES is None:
                continue

            mol = MoleculeResolver.get_from_SMILES(cmp.SMILES)
            if mol is None:
                continue

            # filter isotope containing structures
            if MoleculeResolver.molecule_contains_isotope(mol):
                continue

            if MoleculeResolver.check_SMILES(
                cmp.SMILES, required_formula, required_charge, required_structure_type
            ):
                cmp.SMILES = MoleculeResolver.standardize_SMILES(
                    cmp.SMILES, standardize
                )
                filtered_molecules.append(cmp)
            elif required_structure_type == "salt":
                new_SMILES = MoleculeResolver.try_disconnect_more_metals(
                    cmp.SMILES, standardize
                )
                if MoleculeResolver.check_SMILES(
                    new_SMILES,
                    required_formula,
                    required_charge,
                    required_structure_type,
                ):
                    cmp.SMILES = new_SMILES
                    filtered_molecules.append(cmp)

        return filtered_molecules

    @staticmethod
    def filter_and_combine_molecules(
        molecules: Union[list[Molecule], Molecule],
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = False,
    ) -> Optional[Molecule]:
        filtered_molecules = MoleculeResolver.filter_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

        if not filtered_molecules:
            return None

        grouped_molecules = MoleculeResolver.group_molecules_by_structure(
            filtered_molecules
        )

        final_molecules = []
        for _, cmps_to_combine in grouped_molecules.items():
            final_molecules.append(MoleculeResolver.combine_molecules(cmps_to_combine))

        if len(final_molecules) == 1:
            return final_molecules[0]
        else:
            return None

    @staticmethod
    def filter_and_sort_CAS(synonyms: list[str]) -> list[str]:
        CAS_rns = []

        if synonyms is not None:
            for synonym in synonyms:
                if synonym:
                    synonym = synonym.strip()
                    if MoleculeResolver.is_valid_CAS(synonym):
                        CAS_rns.append(synonym)

        if not CAS_rns:
            return CAS_rns

        counter = collections.Counter(CAS_rns)

        # if all CAS provided are the same, use the first one
        if len(counter) == 1:
            return [CAS_rns[0]]

        # if all CAS provided are as likely as the other use heuristic approach:
        # in my experience the CAS with the lowest first number is mostly the correct one.
        if all([v == 1 for v in counter.values()]):
            CAS_rns = sorted(CAS_rns, key=lambda CAS: int(CAS.split("-")[0]))

        return [counter.most_common()[0][0]]

    @staticmethod
    def check_charge(mol: Chem.rdchem.Mol, required_charge: Optional[int]) -> bool:
        if required_charge is None:
            return True

        MoleculeResolver._check_parameters(required_charges=required_charge)
        calculated_charge = Chem.rdmolops.GetFormalCharge(mol)

        charge_correct = False
        if isinstance(required_charge, str):
            if required_charge == "zero":
                charge_correct = calculated_charge == 0
            elif required_charge == "non_zero":
                charge_correct = calculated_charge != 0
            elif required_charge == "positive":
                charge_correct = calculated_charge > 0
            elif required_charge == "negative":
                charge_correct = calculated_charge < 0
        else:
            charge_correct = int(required_charge) == calculated_charge

        return charge_correct

    @staticmethod
    def formula_to_dictionary(
        formula: str, allow_duplicate_atoms: bool = True
    ) -> dict[str, Union[int, float]]:
        def merge_dictionaries(
            d1: dict[str, Union[int, float]],
            d2: dict[str, Union[int, float]],
            d2_multiplier: float = 1.0,
        ):
            d1 = d1.copy()
            for element, count in d2.items():
                if element not in d1:
                    d1[element] = 0
                d1[element] += count * d2_multiplier
            return d1

        def cast_multiplier(multiplier: str) -> float:
            if len(multiplier) > 0:
                return float(multiplier)
            else:
                return 1.0

        d = {}
        if "*" in formula:
            for formula_part in formula.split("*"):  # for hydrates and solvates
                temp_d = MoleculeResolver.formula_to_dictionary(formula_part)
                d = merge_dictionaries(d, temp_d)
                return d

        formula = formula.replace("[", "(").replace("{", "(")
        formula = formula.replace("]", ")").replace("}", ")")
        brackets_match = MoleculeResolver.formula_bracket_group_regex_compiled.search(
            formula
        )
        if brackets_match:
            rest_of_formula = (
                formula[0 : brackets_match.span()[0]]
                + formula[brackets_match.span()[1] :]
            )
            bracket_content = brackets_match.groups()[0][1:-1]
            bracket_multiplier = cast_multiplier(brackets_match.groups()[1])
            d = merge_dictionaries(
                MoleculeResolver.formula_to_dictionary(
                    rest_of_formula, allow_duplicate_atoms
                ),
                MoleculeResolver.formula_to_dictionary(
                    bracket_content, allow_duplicate_atoms
                ),
                bracket_multiplier,
            )
        else:
            matches = MoleculeResolver.empirical_formula_regex_compiled.findall(formula)
            for match in matches:
                element = match[0]
                multiplier = cast_multiplier(match[1])
                if allow_duplicate_atoms:
                    if element not in d:
                        d[element] = 0
                    d[element] += multiplier
                else:
                    if element in d:
                        raise ValueError("element twice in formula")
                    d[element] = multiplier

        only_integers_available = all([float(v).is_integer() for v in d.values()])
        conversion_function = int if only_integers_available else float
        for k, v in d.items():
            d[k] = conversion_function(v)

        return d

    @staticmethod
    def check_formula(
        mol_or_formula: Union[Chem.rdchem.Mol, str], required_formula: Optional[str]
    ) -> bool:
        if required_formula is None:
            return True

        if isinstance(mol_or_formula, str):
            mol_formula = mol_or_formula
        else:
            temp_mol = MoleculeResolver.standardize_molecule(mol_or_formula)
            mol_formula = rdMolDescriptors.CalcMolFormula(temp_mol)

        formula1 = MoleculeResolver.formula_to_dictionary(mol_formula)
        formula2 = MoleculeResolver.formula_to_dictionary(required_formula)

        return formula1 == formula2

    @staticmethod
    def check_molecular_mass(
        mol: Chem.rdchem.Mol,
        required_molecular_mass: Union[float, int],
        percentage_deviation_allowed: float = 0.001,
    ) -> bool:
        if not isinstance(required_molecular_mass, float) and not isinstance(
            required_molecular_mass, int
        ):
            raise TypeError("required_molecular_mass must be a float or an integer.")

        molecular_mass_from_mol = Descriptors.MolWt(mol)
        minimum_value = min(required_molecular_mass, molecular_mass_from_mol)
        absolute_relative_difference = abs(
            (required_molecular_mass - molecular_mass_from_mol) / minimum_value
        )
        return absolute_relative_difference < percentage_deviation_allowed

    @staticmethod
    def get_structure_type_from_SMILES(SMILES: str) -> str:
        SMILES_parts = SMILES.split(".")
        SMILES_parts_charges = []
        for SMILES_part in SMILES_parts:
            mol = MoleculeResolver.get_from_SMILES(SMILES_part)
            if not mol:
                return "None"
            SMILES_parts_charges.append(Chem.rdmolops.GetFormalCharge(mol))

        total_charge = sum(SMILES_parts_charges)

        if all([charge == 0 for charge in SMILES_parts_charges]):
            if len(SMILES_parts) > 1:
                structure_type = "mixture_neutrals"
            else:
                structure_type = "neutral"
        elif all([charge != 0 for charge in SMILES_parts_charges]):
            if len(SMILES_parts) > 1:
                if total_charge == 0:
                    structure_type = "salt"
                else:
                    structure_type = "mixture_ions"
            else:
                structure_type = "ion"
        else:
            if total_charge == 0:
                structure_type = "mixture_neutrals_salts"
            else:
                structure_type = "mixture_neutrals_ions"

        return structure_type

    @staticmethod
    def check_structure_type(SMILES: str, required_structure_type: str) -> bool:
        if required_structure_type is None:
            return True

        available_structure_types = [
            "mixture_neutrals",
            "mixture_ions",
            "neutral",
            "salt",
            "ion",
            "mixture_neutrals_salts",
            "mixture_neutrals_ions",
        ]
        if required_structure_type not in available_structure_types:
            raise ValueError(
                "required_structure_type must be one of 'mixture_neutrals','mixture_ions', 'neutral', 'salt', 'ion','mixture_neutrals_salts' or 'mixture_neutrals_ions'."
            )

        return (
            MoleculeResolver.get_structure_type_from_SMILES(SMILES)
            == required_structure_type
        )

    @staticmethod
    @cache
    def InChI_to_SMILES(inchi: str, standardize: bool = True) -> Optional[str]:
        mol = MoleculeResolver.get_from_InChI(inchi)

        if mol is None:
            return None

        return MoleculeResolver.standardize_SMILES(Chem.MolToSmiles(mol), standardize)

    @staticmethod
    @cache
    def SMILES_to_InChI(smiles: str, standardize: bool = True) -> Optional[str]:
        mol = MoleculeResolver.get_from_SMILES(smiles)

        if mol is None:
            return None

        return Chem.MolToInchi(MoleculeResolver.standardize_molecule(mol, standardize))

    @staticmethod
    @cache
    def get_from_InChI(inchi: str, addHs: bool = False) -> Optional[Chem.rdchem.Mol]:
        if inchi is None:
            return None
        if inchi == "":
            return None

        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            mol = Chem.MolFromInchi(inchi, sanitize=False)

        if mol is None:
            return None

        try:
            mol.UpdatePropertyCache()
        except Exception as error:
            mol.UpdatePropertyCache(strict=False)

        Chem.GetSymmSSSR(mol)

        if mol is None:
            return None

        if addHs:
            mol = Chem.AddHs(mol)

        return mol

    @staticmethod
    @cache
    def get_from_SMILES(SMILES: str, addHs: bool = False) -> Chem.rdchem.Mol:
        if SMILES is None:
            return None
        SMILES_parts = SMILES.split(".")

        if len(SMILES_parts) > 1:
            mol = Chem.Mol()
            for SMILES in SMILES_parts:
                mol2 = MoleculeResolver.get_from_SMILES(SMILES, addHs)
                if mol2 is None:
                    return None
                mol = Chem.CombineMols(mol, mol2)
            try:
                mol.UpdatePropertyCache()
            except Exception as error:
                mol.UpdatePropertyCache(strict=False)

            Chem.GetSymmSSSR(mol)
            return mol

        mol = Chem.MolFromSmiles(SMILES)
        if mol is None:
            # for some correct SMILES the sanitization process of rdkit does not work properly
            if (
                SMILES.count("F") > 0
                or SMILES.count("Cl") > 0
                or SMILES.count("Br") > 0
                or SMILES.count("B") > 0
            ):
                mol = Chem.MolFromSmiles(SMILES, sanitize=False)
                if mol is None:
                    return None
                mol.UpdatePropertyCache(strict=False)
        else:
            mol.UpdatePropertyCache()

        if mol is None:
            return None

        if addHs:
            mol = Chem.AddHs(mol)

        Chem.GetSymmSSSR(mol)
        return mol

    @staticmethod
    def show_and_pause(
        mol: Chem.rdchem.Mol,
        name: Optional[str] = None,
        size: tuple[int, int] = (1000, 1000),
    ):
        scaling_size = min(size)
        Draw.DrawingOptions.atomLabelFontSize = int(50 * scaling_size / 1000)
        Draw.DrawingOptions.dotsPerAngstrom = int(300 * scaling_size / 1000)
        Draw.DrawingOptions.bondLineWidth = max(float(4.0 * scaling_size / 1000), 1.0)
        img = Draw.MolToImage(mol, size=size)

        title = "charge: " + str(Chem.rdmolops.GetFormalCharge(mol))
        if name is not None:
            title = "name: " + name + "\n" + title

        draw = ImageDraw.Draw(img)
        s = int(size[1] / 30)
        fnt = ImageFont.truetype(font="arial.ttf", size=s)
        draw.multiline_text((10, 10), title, font=fnt, fill=(0, 0, 0, 255))

        img.show()

    @staticmethod
    def to_hill_formula(
        mol_or_dictionary: Union[Chem.rdchem.Mol, dict[str, Union[int, float]]]
    ) -> str:
        if isinstance(mol_or_dictionary, dict):
            temp = mol_or_dictionary
        else:
            mol = Chem.AddHs(mol_or_dictionary)
            temp = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol in temp:
                    temp[symbol] += 1
                else:
                    temp[symbol] = 1

        hill_formula = ""

        if "C" in temp:
            hill_formula = hill_formula + "C"

            if temp["C"] > 1:
                hill_formula = hill_formula + str(temp["C"])

            temp.pop("C")

            if "H" in temp:
                hill_formula = hill_formula + "H"

                if temp["H"] > 1:
                    hill_formula = hill_formula + str(temp["H"])

                temp.pop("H")

        for key in sorted(temp.keys()):
            hill_formula = hill_formula + key

            if temp[key] > 1:
                hill_formula = hill_formula + str(temp[key])

        return hill_formula

    @staticmethod
    def is_valid_SMILES_fast(SMILES: Optional[str]) -> bool:
        if not SMILES:
            return False

        if not isinstance(SMILES, str):
            return False

        if not regex.search(MoleculeResolver.non_generic_SMILES_regex_compiled, SMILES):
            return False

        return True

    @staticmethod
    @cache
    def is_valid_SMILES(SMILES: Optional[str]) -> bool:
        if not MoleculeResolver.is_valid_SMILES_fast(SMILES):
            return False

        return MoleculeResolver.get_from_SMILES(SMILES) != None

    @staticmethod
    @cache
    def is_valid_InChI(InChI: Optional[str]) -> bool:
        if not InChI:
            return False

        if not isinstance(InChI, str):
            return False

        if not regex.search(MoleculeResolver.InChICode_regex_compiled, InChI):
            return False

        return Chem.MolFromInchi(InChI) != None

    @staticmethod
    @cache
    def is_valid_CAS(cas: Union[str, bool, None]) -> bool:
        if not cas:
            return False

        if not isinstance(cas, str):
            return False

        # Takes into account the standard CAS formatting e.g. 7732-18-5
        cas_match = regex.search(MoleculeResolver.CAS_regex_with_groups, cas)
        if cas_match is None:
            return False
        cas_string = cas_match.group(1) + cas_match.group(2) + cas_match.group(3)

        increment = 0
        sum_cas = 0

        # Slices the reversed number string
        for number in reversed(cas_string):
            if increment == 0:
                validate = int(number)
            else:
                sum_cas = sum_cas + (int(number) * increment)

            increment = increment + 1

        # Does the math
        if validate == sum_cas % 10:
            return True
        else:
            return False

    @staticmethod
    @cache
    def are_InChIs_equal(
        InChI1: str,
        InChI2: str,
        standardize: bool = True,
        isomeric: bool = True,
        keep_fixed_Hs: bool = False,
    ) -> bool:
        if not InChI1 or not InChI2:
            return False

        if InChI1 == InChI2:
            return True

        mol1 = MoleculeResolver.get_from_InChI(InChI1)
        mol2 = MoleculeResolver.get_from_InChI(InChI2)
        if standardize:
            mol1 = MoleculeResolver.standardize_molecule(mol1)
            mol2 = MoleculeResolver.standardize_molecule(mol2)

        if not isomeric:
            mol1 = MoleculeResolver.get_from_SMILES(
                Chem.MolToSmiles(mol1, isomericSmiles=False)
            )
            mol2 = MoleculeResolver.get_from_SMILES(
                Chem.MolToSmiles(mol2, isomericSmiles=False)
            )

        options = "/FixedH" if keep_fixed_Hs else ""
        return Chem.MolToInchi(mol1, options=options) == Chem.MolToInchi(
            mol2, options=options
        )

    @staticmethod
    @cache
    def are_SMILES_equal(
        smiles1: str, smiles2: str, standardize: bool = True, isomeric: bool = True
    ) -> bool:
        if standardize:
            smiles1 = MoleculeResolver.standardize_SMILES(smiles1, True)
            smiles2 = MoleculeResolver.standardize_SMILES(smiles2, True)

        return MoleculeResolver.are_equal(
            MoleculeResolver.get_from_SMILES(smiles1),
            MoleculeResolver.get_from_SMILES(smiles2),
            False,
            isomeric,
        )

    @staticmethod
    @cache
    def get_resonance_SMILES(SMILES):
        mol = MoleculeResolver.get_from_SMILES(SMILES)
        rms = Chem.ResonanceMolSupplier(mol)
        rms.SetProgressCallback(EmptyResonanceMolSupplierCallback())
        resonance_mols = list(rms)

        # This workaround is needed in very few cornercases where the
        # resonance mol suppler returns molecules that are None
        # e.g. tetrabromocobaltate(II) [Co+2]([Br-])([Br-])([Br-])([Br-])
        if all([resonance_mol is None for resonance_mol in resonance_mols]):
            return [SMILES]

        return [Chem.MolToSmiles(mol_) for mol_ in resonance_mols if mol_]

    @staticmethod
    def are_equal(
        mol1: Chem.rdchem.Mol,
        mol2: Chem.rdchem.Mol,
        standardize: bool = True,
        isomeric: bool = True,
        check_for_resonance_structures: Optional[bool] = None,
        method: int = 1,
    ) -> bool:
        if mol1 is None or mol2 is None:
            return False

        if standardize:
            mol1 = MoleculeResolver.standardize_molecule(mol1)
            mol2 = MoleculeResolver.standardize_molecule(mol2)

        SMILES1 = Chem.MolToSmiles(mol1, isomericSmiles=isomeric)
        SMILES2 = Chem.MolToSmiles(mol2, isomericSmiles=isomeric)

        if check_for_resonance_structures is None:
            SMILES1_structure_type = MoleculeResolver.get_structure_type_from_SMILES(
                SMILES1
            )
            SMILES2_structure_type = MoleculeResolver.get_structure_type_from_SMILES(
                SMILES2
            )
            if SMILES1_structure_type != SMILES2_structure_type:
                return False
            check_for_resonance_structures = SMILES1_structure_type in [
                "ion",
                "salt",
            ] or SMILES2_structure_type in ["ion", "salt"]

        if method == 1:
            if SMILES1 == SMILES2:
                return True

            if check_for_resonance_structures:
                unique_partial_SMILES1 = set(SMILES1.split("."))
                unique_partial_SMILES2 = set(SMILES2.split("."))

                if len(unique_partial_SMILES1) != len(unique_partial_SMILES2):
                    return False

                matching_unique_partial_SMILES2_found = []
                for unique_partial_SMILES1_ in unique_partial_SMILES1:
                    resonance_SMILES1 = MoleculeResolver.get_resonance_SMILES(
                        unique_partial_SMILES1_
                    )
                    unique_partial_SMILES2 = unique_partial_SMILES2 - set(
                        matching_unique_partial_SMILES2_found
                    )
                    for unique_partial_SMILES2_ in unique_partial_SMILES2:
                        resonance_SMILES2 = MoleculeResolver.get_resonance_SMILES(
                            unique_partial_SMILES2_
                        )
                        if set(resonance_SMILES1) == set(resonance_SMILES2):
                            matching_unique_partial_SMILES2_found.append(
                                unique_partial_SMILES2_
                            )
                            break

                return len(unique_partial_SMILES1) == len(
                    matching_unique_partial_SMILES2_found
                )
                # https://github.com/rdkit/rdkit/discussions/4719
                # return rdMolHash.MolHash(standardized_molecule1, rdMolHash.HashFunction.HetAtomTautomer) == rdMolHash.MolHash(standardized_molecule1, rdMolHash.HashFunction.HetAtomTautomer)

            return False

        if method == 2:  # using the RDKit RegistrationHash

            def get_mol_hash(mol, hashscheme):
                return RegistrationHash.GetMolHash(
                    RegistrationHash.GetMolLayers(mol), hash_scheme=hashscheme
                )

            if isomeric:
                hashscheme = RegistrationHash.HashScheme.ALL_LAYERS
            else:
                hashscheme = RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS

            mol1_hash = get_mol_hash(mol1, hashscheme)
            mol2_hash = get_mol_hash(mol2, hashscheme)

            if mol1_hash != mol2_hash and check_for_resonance_structures:
                unique_partial_SMILES1 = set(SMILES1.split("."))
                unique_partial_SMILES2 = set(SMILES2.split("."))

                if len(unique_partial_SMILES1) != len(unique_partial_SMILES2):
                    return False

                matching_unique_partial_SMILES2_found = []
                for unique_partial_SMILES1_ in unique_partial_SMILES1:
                    unique_partial_SMILES2 = unique_partial_SMILES2 - set(
                        matching_unique_partial_SMILES2_found
                    )
                    for unique_partial_SMILES2_ in unique_partial_SMILES2:
                        resonance_SMILES1 = MoleculeResolver.get_resonance_SMILES(
                            unique_partial_SMILES1_
                        )
                        resonance_SMILES2 = MoleculeResolver.get_resonance_SMILES(
                            unique_partial_SMILES2_
                        )
                        if set(resonance_SMILES1) == set(resonance_SMILES2):
                            matching_unique_partial_SMILES2_found.append(
                                unique_partial_SMILES2_
                            )
                            break

                return len(unique_partial_SMILES1) == len(
                    matching_unique_partial_SMILES2_found
                )

            return mol1_hash == mol2_hash

    @staticmethod
    def find_duplicates_in_molecule_dictionary(
        molecules: dict[Any, Molecule], clustered_molecules=None
    ):
        if not clustered_molecules:
            clustered_molecules = MoleculeResolver.group_molecule_dictionary_by_formula(
                molecules
            )

        yield from MoleculeResolver.intersect_molecule_dictionaries(
            molecules,
            molecules,
            clustered_molecules,
            clustered_molecules,
            report_same_keys=False,
        )

    @staticmethod
    def group_molecule_dictionary_by_formula(
        molecules: dict[Any, Molecule]
    ) -> dict[str, dict[Any, Tuple[Chem.rdchem.Mol, Molecule]]]:
        temporary = {}
        for key, molecule in molecules.items():
            mol_molecule = MoleculeResolver.get_from_SMILES(molecule.SMILES)
            mol_formula = rdMolDescriptors.CalcMolFormula(mol_molecule)
            if mol_formula not in temporary:
                temporary[mol_formula] = {}

            temporary[mol_formula][key] = (mol_molecule, molecule)

        return temporary

    @staticmethod
    def intersect_molecule_dictionaries(
        molecules: dict[Any, Molecule],
        other_molecules: dict[Any, Molecule],
        mode: str = "SMILES",
        clustered_molecules: dict = None,
        clustered_other_molecules: dict = None,
        report_same_keys: bool = True,
    ) -> Tuple[Any, Molecule, Any, Molecule, bool]:
        if not clustered_molecules:
            clustered_molecules = MoleculeResolver.group_molecule_dictionary_by_formula(
                molecules
            )

        if not clustered_other_molecules:
            clustered_other_molecules = (
                MoleculeResolver.group_molecule_dictionary_by_formula(other_molecules)
            )

        for mol_formula, this_formula_molecules in clustered_molecules.items():
            if mol_formula in clustered_other_molecules:
                for key, (mol_molecule, molecule) in this_formula_molecules.items():
                    for other_key, (
                        other_mol_molecule,
                        other_molecule,
                    ) in clustered_other_molecules[mol_formula].items():
                        are_equal = False
                        same_isomer_information = True
                        if mode == "SMILES":
                            if MoleculeResolver.are_equal(
                                mol_molecule, other_mol_molecule
                            ):
                                are_equal = True
                                same_isomer_information = True
                            elif MoleculeResolver.are_equal(
                                mol_molecule, other_mol_molecule, isomeric=False
                            ):
                                are_equal = True
                                same_isomer_information = False
                        elif mode == "inchi":
                            if MoleculeResolver.are_InChIs_equal(
                                Chem.MolToInchi(mol_molecule),
                                Chem.MolToInchi(other_mol_molecule),
                            ):
                                are_equal = True
                                same_isomer_information = True
                            elif MoleculeResolver.are_InChIs_equal(
                                Chem.MolToInchi(mol_molecule),
                                Chem.MolToInchi(other_mol_molecule),
                                isomeric=False,
                            ):
                                are_equal = True
                                same_isomer_information = False
                        else:
                            raise ValueError('Supported modes are ["SMILES", "inchi"].')

                        if are_equal:
                            if report_same_keys or key != other_key:
                                yield (
                                    key,
                                    molecule,
                                    other_key,
                                    other_molecule,
                                    same_isomer_information,
                                )

    @staticmethod
    @cache
    def check_SMILES(
        SMILES: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
    ) -> bool:
        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
        )

        if not MoleculeResolver.is_valid_SMILES_fast(SMILES):
            return False

        mol = MoleculeResolver.get_from_SMILES(SMILES)

        if mol is None:
            return False

        if not MoleculeResolver.check_formula(mol, required_formula):
            return False

        if not MoleculeResolver.check_charge(mol, required_charge):
            return False

        if not MoleculeResolver.check_structure_type(SMILES, required_structure_type):
            return False

        return True

    def get_SMILES_from_Mol_format(
        self,
        *,
        molblock: Optional[str] = None,
        url: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[str]:
        if molblock is None and url is None:
            raise ValueError("molblock and url cannot both be None.")
        if molblock:
            if url is not None:
                raise
            if not isinstance(molblock, str):
                raise TypeError("molblock must be a string.")

        if url:
            if molblock is not None:
                raise
            if not isinstance(url, str):
                raise TypeError("url must be a string.")
            molblock = self._resilient_request(url)
            if molblock is None:
                return None

        # https://github.com/rdkit/rdkit/issues/1361
        def fix_MolBlock(mol_block):
            parts = mol_block.split("\n")
            count_line = parts[3].replace("v2000", "V2000").replace("v3000", "V3000")

            change_was_made = True

            ending = ""
            if count_line[-1] == "\r":
                ending = "\r"
                count_line = count_line.rstrip()

            while change_was_made:
                change_was_made = False
                three_digit_blocks = [
                    count_line[i : i + 3] for i in range(0, len(count_line), 3)
                ]

                for i_block, three_digit_block in enumerate(three_digit_blocks):
                    if three_digit_block != "   " and three_digit_block[-1] == " ":
                        count_line = (
                            count_line[: i_block * 3] + " " + count_line[i_block * 3 :]
                        )
                        change_was_made = True
                        break

            parts[3] = count_line + ending

            return "\n".join(parts)

        molblock = fix_MolBlock(molblock)
        mol = Chem.MolFromMolBlock(molblock)
        if not mol:
            mol = Chem.MolFromMolBlock(molblock, sanitize=False)
        if not mol:
            return None
        if standardize:
            mol = MoleculeResolver.standardize_molecule(mol)
        return Chem.MolToSmiles(mol)

    def get_SMILES_from_image_file(
        self,
        image_path: str,
        engines_order: list[str] = ["osra", "molvec", "imago"],
        mode: str = "single",
        standardize: bool = True,
    ) -> str:
        # usually OSRA works best, this is why it is left as default
        if mode not in ["single", "all"]:
            raise ValueError("The modes single and all are allowed.")

        SMILES = []

        with open(image_path, "rb") as handle:
            try:
                response_text = self._resilient_request(
                    "https://molvec.ncats.io/all",
                    {"data": handle.read()},
                    request_type="post",
                )

                data = json.loads(response_text)
                for engine in engines_order:
                    if data[engine]["status"].lower() in ["success", "cached"]:
                        temp_SMILES = self.get_SMILES_from_Mol_format(
                            molblock=data[engine]["molfile"]
                        )
                        if temp_SMILES is not None:
                            SMILES.append(temp_SMILES)
                            if mode == "single":
                                break
            except Exception:
                pass

        if mode == "single" and len(SMILES) > 0:
            SMILES = SMILES[0]

        SMILES = MoleculeResolver.standardize_SMILES(SMILES, standardize)

        return SMILES

    @staticmethod
    def save_to_PNG(
        mol: Chem.rdchem.Mol,
        filename: str,
        atom_infos=None,
        atom_infos_format_string: str = "%.3f",
        size: set[int, int] = (1000, 1000),
    ):
        new_mol = copy.deepcopy(mol)

        for i, atom in enumerate(new_mol.GetAtoms()):
            label = atom_infos_format_string % atom_infos[i]
            atom.SetProp("atomNote", label)

        Draw.MolToFile(new_mol, filename, size=size)

    @staticmethod
    def normalize_html(html_):
        html_ = regex.sub("\s+", " ", html_)
        html_ = html_.replace("&nbsp;", " ")
        html_ = regex.sub("\s+", " ", html_)
        return html_

    @staticmethod
    def parse_items_from_html(
        html_: str,
        split_tag: str,
        properties_regex: list[Tuple[str, str, list[int]]],
        property_indices_required: Optional[list[int]] = None,
    ):
        html_ = MoleculeResolver.normalize_html(html_)
        # HTMLparser was thought of here, but IMHO it is way to complicated
        # to match the correct opening and closing tags with the python std library
        # I did not want to depend on beautifulsoup at the beginning, but it is
        # a dependency that I am willing to add in the future.
        # The following code is a simple and quick solution that works for the
        # time being.
        if split_tag:
            html_parts = html_.split(split_tag)[1:]
        else:
            html_parts = [html_]

        all_items = []
        for html_part in html_parts:
            text_part = regex.sub(r"<.*?>", "", html_part)
            text_part = regex.sub(r"\s+", " ", text_part)
            item = []
            for source_type, property_regex, allowed_match_number in properties_regex:
                source = html_part if source_type == "html" else text_part
                property_matches = regex.findall(property_regex, source)
                if allowed_match_number:
                    if len(property_matches) not in allowed_match_number:
                        raise RuntimeError(
                            "The webpage either changed its format or the regex is picking up something unwanted."
                        )
                value = property_matches[0].strip() if property_matches else None
                item.append(value)

            add_item = False
            if property_indices_required:
                add_item = all(
                    [
                        item[property_index] is not None
                        for property_index in property_indices_required
                    ]
                )

            if add_item:
                all_items.append(item)

        return all_items

    @staticmethod
    def get_java_path():
        if platform.system().lower() == "windows":
            search_command = "where"
        elif platform.system().lower() in ["linux", "darwin"]:
            search_command = "which"
        else:
            raise NotImplementedError(
                f"For the following OS, getting the full path of the java executable needs to be programmed: {platform.system()}"
            )

        output = subprocess.run([search_command, "java"], capture_output=True)
        java_paths = output.stdout.decode("utf-8").strip().split(os.linesep)
        java_paths = [
            java_path for java_path in java_paths if os.path.exists(java_path)
        ]

        if output.returncode != 0 or len(java_paths) == 0:
            return None
        else:
            return java_paths[0]

    def get_molecule_from_OPSIN(
        self,
        name: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
        allow_warnings: bool = False,
    ) -> Optional[Molecule]:
        SMILES = None
        additional_information = ""
        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            services="opsin",
        )

        with self.query_molecule_cache("opsin", "name", name) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                response_text = self._resilient_request(
                    f"https://opsin.ch.cam.ac.uk/opsin/{urllib.parse.quote(name)}"
                )

                if response_text is not None:
                    temp = json.loads(response_text)
                    if temp["status"] == "SUCCESS" or (
                        allow_warnings and temp["status"] == "WARNING"
                    ):
                        if MoleculeResolver.check_SMILES(
                            temp["smiles"],
                            required_formula,
                            required_charge,
                            required_structure_type,
                        ):
                            SMILES = temp["smiles"]
                            SMILES = MoleculeResolver.standardize_SMILES(
                                SMILES, standardize
                            )
                        else:
                            SMILES_from_InChI = MoleculeResolver.InChI_to_SMILES(
                                temp["inchi"], standardize
                            )
                            if MoleculeResolver.check_SMILES(
                                SMILES_from_InChI,
                                required_formula,
                                required_charge,
                                required_structure_type,
                            ):
                                SMILES = MoleculeResolver.standardize_SMILES(
                                    SMILES_from_InChI, standardize
                                )

                        if "warnings" in temp:
                            additional_information = "WARNINGS: " + ", ".join(
                                temp["warnings"]
                            )

                        if SMILES:
                            molecules.append(
                                Molecule(
                                    SMILES,
                                    service="opsin",
                                    mode="name",
                                    additional_information=additional_information,
                                )
                            )

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def get_molecule_from_OPSIN_batchmode(
        self, names: list[str], standardize: bool = True
    ) -> list[Optional[Molecule]]:
        MoleculeResolver._check_parameters(
            identifiers=names,
            modes="name",
            services="opsin",
            context="get_molecules_batch",
        )

        if not self._java_path:
            raise FileNotFoundError(
                "The java installation could not be found. Either it is not installed or its location has not been added to the path environment variable."
            )

        def download_OPSIN_offline_and_convert_names(names: list[str], tempfolder: str):
            if not os.path.exists(os.path.join(tempfolder, "opsin.jar")):
                response_text = self._resilient_request(
                    "https://api.github.com/repos/dan2097/opsin/releases/latest"
                )
                temp = json.loads(response_text)
                for asset in temp["assets"]:
                    if asset["name"].count("opsin-cli"):
                        download_url = asset["browser_download_url"]
                        with closing(urllib.request.urlopen(download_url)) as r:
                            with open(os.path.join(tempfolder, "opsin.jar"), "wb") as f:
                                for chunk in r:
                                    f.write(chunk)

            unique_id = str(uuid.uuid4())  # needed for multiple runs in parallel
            input_file = os.path.join(tempfolder, f"input_{unique_id}.txt")
            with open(input_file, "w", encoding="utf8") as f:
                f.write("\n".join(names))

            # in the future adding the parameter -s would allow getting more structures if the stereo information opsin cannot interpret can be ignored
            _ = subprocess.run(
                [
                    self._java_path,
                    "-jar",
                    "opsin.jar",
                    "-osmi",
                    f"input_{unique_id}.txt",
                    f"output_{unique_id}.txt",
                ],
                capture_output=True,
                cwd=tempfolder,
            )

            with open(os.path.join(tempfolder, f"output_{unique_id}.txt"), "r") as f:
                output = f.read()

            output_values = output.split("\n")
            # I don't know in what situations, but in some OPSIN adds an empty line at the end of the output file.
            if len(output_values) == len(names) + 1 and output_values[-1] == "":
                output_values = output_values[:-1]

            SMILES = [value.strip() for value in output_values]
            if len(SMILES) != len(names):
                raise RuntimeError(
                    "There was a problem parsing the OPSIN offline output file."
                )

            return SMILES

        with self.query_molecule_cache_batchmode("opsin", "name", names) as (
            identifiers_to_search,
            indices_of_identifiers_to_search,
            results,
        ):
            if len(identifiers_to_search) == 0:
                return results

            if self._OPSIN_tempfolder is None or not os.path.exists(
                self._OPSIN_tempfolder.name
            ):
                with tempfile.TemporaryDirectory() as tempfolder:
                    SMILES = download_OPSIN_offline_and_convert_names(
                        identifiers_to_search, tempfolder
                    )
            else:
                SMILES = download_OPSIN_offline_and_convert_names(
                    identifiers_to_search, self._OPSIN_tempfolder.name
                )

            for molecule_index, name, smi in zip(
                indices_of_identifiers_to_search,
                identifiers_to_search,
                SMILES,
            ):
                smi = MoleculeResolver.standardize_SMILES(smi, standardize)
                if smi:
                    results[molecule_index] = [
                        Molecule(
                            SMILES=smi,
                            synonyms=[name],
                            mode="name",
                            service="opsin",
                            identifier=name,
                        )
                    ]

            return results

    @cache
    def get_molecule_for_ion_from_partial_pubchem_search(
        self,
        identifier: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        standardize: bool = True,
    ) -> Optional[list]:
        required_structure_type = "ion"
        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
        )

        search_response_text = self._resilient_request(
            f'https://www.ncbi.nlm.nih.gov/pccompound/?term={urllib.parse.quote(identifier, safe="")}'
        )
        found_SMILES = []
        if search_response_text is not None:
            items = MoleculeResolver.parse_items_from_html(
                search_response_text,
                '<div class="rprt">',
                [("text", "CID: (\d+)", [1])],
                [0],
            )
            for (cid,) in items:
                molecule = self.get_molecule_from_pubchem(cid, "cid")
                is_salt_or_mixture = molecule.SMILES.count(".") > 0

                if is_salt_or_mixture:
                    SMILES_list = molecule.SMILES.split(".")
                    for temptative_SMILES in SMILES_list:
                        temp_mol = MoleculeResolver.get_from_SMILES(temptative_SMILES)
                        if MoleculeResolver.check_formula(
                            temp_mol, required_formula
                        ) and MoleculeResolver.check_charge(temp_mol, required_charge):
                            SMILES = Chem.MolToSmiles(temp_mol)
                            SMILES = MoleculeResolver.standardize_SMILES(
                                SMILES, standardize
                            )
                            found_SMILES.append(SMILES)

        if len(found_SMILES) == 0:
            return None

        # sort to most likely candidate
        counts = collections.Counter(found_SMILES)

        return [
            (Molecule(smi, service="pubchem"), y) for smi, y in counts.most_common()
        ]

    def get_molecule_from_ChEBI(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        if required_formula is None:
            if mode == "formula":
                required_formula = identifier

        MoleculeResolver._check_parameters(
            services="chebi",
            modes=mode,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            required_formulas=required_formula,
        )

        with self.query_molecule_cache("chebi", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                CHEBI_URL = "https://www.ebi.ac.uk/webservices/chebi/2.0/test/"
                mode_mapping = {
                    "name": "ALL+NAMES",
                    "cas": "REGISTRY+NUMBERS",
                    "formula": "FORMULA",
                    "smiles": "SMILES",
                    "inchi": "INCHI%2FINCHI+KEY",
                    "inchikey": "INCHI%2FINCHI+KEY",
                }
                maximumResults = 5 if mode == "formula" else 1
                search_response_text = self._resilient_request(
                    f'{CHEBI_URL}getLiteEntity?search={urllib.parse.quote(identifier, safe="")}&searchCategory={mode_mapping[mode]}&maximumResults={maximumResults}&starsCategory=ALL'
                )

                SMILES = None
                synonyms = []
                CAS = []
                ChEBI_id = None

                if search_response_text is not None:
                    root = xmltodict.parse(search_response_text)
                    if "getLiteEntityResponse" not in root["S:Envelope"]["S:Body"]:
                        return None

                    temp = root["S:Envelope"]["S:Body"]["getLiteEntityResponse"][
                        "return"
                    ]
                    if not temp:
                        return None

                    temp_list = temp["ListElement"]
                    if not isinstance(temp["ListElement"], list):
                        temp_list = [temp["ListElement"]]

                    for ChEBI_id in [item["chebiId"] for item in temp_list]:
                        molecule_response_text = self._resilient_request(
                            f"{CHEBI_URL}getCompleteEntity?chebiId={ChEBI_id}"
                        )

                        if molecule_response_text is not None:
                            root = xmltodict.parse(molecule_response_text)
                            molecule = root["S:Envelope"]["S:Body"][
                                "getCompleteEntityResponse"
                            ]["return"]

                            if "smiles" not in molecule:
                                continue
                            SMILES = molecule["smiles"]

                            if "IupacNames" in molecule:
                                if isinstance(molecule["IupacNames"], list):
                                    synonyms.extend(
                                        [
                                            synonym["data"]
                                            for synonym in molecule["IupacNames"]
                                        ]
                                    )
                                elif isinstance(molecule["IupacNames"], dict):
                                    synonyms.append(molecule["IupacNames"]["data"])
                            if "Synonyms" in molecule:
                                temp = molecule["Synonyms"]
                                if not isinstance(molecule["Synonyms"], list):
                                    temp = [temp]
                                synonyms.extend([synonym["data"] for synonym in temp])

                            synonyms = MoleculeResolver.filter_and_sort_synonyms(
                                synonyms
                            )

                            if "RegistryNumbers" in molecule:
                                temp = molecule["RegistryNumbers"]
                                if not isinstance(molecule["RegistryNumbers"], list):
                                    temp = [molecule["RegistryNumbers"]]
                                CAS = MoleculeResolver.filter_and_sort_CAS(
                                    [synonym["data"] for synonym in temp]
                                )

                            molecule = Molecule(
                                SMILES,
                                synonyms,
                                CAS,
                                ChEBI_id.split(":")[1],
                                mode,
                                service="chebi",
                            )
                            molecules.append(molecule)

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def get_molecules_using_batchmode_from(
        self,
        identifiers: str,
        modes: str,
        service: str,
        batch_size: int = 1000,
        progressbar: bool = False,
        ignore_exceptions: bool = False,
    ) -> tuple[dict[str, list[Optional[Molecule]]], list[str]]:
        if service not in self._available_services_with_batch_capabilities:
            raise ValueError(
                f"The service {service} does not have batch capabilities or is not supported."
            )

        all_supported_modes = []
        new_modes = []
        for _modes in modes:
            _modes = [mode.lower() for mode in _modes]
            new_modes.append(_modes)
            [
                all_supported_modes.append(mode)
                for mode in _modes
                if mode not in all_supported_modes
                and mode in MoleculeResolver._supported_modes_by_services[service]
            ]

        modes = new_modes

        MoleculeResolver._check_parameters(
            identifiers=identifiers, modes=modes, context="batch"
        )

        identifier_sets = {}
        identifier_indices_sets = {}
        for mode in all_supported_modes:
            identifier_sets[mode] = []
            identifier_indices_sets[mode] = []

        for i_molecule, (i_molecule_identifiers, i_molecule_modes) in enumerate(
            zip(identifiers, modes)
        ):
            if len(i_molecule_identifiers) != len(i_molecule_modes):
                raise ValueError(
                    "The list of identifiers and modes of a molecule must be of the same length."
                )
            (
                flattened_i_molecule_identifiers,
                flattened_i_molecule_modes,
                _,
                _,
                _,
            ) = MoleculeResolver._check_and_flatten_identifiers_and_modes(
                i_molecule_identifiers, i_molecule_modes
            )
            for identifier, mode in zip(
                flattened_i_molecule_identifiers, flattened_i_molecule_modes
            ):
                if mode in all_supported_modes:
                    identifier_sets[mode].append(identifier)
                    identifier_indices_sets[mode].append(i_molecule)

        results = {}
        for mode in all_supported_modes:
            this_mode_results = {}

            _identifiers = identifier_sets[mode]
            _identifier_indices = identifier_indices_sets[mode]

            chunks = list(
                MoleculeResolver.chunker(
                    list(zip(_identifiers, _identifier_indices)), batch_size
                )
            )

            if chunks:
                for chunk in tqdm(chunks, disable=not progressbar):
                    identifier_indices_chunk = [val[1] for val in chunk]
                    identifier_chunk = [val[0] for val in chunk]

                    if service == "pubchem":
                        results_for_this_chunk = (
                            self.get_molecules_from_pubchem_batchmode(
                                identifier_chunk, mode, True
                            )
                        )
                    elif service == "srs":
                        results_for_this_chunk = self.get_molecules_from_SRS_batchmode(
                            identifier_chunk, mode, True
                        )
                    elif service == "comptox":
                        results_for_this_chunk = (
                            self.get_molecules_from_CompTox_batchmode(
                                identifier_chunk, mode, True
                            )
                        )
                    elif service == "opsin":
                        results_for_this_chunk = self.get_molecule_from_OPSIN_batchmode(
                            identifier_chunk, True
                        )
                    else:
                        raise NotImplementedError(
                            "The batch functionality for this service is either not implemented or the service does not support batch searches."
                        )

                    for index, result in zip(
                        identifier_indices_chunk, results_for_this_chunk
                    ):
                        if index not in this_mode_results:
                            this_mode_results[index] = None

                        if result and this_mode_results[index] is None:
                            this_mode_results[index] = result

                    if (
                        len(
                            set(identifier_indices_chunk)
                            - set(this_mode_results.keys())
                        )
                        != 0
                    ):
                        raise RuntimeError(
                            "The number of returned values should be equal to the unique values asked for."
                        )

            results[mode] = []
            for index in range(len(identifiers)):
                if index in _identifier_indices:
                    if index in this_mode_results:
                        results[mode].append(this_mode_results[index])
                    else:
                        results[mode].append(None)

        return results, all_supported_modes

    def get_molecule_from_CompTox(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        # new API, at some point we need to change
        # https://api-ccte.epa.gov/docs/chemical.html#/

        MoleculeResolver._check_parameters(
            modes=mode,
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            services="comptox",
        )

        with self.query_molecule_cache("comptox", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                COMPTOX_URL = "https://comptox.epa.gov/dashboard-api/"

                response_text = self._resilient_request(
                    f'{COMPTOX_URL}ccdapp1/search/chemical/equal/{urllib.parse.quote(identifier, safe="")}',
                    rejected_status_codes=[400, 404],
                )

                if response_text is not None:
                    original_response = json.loads(response_text)
                    # sort and filter best results
                    original_response = sorted(
                        original_response, key=lambda x: x["rank"]
                    )
                    original_response = list(
                        filter(
                            lambda x: x["rank"] == original_response[0]["rank"],
                            original_response,
                        )
                    )

                    def process_dtxsid(temp_dtxsid):
                        substance_response_text = self._resilient_request(
                            f"{COMPTOX_URL}ccdapp2/chemical-detail/search/by-dsstoxsid/?id={urllib.parse.quote(temp_dtxsid)}"
                        )

                        if substance_response_text is not None:
                            temp_substance = json.loads(substance_response_text)

                            temp_SMILES = temp_substance["smiles"]
                            temp_SMILES = MoleculeResolver.standardize_SMILES(
                                temp_SMILES, standardize
                            )
                            if not temp_SMILES:
                                return
                            temp_synonyms = []
                            temp_CAS = []

                            if "casrn" in temp_substance:
                                temp_CAS.append(temp_substance["casrn"])

                            if "preferredName" in temp_substance:
                                if temp_substance["preferredName"] is not None:
                                    temp_synonyms.append(
                                        temp_substance["preferredName"]
                                    )
                            if "acdIupacName" in temp_substance:
                                if temp_substance["acdIupacName"] is not None:
                                    temp_synonyms.append(temp_substance["acdIupacName"])
                            if "wikipediaName" in temp_substance:
                                if temp_substance["wikipediaName"] is not None:
                                    temp_synonyms.append(
                                        temp_substance["wikipediaName"]
                                    )

                            if mode == "name":
                                temp_synonyms.append(identifier)

                            QC_LEVEL_str = ""
                            if "qcLevel" in temp_substance:
                                QC_LEVEL_str = (
                                    f'|QC_LEVEL:{float(temp_substance["qcLevel"])}'
                                )
                            temp_synonyms = MoleculeResolver.filter_and_sort_synonyms(
                                temp_synonyms
                            )
                            molecules.append(
                                Molecule(
                                    temp_SMILES,
                                    temp_synonyms,
                                    temp_CAS,
                                    temp_dtxsid + QC_LEVEL_str,
                                    mode,
                                    service="comptox",
                                )
                            )

                    if isinstance(original_response, list):
                        for i in range(len(original_response)):
                            temp = original_response[i]
                            temp_dtxsid = temp["dtxsid"]
                            process_dtxsid(temp_dtxsid)

                    else:
                        temp_dtxsid = original_response["dtxsid"]
                        process_dtxsid(temp_dtxsid)

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def get_molecule_from_CTS(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        MoleculeResolver._check_parameters(
            identifiers=identifier,
            modes=mode,
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            services="cts",
        )

        with self.query_molecule_cache("cts", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                if "CTS_is_down" in self._message_slugs_shown:
                    return None

                cts_modes = {
                    "name": "Chemical Name",
                    "cas": "CAS",
                    "inchi": "smiles",
                    "smiles": "smiles",
                }

                SMILES = None
                CAS = []
                synonyms = []

                mode_used = mode
                if mode == "smiles":
                    SMILES = identifier

                elif mode == "inchi":
                    SMILES = MoleculeResolver.InChI_to_SMILES(identifier, standardize)
                    identifier = SMILES
                    mode_used = "smiles calculated from inchi"

                elif mode == "cas":
                    CAS.append(identifier)

                # unfortunatedly sometimes the server returns a status code 500 for valid names
                try:
                    CTS_URL = "https://cts.fiehnlab.ucdavis.edu/rest/convert/"
                    response_text = self._resilient_request(
                        f'{CTS_URL}{urllib.parse.quote(cts_modes[mode])}/Chemical%20Name/{urllib.parse.quote(identifier, safe="")}',
                        kwargs={"timeout": 5},
                        rejected_status_codes=[404, 500],
                        max_retries = 2
                    )
                except requests.exceptions.ConnectionError:
                    # I don't know why, but somtimes CTS is offline. This would make the module much slower as
                    # it tries to connect multiple times anyway. Instead we give a warning and skip CTS.
                    if "CTS_is_down" not in self._message_slugs_shown:
                        self._message_slugs_shown.append("CTS_is_down")
                        warnings.warn(
                            "CTS seems to be down, to continue working this instance of MoleculeResolver will skip CTS."
                        )

                if "CTS_is_down" not in self._message_slugs_shown:
                    # the search for synonyms returns a lot of unusable data
                    # only use the three most sensible ones
                    if response_text is not None:
                        temp = json.loads(response_text)[0]
                        synonyms.extend(
                            self.filter_and_sort_synonyms(temp["results"], 3)
                        )

                    if mode != "cas":
                        response_text = self._resilient_request(
                            f'{CTS_URL}{urllib.parse.quote(cts_modes[mode])}/CAS/{urllib.parse.quote(identifier, safe="")}',
                            kwargs={"timeout": 10},
                            rejected_status_codes=[404, 500],
                        )
                        if response_text is not None:
                            CAS_rns = json.loads(response_text)[0]["results"]
                            CAS = MoleculeResolver.filter_and_sort_CAS(CAS_rns)

                    if mode not in ["inchi", "smiles"]:
                        response_text = self._resilient_request(
                            f'{CTS_URL}{urllib.parse.quote(cts_modes[mode])}/InChI%20Code/{urllib.parse.quote(identifier, safe="")}',
                            kwargs={"timeout": 10},
                            rejected_status_codes=[404, 500],
                        )
                        if response_text is not None:
                            temp = json.loads(response_text)[0]
                            found_InChIs = temp["results"]
                            accepted_SMILES = []
                            for InChI in found_InChIs:
                                this_SMILES = MoleculeResolver.InChI_to_SMILES(
                                    InChI, standardize
                                )
                                if not this_SMILES:
                                    continue
                                if (
                                    not required_structure_type
                                    or "mixture" not in required_structure_type
                                ):
                                    if "." in this_SMILES:
                                        continue
                                # try filtering radicals as CTS has them in some cases
                                try:
                                    if (
                                        Descriptors.NumRadicalElectrons(
                                            self.get_from_SMILES(this_SMILES)
                                        )
                                        == 0
                                    ):
                                        accepted_SMILES.append(this_SMILES)
                                except Exception:
                                    continue

                            if len(accepted_SMILES) == 1:
                                if mode == "name":
                                    synonyms.insert(0, identifier)
                                    synonyms = self.filter_and_sort_synonyms(synonyms)
                                SMILES = accepted_SMILES[0]
                    if SMILES is not None:
                        molecules.append(
                            Molecule(SMILES, synonyms, CAS, None, mode_used, "cts")
                        )

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    @staticmethod
    def get_CompTox_request_unique_id() -> str:
        def base36encode(
            number, alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ".lower()
        ):
            """Converts an integer to a base36 string."""
            if not isinstance(number, (int)):
                raise TypeError("number must be an integer")

            base36 = ""
            sign = ""

            if number < 0:
                sign = "-"
                number = -number

            if 0 <= number < len(alphabet):
                return sign + alphabet[number]

            while number != 0:
                number, i = divmod(number, len(alphabet))
                base36 = alphabet[i] + base36

            return sign + base36

        return base36encode(int(time.time() * 1000))

    def get_molecules_from_CompTox_batchmode(
        self, identifiers: list[str], mode: str, standardize: bool = True
    ) -> list[Optional[Molecule]]:
        MoleculeResolver._check_parameters(
            identifiers=identifiers,
            modes=mode,
            services="comptox",
            context="get_molecules_batch",
        )
        if not all([isinstance(identifier, str) for identifier in identifiers]):
            raise TypeError("All identifiers must be strings.")

        with self.query_molecule_cache_batchmode(
            "comptox", mode, identifiers, save_not_found=False
        ) as (identifiers_to_search, indices_of_identifiers_to_search, results):
            if len(identifiers_to_search) == 0:
                return results

            COMPTOX_URL = "https://comptox.epa.gov/dashboard-api/batchsearch/export/"
            comptox_modes = {
                "name": "chemical_name",
                "cas": "CASRN",
                "inchikey": "INCHIKEY",
            }

            postdata = {
                "downloadItems": [
                    "SYNONYM_IDENTIFIER",
                    "QC_LEVEL",
                    "CASRN",
                    "IUPAC_NAME",
                    "SMILES",
                ],
                "downloadType": "EXCEL",
                "identifierTypes": [comptox_modes[mode]],
                "inputType": "IDENTIFIER",
                "massError": 0,
                "searchItems": "\n".join(identifiers_to_search)
                #'qc_level', 'expocast', 'data_sources', 'toxvaldata', 'assays', 'number_of_pubmed_articles', 'number_of_pubchem_data_sources', 'number_of_cpdat_sources', 'in_iris_list', 'in_pprtv_list', 'in_wikipedia_list', 'qsar_ready_smiles', 'ms_ready_smiles', 'synonym_identifier'
            }

            def poll_request(job_id):
                download_url = None
                n_try = 0
                while True:
                    poll_response_text = self._resilient_request(
                        f"{COMPTOX_URL}/status/{job_id}/?{self.get_CompTox_request_unique_id()}"
                    )

                    if poll_response_text.strip().lower() == "true":
                        download_url = f"{COMPTOX_URL}content/{job_id}"
                        break

                    if n_try > 15:
                        break

                    n_try += 1

                    time.sleep(2)

                return download_url

            post_request_text = self._resilient_request(
                f"{COMPTOX_URL}/?{MoleculeResolver.get_CompTox_request_unique_id()}",
                json=postdata,
                request_type="post",
                accepted_status_codes=[202],
            )  # allow_redirects=True)

            download_url = poll_request(post_request_text.strip())

            # QC_Level : https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0247-6/tables/1
            # 1 Expert curated: highest confidence in accuracy and consistency of unique chemical identifiers
            # 2 Expert curated: unique chemical identifiers confirmed using multiple public sources
            # 3 Programmatically curated from high quality EPA source(s) and unique chemical identifiers have no conflicts in ChemIDPlus and PubChem
            # 4 Programmatically curated from ChemIDPlus. Unique chemical identifiers have no conflicts in PubChem
            # 5 Programmatically curated from ACToR or PubChem. Unique chemical identifiers have low confidence and have a single public source

            # in some rare occasions, this times out and returns None
            if download_url:
                request_response = self._resilient_request(
                    download_url, return_response=True
                )
                with (
                    tempfile.TemporaryDirectory() as temp_dir,
                    warnings.catch_warnings(),
                ):
                    warnings.simplefilter("ignore")

                    temp_path = os.path.join(
                        temp_dir, "moleculeresolver_comptox_postdata.xlsx"
                    )
                    with open(temp_path, "wb") as f:
                        f.write(request_response.content)

                    wb = openpyxl.load_workbook(temp_path)

                    synonym_sheet = wb.worksheets[2]
                    synonyms_by_identifier_lower = {}
                    for i_row, row in enumerate(synonym_sheet.iter_rows()):
                        if i_row == 0:
                            continue

                        row_values = [cell.value for cell in row]
                        synonyms = []
                        identifier = row_values[0].strip()
                        if row_values[1]:
                            synonyms = [v.strip() for v in row_values[1].split("|")]

                        synonyms_by_identifier_lower[identifier.lower()] = synonyms

                    temp_results_by_identifier = {}

                    def parse_CompTox_result(row):
                        identifier = str(row[0]).strip()
                        found_by = row[1].lower().strip()

                        SMILES = None
                        CAS = []
                        synonyms = []

                        if found_by.count("found 0 results") == 0:
                            SMILES = row[7].strip()
                            if not SMILES:
                                return
                            dtxsid = row[2].strip()
                            preferred_name = row[3].strip()
                            qc_level = row[4]
                            CAS = str(row[5]).strip()
                            iupac_name = row[6].strip()

                            if iupac_name and iupac_name != "N/A":
                                if iupac_name not in synonyms:
                                    synonyms.insert(0, iupac_name)

                            if preferred_name and preferred_name != "N/A":
                                if preferred_name not in synonyms:
                                    synonyms.insert(0, preferred_name)

                            if not MoleculeResolver.is_valid_CAS(CAS):
                                CAS = []
                            else:
                                CAS = [CAS]

                            if "N/A" in SMILES or not MoleculeResolver.is_valid_SMILES(
                                SMILES
                            ):
                                return

                            SMILES = MoleculeResolver.standardize_SMILES(
                                SMILES, standardize
                            )

                            # the rating is needed because for some substances it
                            # returns more than one row. It is used to get the best result later.
                            if "expert" in found_by:
                                own_rating = 10
                            elif "approved" in found_by:
                                own_rating = 8
                            elif "valid source" in found_by:
                                own_rating = 6
                            elif "systematic name" in found_by:
                                own_rating = 4
                            else:
                                own_rating = 2

                            if CAS:
                                own_rating += 1

                            if identifier not in temp_results_by_identifier:
                                temp_results_by_identifier[identifier] = []

                            if SMILES:
                                temp_results_by_identifier[identifier].append(
                                    (
                                        qc_level,
                                        own_rating,
                                        SMILES,
                                        synonyms,
                                        CAS,
                                        dtxsid + "|QC_LEVEL:" + str(qc_level),
                                    )
                                )

                        else:
                            temp_results_by_identifier[identifier] = None

                    results_sheet = wb.worksheets[1]
                    for i_row, row in enumerate(results_sheet.iter_rows()):
                        if i_row == 0:
                            continue

                        row_values = [cell.value for cell in row]
                        parse_CompTox_result(row_values)

                    for molecule_index, identifier in zip(
                        indices_of_identifiers_to_search, identifiers_to_search
                    ):
                        if identifier in temp_results_by_identifier:
                            temp_result = temp_results_by_identifier[identifier]
                            if temp_result:
                                if len(temp_result) > 1:
                                    best_results = temp_result
                                    try:
                                        best_qc_level = min(
                                            temp_result, key=lambda x: x[0]
                                        )[0]
                                        best_results = list(
                                            filter(
                                                lambda x: x[0] == best_qc_level,
                                                temp_result,
                                            )
                                        )
                                    except Exception:
                                        pass

                                    if len(best_results) > 1:
                                        best_own_rating = max(
                                            best_results, key=lambda x: x[1]
                                        )[1]
                                        best_results = list(
                                            filter(
                                                lambda x: x[1] == best_own_rating,
                                                best_results,
                                            )
                                        )

                                    temp_results_by_identifier[
                                        identifier
                                    ] = best_results
                                results[molecule_index] = []

                                for (
                                    _,
                                    _,
                                    SMILES,
                                    synonyms,
                                    CAS,
                                    additional_information,
                                ) in temp_results_by_identifier[identifier]:
                                    this_synonyms = synonyms.copy()
                                    for synonym in synonyms:
                                        if (
                                            synonym.lower()
                                            in synonyms_by_identifier_lower
                                        ):
                                            for (
                                                new_synonym
                                            ) in synonyms_by_identifier_lower[
                                                synonym.lower()
                                            ]:
                                                if new_synonym not in this_synonyms:
                                                    this_synonyms.append(new_synonym)
                                    synonyms = (
                                        MoleculeResolver.filter_and_sort_synonyms(
                                            this_synonyms
                                        )
                                    )
                                    results[molecule_index].append(
                                        Molecule(
                                            SMILES,
                                            synonyms,
                                            CAS,
                                            additional_information,
                                            mode,
                                            "comptox",
                                            identifier=identifier,
                                        )
                                    )

            return results

    def get_molecule_from_Chemeo(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        MoleculeResolver._check_parameters(
            identifiers=identifier,
            modes=mode,
            required_charges=required_charge,
            required_formulas=required_formula,
            required_structure_types=required_structure_type,
            services="chemeo",
        )

        with self.query_molecule_cache("chemeo", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                valid_token_found = (
                    self.available_service_API_keys["chemeo"] is not None
                )
                if valid_token_found:
                    valid_token_found = isinstance(
                        self.available_service_API_keys["chemeo"], str
                    )
                    valid_token_found = (
                        MoleculeResolver.chemeo_API_token_regex_compiled.match(
                            self.available_service_API_keys["chemeo"]
                        )
                    )

                if not valid_token_found:
                    if "chemeo_API_token_missing" not in self._message_slugs_shown:
                        self._message_slugs_shown.append("chemeo_API_token_missing")
                        warnings.warn(
                            "chemeo requires a valid API token to allow search. After registering please insert the token in the corresponding variable when initializing the class instance."
                        )

                    return None

                CHEMEO_URL = "https://www.chemeo.com/api/v1/"
                API_bearer_headers = {}
                API_bearer_headers["Accept"] = "application/json"
                API_bearer_headers[
                    "Authorization"
                ] = f'Bearer {self.available_service_API_keys["chemeo"]}'

                request_text = self._resilient_request(
                    f'{CHEMEO_URL}convert/{mode}/{urllib.parse.quote(identifier, safe="")}',
                    kwargs={"headers": API_bearer_headers},
                    rejected_status_codes=[403, 404],
                )

                if (
                    request_text is not None
                    and request_text != '{"message":"Compound not found."}'
                ):
                    temp = json.loads(request_text)
                    request_text = self._resilient_request(
                        f'{CHEMEO_URL}search?q={temp["cid"]}',
                        kwargs={"headers": API_bearer_headers},
                        rejected_status_codes=[403, 404],
                    )

                    temp = json.loads(request_text)

                    for i in range(temp["total"]):
                        result = temp["comps"][i]

                        temp_synonyms = []
                        temp_SMILES = None

                        if "other_names" in result:
                            temp_synonyms.extend(result["other_names"])
                        temp_synonyms = MoleculeResolver.filter_and_sort_synonyms(
                            temp_synonyms
                        )

                        temp_CAS = []
                        if "cas" in result:
                            temp_CAS = result["cas"]
                            if not MoleculeResolver.is_valid_CAS(temp_CAS):
                                temp_CAS = []
                            else:
                                temp_CAS = [temp_CAS]

                        # using the mol block first instead of the readily available SMILES,
                        # because the later one is not an isomeric SMILES
                        if "mol3d" in result:
                            temp_SMILES = self.get_SMILES_from_Mol_format(
                                molblock=result["mol3d"], standardize=standardize
                            )

                        if temp_SMILES is None and "mol2d" in result:
                            temp_SMILES = self.get_SMILES_from_Mol_format(
                                molblock=result["mol2d"], standardize=standardize
                            )

                        if temp_SMILES is None and "inchi" in result:
                            temp_SMILES = MoleculeResolver.InChI_to_SMILES(
                                result["inchi"], standardize
                            )

                        if temp_SMILES is None and "smiles" in result:
                            temp_SMILES = result["smiles"]

                        if temp_SMILES is None:
                            raise ValueError("Structure could not be found.")

                        molecules.append(
                            Molecule(
                                temp_SMILES,
                                temp_synonyms,
                                temp_CAS,
                                result["id"],
                                mode,
                                "chemeo",
                            )
                        )

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def get_molecules_from_pubchem_batchmode(
        self, original_identifiers: list[str], mode: str, standardize: bool = True
    ) -> list[Optional[list[Molecule]]]:
        # https://pubchem.ncbi.nlm.nih.gov/docs/power-user-gateway
        # https://pubchem.ncbi.nlm.nih.gov/docs/identifier-exchange-service
        MoleculeResolver._check_parameters(
            identifiers=original_identifiers,
            modes=mode,
            services="pubchem",
            context="get_molecules_batch",
        )
        if not all([type(identifier) is str for identifier in original_identifiers]):
            raise TypeError("All identifiers must be strings.")

        def clean_identifier(identifier):
            if mode == "name":
                greek_letters = [
                    "α",
                    "β",
                    "γ",
                    "δ",
                    "ε",
                    "ϵ",
                    "ζ",
                    "η",
                    "θ",
                    "ι",
                    "κ",
                    "λ",
                    "μ",
                    "ν",
                    "ξ",
                    "ο",
                    "π",
                    "ρ",
                    "σ",
                    "τ",
                    "υ",
                    "φ",
                    "χ",
                    "ψ",
                    "ω",
                    "Α",
                    "Β",
                    "Γ",
                    "Δ",
                    "Ε",
                    "Ζ",
                    "Η",
                    "Θ",
                    "Ι",
                    "Κ",
                    "Λ",
                    "Μ",
                    "Ν",
                    "Ξ",
                    "Ο",
                    "Π",
                    "Ρ",
                    "Σ",
                    "Τ",
                    "Υ",
                    "Φ",
                    "Χ",
                    "Ψ",
                    "Ω",
                    "µ",
                    "∆",
                ]
                spelled_out_versions = [
                    "alpha",
                    "beta",
                    "gamma",
                    "delta",
                    "epsilon",
                    "epsilon",
                    "zeta",
                    "eta",
                    "theta",
                    "iota",
                    "kappa",
                    "lambda",
                    "mu",
                    "nu",
                    "xi",
                    "omicron",
                    "pi",
                    "rho",
                    "sigma",
                    "tau",
                    "upsilon",
                    "phi",
                    "chi",
                    "psi",
                    "omega",
                    "alpha",
                    "beta",
                    "gamma",
                    "delta",
                    "epsilon",
                    "zeta",
                    "eta",
                    "theta",
                    "iota",
                    "kappa",
                    "lambda",
                    "mu",
                    "nu",
                    "xi",
                    "omicron",
                    "pi",
                    "rho",
                    "sigma",
                    "tau",
                    "upsilon",
                    "phi",
                    "chi",
                    "psi",
                    "omega",
                    "mu",
                    "delta",
                ]

                map_to_replace = [
                    ("’", "'"),
                    ("′", "'"),
                    ("±", "+-"),
                    ("→", "-->"),
                    ("≥", ">="),
                    ("≤", "<="),
                    ("·", "."),
                    ("#", "no. "),
                    ("«", ""),
                    ("»", ""),
                ]

                for greek_letter, spelled_out_version in zip(
                    greek_letters, spelled_out_versions
                ):
                    map_to_replace.append((greek_letter, spelled_out_version))

                old = identifier
                identifier = html.unescape(identifier)
                for old, new in map_to_replace:
                    identifier = identifier.replace(old, new)

                nfkd_form = unicodedata.normalize("NFKD", identifier)
                identifier = "".join(
                    [c for c in nfkd_form if not unicodedata.combining(c)]
                )

                if identifier.endswith(", (+-)-"):
                    identifier = "".join(identifier.split(", (+-)-")[:-1])
                if identifier.startswith("(+-)-"):
                    identifier = "".join(identifier.split("(+-)-")[1:])

                return identifier
            else:
                return identifier

        PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
        pubchem_xml_mode = {
            "name": "synonyms",
            "cas": "synonyms",
            "formula": "synonyms",
            "smiles": "smiles",
            "inchi": "inchis",
            "inchikey": "inchi-keys",
        }

        def parse_request_id(request_text):
            matches = regex.findall(
                "PCT-Waiting_reqid>(.*)</PCT-Waiting_reqid", request_text
            )
            if len(matches) != 1:
                return None
            request_id = matches[0]
            return request_id

        def cid_search_request():
            query_Uids = []
            for cleaned_unique_identifier in cleaned_unique_identifiers:
                query_Uids.append(
                    f"<PCT-QueryUids_{pubchem_xml_mode[mode]}_E>{cleaned_unique_identifier}</PCT-QueryUids_{pubchem_xml_mode[mode]}_E>"
                )
            temp = "".join(query_Uids)
            root_query_Uids = f"<PCT-QueryUids_{pubchem_xml_mode[mode]}>{temp}</PCT-QueryUids_{pubchem_xml_mode[mode]}>"

            xml_template = f"""<?xml version="1.0"?>
                                <!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "NCBI_PCTools.dtd">
                                <PCT-Data>
                                <PCT-Data_input>
                                    <PCT-InputData>
                                    <PCT-InputData_query>
                                        <PCT-Query>
                                        <PCT-Query_type>
                                            <PCT-QueryType>
                                            <PCT-QueryType_id-exchange>
                                                <PCT-QueryIDExchange>
                                                <PCT-QueryIDExchange_input>
                                                    <PCT-QueryUids>
                                                    {root_query_Uids}
                                                    </PCT-QueryUids>
                                                </PCT-QueryIDExchange_input>
                                                <PCT-QueryIDExchange_operation-type value="same"/>
                                                <PCT-QueryIDExchange_output-type value="cid"/>
                                                <PCT-QueryIDExchange_output-method value="file-pair"/>
                                                <PCT-QueryIDExchange_compression value="gzip"/>
                                                </PCT-QueryIDExchange>
                                            </PCT-QueryType_id-exchange>
                                            </PCT-QueryType>
                                        </PCT-Query_type>
                                        </PCT-Query>
                                    </PCT-InputData_query>
                                    </PCT-InputData>
                                </PCT-Data_input>
                                </PCT-Data>"""

            request_id = None
            n_try = 0
            while request_id is None and n_try < 5:
                request_text = self._resilient_request(
                    PUBCHEM_URL,
                    {
                        "data": xml_template.encode("utf-8"),
                        "headers": {"Content-type": "application/xml; charset=utf-8"},
                    },
                    request_type="post",
                )

                if (
                    request_text is not None
                    and request_text.count("Result set is empty.") == 0
                    and request_text.count("server-error") == 0
                ):
                    request_id = parse_request_id(request_text)

                n_try += 1
                time.sleep(3)

            return request_id

        def info_request(cids_to_request, info_name):
            xml_for_cids_to_request = [
                f"<PCT-ID-List_uids_E>{cid}</PCT-ID-List_uids_E>"
                for cid in cids_to_request
            ]
            xml_template = f"""<?xml version="1.0"?>
                            <!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "NCBI_PCTools.dtd">
                            <PCT-Data>
                            <PCT-Data_input>
                                <PCT-InputData>
                                <PCT-InputData_query>
                                    <PCT-Query>
                                    <PCT-Query_type>
                                        <PCT-QueryType>
                                        <PCT-QueryType_id-exchange>
                                            <PCT-QueryIDExchange>
                                            <PCT-QueryIDExchange_input>
                                                <PCT-QueryUids>
                                                <PCT-QueryUids_ids>
                                                    <PCT-ID-List>
                                                    <PCT-ID-List_db>pccompound</PCT-ID-List_db>
                                                    <PCT-ID-List_uids>
                                                        {''.join(xml_for_cids_to_request)}
                                                    </PCT-ID-List_uids>
                                                    </PCT-ID-List>
                                                </PCT-QueryUids_ids>
                                                </PCT-QueryUids>
                                            </PCT-QueryIDExchange_input>
                                            <PCT-QueryIDExchange_operation-type value="same"/>
                                            <PCT-QueryIDExchange_output-type value="{info_name}"/>
                                            <PCT-QueryIDExchange_output-method value="file-pair"/>
                                            <PCT-QueryIDExchange_compression value="gzip"/>
                                            </PCT-QueryIDExchange>
                                        </PCT-QueryType_id-exchange>
                                        </PCT-QueryType>
                                    </PCT-Query_type>
                                    </PCT-Query>
                                </PCT-InputData_query>
                                </PCT-InputData>
                            </PCT-Data_input>
                            </PCT-Data>"""

            request_id = None
            n_try = 0
            while request_id is None or n_try < 5:
                request_text = self._resilient_request(
                    PUBCHEM_URL,
                    {
                        "data": xml_template.encode("utf-8"),
                        "headers": {"Content-type": "application/xml; charset=utf-8"},
                    },
                    request_type="post",
                )

                if (
                    request_text is not None
                    and request_text.count("Result set is empty.") == 0
                    and request_text.count("server-error") == 0
                ):
                    request_id = parse_request_id(request_text)

                n_try += 1
                time.sleep(3)

            return request_id

        def poll_request(_request_id):
            download_url = None
            poll_request_xml = f"""<?xml version="1.0"?>
                        <!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">
                        <PCT-Data>
                        <PCT-Data_input>
                            <PCT-InputData>
                            <PCT-InputData_request>
                                <PCT-Request>
                                <PCT-Request_reqid>{_request_id}</PCT-Request_reqid>
                                <PCT-Request_type value="status"/>
                                </PCT-Request>
                            </PCT-InputData_request>
                            </PCT-InputData>
                        </PCT-Data_input>
                        </PCT-Data>"""

            n_try = 0
            while True:
                poll_response = self._resilient_request(
                    PUBCHEM_URL,
                    {
                        "data": poll_request_xml.encode("utf-8"),
                        "headers": {
                            "Content-type": "application/xml; charset=utf-8",
                            "timeout": str(30),
                        },
                    },
                    request_type="post",
                )
                if poll_response:
                    matches = regex.findall(
                        "PCT-Download-URL_url>(.*)</PCT-Download-URL_url", poll_response
                    )

                    if matches:
                        if len(matches) != 1:
                            raise
                        download_url = matches[0]
                        if download_url.startswith("ftp://"):
                            download_url = "https" + download_url[3:]
                            break

                if n_try > 15:
                    break

                time.sleep(2)

            return download_url

        def get_results(download_url):
            found_results = {}
            n_try = 5
            while True:
                try:
                    with closing(urllib.request.urlopen(download_url)) as r:
                        zipped_content = gzip.GzipFile(fileobj=r)
                        content = zipped_content.read()
                        content_str = content.decode("utf-8")
                        results = content_str.split("\n")
                        found_results = [
                            result
                            for result in results
                            if result != "" and result != "Result set is empty."
                        ]
                        break
                except Exception as e:
                    if n_try > 5:
                        raise e
                    n_try += 1
                    time.sleep(3)

            return found_results

        with self.query_molecule_cache_batchmode(
            "pubchem", mode, original_identifiers, save_not_found=False
        ) as (
            original_identifiers_to_search,
            indices_of_identifiers_to_search,
            results,
        ):
            if len(original_identifiers_to_search) == 0:
                return results

            cleaned_identifiers = [
                clean_identifier(identifier)
                for identifier in original_identifiers_to_search
            ]  # cleaned identifiers
            cleaned_unique_identifiers = set(cleaned_identifiers)

            request_id = cid_search_request()

            if request_id is not None:
                download_url = poll_request(request_id)

                found_results_by_identifier = {}
                found_results_by_cid_identifier = {}
                found_results_by_cid = {}
                SMILES = {}
                CAS = {}

                results_with_name_errors = []
                for result in get_results(download_url):
                    cleaned_unique_identifier, cid = result.split("\t")

                    if cleaned_unique_identifier not in found_results_by_cid_identifier:
                        found_results_by_cid_identifier[cleaned_unique_identifier] = []

                    if cid != "":
                        found_results_by_cid_identifier[
                            cleaned_unique_identifier
                        ].append(int(cid))
                    else:
                        found_results_by_cid_identifier[
                            cleaned_unique_identifier
                        ].append(None)

                    if "#" in cleaned_unique_identifier:
                        results_with_name_errors.append(cleaned_unique_identifier)

                if results_with_name_errors:
                    temp = {}
                    for name_with_error in results_with_name_errors:
                        parts = name_with_error.split("#")
                        longest_part = max(parts, key=len)
                        candidates = [
                            v for v in cleaned_identifiers if longest_part in v
                        ]
                        temp[name_with_error] = candidates
                    print(
                        "The following identifiers have errors: {}".format(
                            ", ".join(temp)
                        )
                    )

                if len(found_results_by_cid_identifier) > 0:
                    for original_identifier, cleaned_identifier in zip(
                        original_identifiers, cleaned_identifiers
                    ):
                        if cleaned_identifier in found_results_by_cid_identifier:
                            cids = found_results_by_cid_identifier[cleaned_identifier]
                            cids = [cid for cid in cids if cid is not None]
                            if len(cids) == 0:
                                found_results_by_identifier[original_identifier] = None
                            else:
                                cid = cids[
                                    0
                                ]  # first cid seems always to be the best match
                                found_results_by_identifier[original_identifier] = cid
                                found_results_by_cid[cid] = original_identifier

                    cids_to_request = found_results_by_cid.keys()

                    synonyms = {}
                    for cid in cids_to_request:
                        if cid not in synonyms:
                            synonyms[cid] = []
                            CAS[cid] = []

                        if mode == "name":
                            this_cid_requested_name = found_results_by_cid[cid]
                            synonyms[cid].append(this_cid_requested_name)

                    request_id = info_request(cids_to_request, "iupac")
                    download_url = poll_request(request_id)
                    iupac_results = get_results(download_url)
                    if not len(cids_to_request) <= len(iupac_results):
                        raise

                    for result in iupac_results:
                        cid, iupac_name = result.split("\t")

                        if len(iupac_name) > 1:
                            synonyms[int(cid)].append(iupac_name)

                    request_id = info_request(cids_to_request, "synonyms")
                    download_url = poll_request(request_id)
                    synonym_results = get_results(download_url)

                    for result in synonym_results:
                        cid, synonym = result.split("\t")
                        synonyms[int(cid)].append(synonym)

                    request_id = info_request(cids_to_request, "smiles")
                    download_url = poll_request(request_id)
                    SMILES_results = get_results(download_url)
                    if not len(cids_to_request) <= len(SMILES_results):
                        raise

                    for result in SMILES_results:
                        cid, this_SMILES = result.split("\t")
                        this_SMILES = MoleculeResolver.standardize_SMILES(
                            this_SMILES, standardize
                        )
                        SMILES[int(cid)] = this_SMILES

                for molecule_index, original_identifier in zip(
                    indices_of_identifiers_to_search,
                    original_identifiers_to_search,
                ):
                    cid = (
                        found_results_by_identifier[original_identifier]
                        if original_identifier in found_results_by_identifier
                        else None
                    )
                    if cid is not None:
                        results[molecule_index] = [
                            Molecule(
                                SMILES[cid],
                                MoleculeResolver.filter_and_sort_synonyms(
                                    synonyms[cid]
                                ),
                                MoleculeResolver.filter_and_sort_CAS(synonyms[cid]),
                                cid,
                                mode,
                                "pubchem",
                                identifier=original_identifier,
                            )
                        ]

            return results

    def get_molecule_from_pubchem(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        if required_formula is None:
            if mode == "formula":
                required_formula = identifier

        MoleculeResolver._check_parameters(
            identifiers=identifier,
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            modes=mode,
            services="pubchem",
        )

        with self.query_molecule_cache("pubchem", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                SMILES = None
                synonyms = []
                CAS = []
                cid = None

                PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
                pubchem_mode = mode
                if mode == "cas":
                    pubchem_mode = "name"

                # here we use rejected_status_codes = [400, 404] 400 means the request was not correctly formatted, but
                # we can rule this out and it is probably to the identifier causing some URL issues. Because of this it is rejected.
                accepted_status_codes = (
                    [200, 202] if pubchem_mode == "formula" else [200]
                )
                request_text = self._resilient_request(
                    f'{PUBCHEM_URL}{pubchem_mode}/json?{pubchem_mode}={urllib.parse.quote(identifier, safe="")}',
                    rejected_status_codes=[400, 404],
                    accepted_status_codes=accepted_status_codes,
                )

                accepted_results = []

                if request_text is not None:
                    temp = json.loads(request_text)

                    if pubchem_mode == "formula":
                        listkey = temp["Waiting"]["ListKey"]
                        cids_response_text = self._resilient_request(
                            f"{PUBCHEM_URL}listkey/{listkey}/cids/json",
                            rejected_status_codes=[400, 404],
                            accepted_status_codes=[200],
                            sleep_time=7,
                        )

                        if cids_response_text is not None:
                            cid_list = json.loads(cids_response_text)["IdentifierList"][
                                "CID"
                            ]

                            for cid_temp in cid_list:
                                cmp_temp = self.get_molecule_from_pubchem(
                                    str(cid_temp),
                                    "cid",
                                    required_formula,
                                    required_charge,
                                    required_structure_type,
                                    standardize,
                                )
                                if cmp_temp:
                                    cmp_temp.mode = mode
                                    accepted_results.append(cmp_temp)

                        if not accepted_results:
                            # if nothing was foung, try with the hill formula, as pubchem sometimes finds results using it
                            hill_formula = MoleculeResolver.to_hill_formula(
                                MoleculeResolver.formula_to_dictionary(identifier)
                            )
                            if identifier != hill_formula:
                                return self.get_molecule_from_pubchem(
                                    hill_formula,
                                    mode,
                                    required_formula,
                                    required_charge,
                                    required_structure_type,
                                    standardize,
                                )

                    else:

                        def get_prop_value(
                            _compound, label, name, conversion_funtion, all_names=False
                        ):
                            props_found = [
                                prop["value"]
                                for prop in _compound["props"]
                                if prop["urn"]["label"] == label
                                and prop["urn"]["name"] == name
                            ]

                            if all_names:
                                return {
                                    prop["urn"]["name"]: conversion_funtion(
                                        prop["value"]["sval"]
                                    )
                                    for prop in _compound["props"]
                                    if prop["urn"]["label"] == label
                                }
                            else:
                                if len(props_found) != 1:
                                    raise
                                prop_vals = list(props_found[0].values())

                                return conversion_funtion(prop_vals[0])

                        # for i in range(len(temp['PC_Compounds'])):
                        # first cid seems always to be the best match
                        compound = temp["PC_Compounds"][0]

                        if len(compound["id"]) == 0:
                            molecules.clear()
                            return

                        cid = compound["id"]["id"]["cid"]
                        SMILES = MoleculeResolver.standardize_SMILES(
                            get_prop_value(compound, "SMILES", "Isomeric", str),
                            standardize,
                        )
                        SMILES_from_InChI = MoleculeResolver.InChI_to_SMILES(
                            get_prop_value(compound, "InChI", "Standard", str),
                            standardize,
                        )

                        if MoleculeResolver.check_SMILES(
                            SMILES,
                            required_formula,
                            required_charge,
                            required_structure_type,
                        ):
                            pass
                        elif MoleculeResolver.check_SMILES(
                            SMILES_from_InChI,
                            required_formula,
                            required_charge,
                            required_structure_type,
                        ):
                            SMILES = SMILES_from_InChI
                        else:
                            molecules.clear()
                            return

                        synonym_response_text = self._resilient_request(
                            f"{PUBCHEM_URL}cid/{cid}/synonyms/TXT"
                        )
                        if synonym_response_text is None:
                            temp_synonyms = []
                        else:
                            temp_synonyms = synonym_response_text.split("\n")
                        IUPAC_names = get_prop_value(
                            compound, "IUPAC Name", "Preferred", str, all_names=True
                        )

                        if "Preferred" in IUPAC_names:
                            temp_synonyms.insert(0, IUPAC_names["Preferred"])
                        elif "CAS-like Style" in IUPAC_names:
                            temp_synonyms.insert(0, IUPAC_names["CAS-like Style"])
                        elif "Systematic" in IUPAC_names:
                            temp_synonyms.insert(0, IUPAC_names["Systematic"])
                        elif "Traditional" in IUPAC_names:
                            temp_synonyms.insert(0, IUPAC_names["Traditional"])

                        CAS = MoleculeResolver.filter_and_sort_CAS(temp_synonyms)
                        synonyms = MoleculeResolver.filter_and_sort_synonyms(
                            temp_synonyms
                        )

                        molecules.append(
                            Molecule(SMILES, synonyms, CAS, cid, mode, "pubchem")
                        )

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def get_molecule_from_CAS_registry(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            modes=mode,
            services="cas_registry",
        )

        with self.query_molecule_cache("cas_registry", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                CAS_URL = "https://commonchemistry.cas.org/api/"
                SMILES = None
                synonyms = []
                CAS = []
                mode_used = mode

                search_response_text = self._resilient_request(
                    f'{CAS_URL}search?q={urllib.parse.quote(identifier, safe="")}',
                    {"headers": {"accept": "application/json"}},
                    rejected_status_codes=[403, 404],
                )

                if search_response_text is not None:
                    results = json.loads(search_response_text)
                    results = results["results"]

                    for result in results:
                        CAS = result["rn"]
                        detail_response_text = self._resilient_request(
                            f"{CAS_URL}detail?cas_rn={urllib.parse.quote(CAS)}",
                            {"headers": {"accept": "application/json"}},
                        )
                        if detail_response_text is not None:
                            details = json.loads(detail_response_text)

                            synonyms = []
                            if "name" in details:
                                synonyms.append(details["name"])
                            synonyms.extend(details["synonyms"])

                            SMILES = details["smile"]
                            SMILES = MoleculeResolver.standardize_SMILES(
                                SMILES, standardize
                            )

                            inchi = details["inchi"]

                            if MoleculeResolver.check_SMILES(
                                SMILES,
                                required_formula,
                                required_charge,
                                required_structure_type,
                            ):
                                molecules.append(
                                    Molecule(
                                        SMILES,
                                        MoleculeResolver.filter_and_sort_synonyms(
                                            synonyms
                                        ),
                                        [CAS],
                                        None,
                                        mode_used,
                                        "cas_registry",
                                    )
                                )
                            elif inchi != "":
                                SMILES_from_InChI = MoleculeResolver.InChI_to_SMILES(
                                    inchi, standardize
                                )
                                if MoleculeResolver.check_SMILES(
                                    SMILES_from_InChI,
                                    required_formula,
                                    required_charge,
                                    required_structure_type,
                                ):
                                    molecules.append(
                                        Molecule(
                                            SMILES_from_InChI,
                                            MoleculeResolver.filter_and_sort_synonyms(
                                                synonyms
                                            ),
                                            [CAS],
                                            None,
                                            mode_used,
                                            "cas_registry",
                                        )
                                    )

                # the search by SMILES on the CAS common chemistry service is not doing any kind of canonicalization
                # prior to searching, so in a lot of cases it does not deliver a result even though the molecule is
                # in the database, the service however does not suffer from this issue when searching by InChI
                if mode == "smiles" and len(molecules) == 0:
                    inchi = MoleculeResolver.SMILES_to_InChI(identifier)
                    cmp = self.get_molecule_from_CAS_registry(
                        inchi,
                        "inchi",
                        required_formula,
                        required_charge,
                        required_structure_type,
                        standardize,
                    )
                    if cmp is not None:
                        cmp.mode = "inchi calculated from smiles"
                    molecules.append(cmp)

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def _match_SRS_results_to_identifiers(
        self, identifiers: list[str], mode: str, standardize: bool, results: list[dict]
    ):
        infos_by_ITN = {}
        ITNs_by_primary_name = {}
        ITNs_by_synonym = {}
        ITNs_by_all_names = {}
        ITNs_by_CAS = {}
        for result in results:
            primary_names = []
            epaName = result["epaName"]
            if epaName:
                primary_names.append(epaName.strip())
            iupacName = result["iupacName"]
            if iupacName:
                primary_names.append(iupacName.strip())
            systematicName = result["systematicName"]
            if systematicName:
                primary_names.append(systematicName.strip())

            CAS = []
            if "currentCasNumber" in result:
                if result["currentCasNumber"]:
                    CAS.append(result["currentCasNumber"].strip())

            synonyms = []
            if "synonyms" in result:
                synonyms.extend(
                    [
                        synonym["synonymName"].strip()
                        for synonym in result["synonyms"]
                        if synonym["synonymName"]
                    ]
                )

            ITN = result["internalTrackingNumber"]
            SMILES = result["smilesNotation"]
            SMILES = self.standardize_SMILES(SMILES, standardize)

            if SMILES is None:
                if isinstance(result["inchiNotation"], str):
                    inchi = result["inchiNotation"]
                    if not inchi.lower().startswith("inchi="):
                        inchi = "InChI=" + inchi
                    SMILES = MoleculeResolver.InChI_to_SMILES(inchi, standardize)

            if not SMILES:
                continue

            if ITN in infos_by_ITN:
                raise RuntimeError(
                    "ITN already exists in infos_by_ITN, ITNs are expected to be returned only once."
                )

            all_synonyms_lower = []
            for name in primary_names:
                nl = name.lower()
                all_synonyms_lower.append(nl)
                if nl not in ITNs_by_primary_name:
                    ITNs_by_primary_name[nl] = []

                if ITN not in ITNs_by_primary_name[nl]:
                    ITNs_by_primary_name[nl].append(ITN)

                if nl not in ITNs_by_all_names:
                    ITNs_by_all_names[nl] = []

                if ITN not in ITNs_by_all_names[nl]:
                    ITNs_by_all_names[nl].append(ITN)

            for name in synonyms:
                nl = name.lower()
                all_synonyms_lower.append(nl)
                if nl not in ITNs_by_synonym:
                    ITNs_by_synonym[nl] = []

                if ITN not in ITNs_by_synonym[nl]:
                    ITNs_by_synonym[nl].append(ITN)

                if nl not in ITNs_by_all_names:
                    ITNs_by_all_names[nl] = []

                if ITN not in ITNs_by_all_names[nl]:
                    ITNs_by_all_names[nl].append(ITN)

            for CAS_ in CAS:
                if CAS_:
                    if CAS_ not in ITNs_by_CAS:
                        ITNs_by_CAS[CAS_] = []

                    if ITN not in ITNs_by_CAS[CAS_]:
                        ITNs_by_CAS[CAS_].append(ITN)

            infos_by_ITN[ITN] = (
                SMILES,
                primary_names,
                synonyms,
                CAS,
                ITN,
                all_synonyms_lower,
            )

        SRS_results_by_identifier = {}
        uniquely_matched_ITNs = {}

        last_SRS_results_by_identifier_length = -1
        i_iteration = 0
        max_iterations = len(results) * 2
        while len(SRS_results_by_identifier) != last_SRS_results_by_identifier_length:
            i_iteration += 1
            if i_iteration > max_iterations:
                break
            last_SRS_results_by_identifier_length = len(SRS_results_by_identifier)
            for identifier in identifiers:
                il = identifier.lower()
                if identifier not in SRS_results_by_identifier:
                    ITNs_for_this_identifier = []
                    if mode == "cas":
                        if identifier in ITNs_by_CAS:
                            ITNs_for_this_identifier = ITNs_by_CAS[identifier]

                    elif mode == "name":
                        if il in ITNs_by_primary_name:
                            ITNs_for_this_identifier = ITNs_by_primary_name[il]
                        elif il in ITNs_by_synonym:
                            ITNs_for_this_identifier = ITNs_by_synonym[il]

                    best_ITNs = None

                    if ITNs_for_this_identifier:
                        if len(ITNs_for_this_identifier) == 1:
                            best_ITNs = ITNs_for_this_identifier
                        else:
                            temptative_ITNs = set(ITNs_for_this_identifier) - set(
                                uniquely_matched_ITNs.values()
                            )
                            if len(temptative_ITNs) == 1:
                                best_ITNs = list(temptative_ITNs)
                            else:
                                unique_SMILES = set(
                                    [
                                        self.standardize_SMILES(
                                            infos_by_ITN[ITN][0], True
                                        )
                                        for ITN in temptative_ITNs
                                    ]
                                )
                                if len(unique_SMILES) == 1:
                                    # if same structure, use the one with the smallest ITN
                                    best_ITNs = [
                                        str(min([int(ITN) for ITN in temptative_ITNs]))
                                    ]
                                else:
                                    if len(temptative_ITNs) > 1:
                                        number_of_synonyms_found_in_ITN_info = []
                                        for ITN in temptative_ITNs:
                                            number_of_synonyms_found_in_ITN_info.extend(
                                                infos_by_ITN[ITN][5].count(il) * [ITN]
                                            )

                                        most_common = collections.Counter(
                                            number_of_synonyms_found_in_ITN_info
                                        ).most_common(2)
                                        if most_common[0][1] > most_common[1][1]:
                                            best_ITNs = [most_common[0][0]]

                    if best_ITNs:
                        if len(best_ITNs) == 1:
                            uniquely_matched_ITNs[identifier] = best_ITNs[0]

                        SRS_results_by_identifier[identifier] = [
                            infos_by_ITN[ITN] for ITN in best_ITNs
                        ]

        return SRS_results_by_identifier

    def get_molecules_from_SRS_batchmode(
        self, identifiers: list[str], mode: str, standardize: bool = True
    ) -> list[Optional[Molecule]]:
        # https://www.postman.com/api-evangelist/workspace/environmental-protection-agency-epa/collection/35240-6b84cc71-ce77-48b8-babd-323eb8d670bd
        # new api https://cdxappstest.epacdx.net/oms-substance-registry-services/swagger-ui/

        MoleculeResolver._check_parameters(
            identifiers=identifiers,
            modes=mode,
            services="srs",
            context="get_molecules_batch",
        )
        if not all([type(identifier) is str for identifier in identifiers]):
            raise TypeError("All identifiers must be strings.")

        with self.query_molecule_cache_batchmode(
            "srs", mode, identifiers, save_not_found=False
        ) as (identifiers_to_search, indices_of_identifiers_to_search, results):
            if len(identifiers_to_search) == 0:
                return results

            SRS_URL = "https://cdxappstest.epacdx.net/oms-substance-registry-services/rest-api/substances"
            chunks_identifiers = []
            chunks_identifer_indices = []
            this_chunk_identifiers = []
            this_chunk_identifier_indices = []
            for identifier_index, identifier in zip(
                indices_of_identifiers_to_search, identifiers_to_search
            ):
                this_chunk_identifiers.append(identifier)
                this_chunk_identifier_indices.append(identifier_index)
                if (
                    len(
                        f'{SRS_URL}/{mode}?{mode}List={urllib.parse.quote("|".join(this_chunk_identifiers))}&qualifier=exact'
                    )
                    >= 2000
                ):
                    chunks_identifiers.append(this_chunk_identifiers)
                    chunks_identifer_indices.append(this_chunk_identifier_indices)
                    this_chunk_identifiers = []
                    this_chunk_identifier_indices = []

            if this_chunk_identifiers:
                chunks_identifiers.append(this_chunk_identifiers)
                chunks_identifer_indices.append(this_chunk_identifier_indices)

            for chunk_identifier_indices, chunk_identifiers in zip(
                chunks_identifer_indices, chunks_identifiers
            ):
                search_response_text = self._resilient_request(
                    f'{SRS_URL}/{mode}?{mode}List={urllib.parse.quote("|".join(chunk_identifiers))}&qualifier=exact',
                    rejected_status_codes=[404, 500],
                    kwargs={"timeout": 60},
                )

                if search_response_text is not None:
                    infos_by_identifier = self._match_SRS_results_to_identifiers(
                        chunk_identifiers,
                        mode,
                        standardize,
                        json.loads(search_response_text),
                    )

                    for molecule_index, identifier in zip(
                        chunk_identifier_indices, chunk_identifiers
                    ):
                        if identifier in infos_by_identifier:
                            this_molecules = []
                            for (
                                SMILES,
                                primary_names,
                                synonyms,
                                CAS,
                                ITN,
                                _,
                            ) in infos_by_identifier[identifier]:
                                temp_synonyms = (
                                    MoleculeResolver.filter_and_sort_synonyms(
                                        primary_names + synonyms
                                    )
                                )
                                this_molecules.append(
                                    Molecule(
                                        SMILES,
                                        temp_synonyms,
                                        CAS,
                                        ITN,
                                        mode,
                                        "srs",
                                        identifier=identifier,
                                    )
                                )

                            results[molecule_index] = this_molecules

            return results

    def get_molecule_from_SRS(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            services="srs",
            modes=mode,
        )

        with self.query_molecule_cache("srs", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                SRS_URL = "https://cdxappstest.epacdx.net/oms-substance-registry-services/rest-api/substance"
                search_response_text = self._resilient_request(
                    f'{SRS_URL}/{mode}/{urllib.parse.quote(identifier, safe="")}',
                    rejected_status_codes=[400, 404, 500],
                    kwargs={"timeout": 10},
                )

                if search_response_text is not None:
                    results = json.loads(search_response_text)
                    infos_by_identifier = self._match_SRS_results_to_identifiers(
                        [identifier], mode, standardize, results
                    )

                    if identifier in infos_by_identifier:
                        if len(infos_by_identifier[identifier]) != 1:
                            raise RuntimeError(
                                "More than one molecule found for the given identifier."
                            )
                        (
                            SMILES,
                            primary_names,
                            synonyms,
                            CAS,
                            ITN,
                            _,
                        ) = infos_by_identifier[identifier][0]
                        synonyms = MoleculeResolver.filter_and_sort_synonyms(
                            primary_names + synonyms
                        )
                        if MoleculeResolver.check_SMILES(
                            SMILES,
                            required_formula,
                            required_charge,
                            required_structure_type,
                        ):
                            molecules.append(
                                Molecule(
                                    SMILES,
                                    synonyms,
                                    CAS,
                                    ITN,
                                    mode,
                                    "srs",
                                    identifier=identifier,
                                )
                            )

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    @cache
    def _get_info_from_CIR(
        self,
        structure_identifier,
        representation,
        resolvers_to_use: Tuple[str],
        expected_number_of_results: Optional[int] = None,
    ) -> Optional[list[str]]:
        # got API info from:
        # https://cactus.nci.nih.gov/chemical/structure_documentation
        # https://search.r-project.org/CRAN/refmans/webchem/html/cir_query.html
        # https://github.com/mcs07/CIRpy

        if "CIR_is_down" in self._message_slugs_shown:
            return None

        CIR_URL = "https://cactus.nci.nih.gov/chemical/structure/"
        # info can be e.g. smiles, iupac_name
        try:
            resolver_info = ""
            if resolvers_to_use:
                resolver_info = f'?resolver={",".join(resolvers_to_use)}'
            # although the documentation says otherwise it returns a 500 response even if it should return 404
            response_text = self._resilient_request(
                f"{CIR_URL}{urllib.parse.quote(structure_identifier)}/{representation}{resolver_info}",
                rejected_status_codes=[404, 500],
            )
            if not response_text:
                return None
            response_values = response_text.split("\n")
            if response_values:
                if expected_number_of_results is not None:
                    if len(response_values) != expected_number_of_results:
                        raise RuntimeError(
                            "More than one iupac_name was found. This was unexpected."
                        )
                return response_values
        except requests.exceptions.ConnectionError:
            # I don't know why, but somtimes CIR is offline. This would make the module much slower as
            # it tries to connect multiple times anyway. Instead we give a warning and skip CIR.
            if "CIR_is_down" not in self._message_slugs_shown:
                self._message_slugs_shown.append("CIR_is_down")
                warnings.warn(
                    "CIR seems to be down, to continue working this instance of MoleculeResolver will skip CIR."
                )

        return None

    @cache
    def get_iupac_name_from_CIR(self, SMILES: str):
        result = self._get_info_from_CIR(SMILES, "iupac_name", ("smiles",))
        return result[0] if result else None

    def get_molecule_from_CIR(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        if required_formula is None:
            if mode == "formula":
                required_formula = identifier

        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            services="cir",
            modes=mode,
        )

        with self.query_molecule_cache("cir", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                resolvers_by_mode = {
                    "formula": None,
                    "name": ("name_by_cir",),
                    "smiles": ("smiles",),
                    "inchi": ("stdinchi",),
                    "inchikey": ("stdinchikey",),
                    "cas": ("cas_number",),
                }

                SMILES = self._get_info_from_CIR(
                    identifier, "smiles", resolvers_by_mode[mode], 1
                )
                if not SMILES:
                    return None
                else:
                    SMILES = MoleculeResolver.standardize_SMILES(SMILES[0], standardize)

                if SMILES:
                    CIR_names = self._get_info_from_CIR(
                        identifier, "names", resolvers_by_mode[mode]
                    )
                    synonyms = MoleculeResolver.filter_and_sort_synonyms(
                        CIR_names if CIR_names else []
                    )
                    CAS = MoleculeResolver.filter_and_sort_CAS(
                        CIR_names if CIR_names else []
                    )
                    molecules.append(
                        Molecule(SMILES, synonyms, CAS, mode=mode, service="cir")
                    )

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def get_molecule_from_NIST(
        self,
        identifier: str,
        mode: str,
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        if required_formula is None:
            if mode == "formula":
                required_formula = identifier

        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            services="nist",
            modes=mode,
        )

        with self.query_molecule_cache("nist", mode, identifier) as (
            entry_available,
            molecules,
        ):
            if not entry_available:
                nist_modes = {
                    "formula": "Formula",
                    "name": "Name",
                    "cas": "ID",
                    "inchi": "InChI",
                    "smiles": "InChI",
                }
                NIST_Webbook_ID_regex = (
                    r'<a\s+href\s*=\s*"/cgi/cbook\.cgi\?ID=(.\d+).*">'
                )

                mode_used = mode
                if mode == "smiles":
                    identifier = MoleculeResolver.SMILES_to_InChI(
                        identifier, standardize
                    )
                    mode_used = "inchi calculated from smiles"

                response_text = self._resilient_request(
                    f'https://webbook.nist.gov/cgi/cbook.cgi?{urllib.parse.quote(nist_modes[mode])}={urllib.parse.quote(identifier, safe="")}'
                )

                def parse_molecule(temp_content):
                    items = MoleculeResolver.parse_items_from_html(
                        temp_content,
                        None,
                        [
                            ("html", '<h1 id="Top">(.*?)</h1>', [1]),
                            ("html", NIST_Webbook_ID_regex, [0, 1]),
                            (
                                "text",
                                rf"InChI: ({MoleculeResolver.InChI_regex_compiled.pattern[1:-1]})",
                                [0, 1],
                            ),
                            (
                                "text",
                                f"CAS Registry Number: {MoleculeResolver.CAS_regex}",
                                [0, 1],
                            ),
                            (
                                "html",
                                r"<li>\s*<strong>\s*Other\s+names:\s*</strong>(.*?)<\/li>",
                                [0, 1],
                            ),
                        ],
                        [0, 1],
                    )
                    if len(items) == 0:
                        return
                    if len(items) != 1:
                        raise RuntimeError("The webpage probably changed the format.")
                    name, additional_informationm, inchi, CAS, synonyms = items[0]

                    if not inchi:
                        return

                    if not synonyms:
                        synonyms = ""

                    if CAS:
                        CAS = [CAS]
                    else:
                        CAS = []

                    synonyms = [
                        name.strip()
                        for name in html.unescape(synonyms).split(";")
                        if name.strip()
                    ]
                    synonyms.insert(0, name)
                    synonyms = MoleculeResolver.filter_and_sort_synonyms(synonyms)
                    CAS = MoleculeResolver.filter_and_sort_CAS(CAS)
                    SMILES = MoleculeResolver.InChI_to_SMILES(inchi, standardize)

                    if MoleculeResolver.check_SMILES(
                        SMILES,
                        required_formula,
                        required_charge,
                        required_structure_type,
                    ):
                        molecules.append(
                            Molecule(
                                SMILES,
                                synonyms,
                                CAS,
                                additional_informationm,
                                mode_used,
                                "nist",
                                1,
                                identifier,
                            )
                        )

                if response_text is not None:
                    relevant_response_text = MoleculeResolver.normalize_html(
                        response_text
                    )
                    relevant_response_text = regex.findall(
                        '<main id="main">(.*)</main>', relevant_response_text
                    )
                    if len(relevant_response_text) != 1:
                        raise RuntimeError("The format of the page might have changed.")
                    relevant_response_text = relevant_response_text[0].strip()

                    if (
                        not regex.search(
                            "<h1>.*(Not Found|no)\s*</h1>", relevant_response_text
                        )
                        and not regex.search(
                            "<h1>.*No Matching Names Found\s*</h1>",
                            relevant_response_text,
                        )
                        and not regex.search(
                            "<h1>.*No Matching Species Found\s*</h1>",
                            relevant_response_text,
                        )
                        and "no matching entries" not in relevant_response_text
                    ):
                        if "Search Results" not in relevant_response_text:
                            parse_molecule(relevant_response_text)
                        else:
                            # check all search results using check_SMILES
                            ids = regex.findall(NIST_Webbook_ID_regex, response_text)

                            for id in ids:
                                temp_result = self.get_molecule_from_NIST(
                                    id,
                                    "cas",
                                    required_formula,
                                    required_charge,
                                    required_structure_type,
                                    standardize,
                                )
                                if temp_result is not None:
                                    temp_result.mode = mode_used
                                    molecules.append(temp_result)

        return MoleculeResolver.filter_and_combine_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

    def find_salt_molecules(
        self,
        identifiers: list[str],
        modes: list[str] = ["name"],
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        services_to_use: Optional[list[str]] = None,
        search_iupac_name: bool = False,
        interactive: bool = False,
        minimum_number_of_cross_checks: Optional[int] = 1,
        ignore_exceptions: bool = False,
    ) -> tuple[list, list[int]]:
        (
            flattened_identifiers,
            flattened_modes,
            synonyms,
            CAS,
            given_SMILES,
        ) = MoleculeResolver._check_and_flatten_identifiers_and_modes(
            identifiers, modes
        )
        CAS = list(CAS)

        if given_SMILES:
            if MoleculeResolver.get_structure_type_from_SMILES(given_SMILES) != "salt":
                raise ValueError("The given SMILES is not a salt.")

        if required_charge is None:
            required_charge = 0

        if required_structure_type is None:
            required_structure_type = "salt"

        if minimum_number_of_cross_checks is None:
            minimum_number_of_cross_checks = 1

        salt_info = self.find_single_molecule_cross_checked(
            identifiers,
            modes=modes,
            required_formula=required_formula,
            required_charge=required_charge,
            required_structure_type=required_structure_type,
            services_to_use=services_to_use,
            search_iupac_name=search_iupac_name,
            minimum_number_of_cross_checks=minimum_number_of_cross_checks,
            ignore_exceptions=ignore_exceptions,
        )

        if salt_info[0] is None:
            salt_info = (
                given_SMILES,
                synonyms,
                CAS,
                "given identifiers: " + ", ".join(flattened_identifiers),
                ", ".join(flattened_modes),
                "manual",
                1,
            )
            if interactive:
                salt_info = self.find_single_molecule(
                    identifiers,
                    modes=modes,
                    required_formula=required_formula,
                    required_charge=required_charge,
                    required_structure_type=required_structure_type,
                    services_to_use=[],
                    search_iupac_name=search_iupac_name,
                    interactive=interactive,
                    ignore_exceptions=ignore_exceptions,
                )

        all_molecules = []
        stoichometric_coefficients = []
        if salt_info[0] is not None:
            all_molecules.append(salt_info)
            stoichometric_coefficients.append(-1)
            SMILES = salt_info[0]

            if MoleculeResolver.is_valid_SMILES(SMILES):
                if SMILES.count(".") > 0:
                    ionic_SMILES_list = SMILES.split(".")
                    ionic_SMILES_list = [
                        MoleculeResolver.standardize_SMILES(smi, True)
                        for smi in ionic_SMILES_list
                    ]
                    SMILES = ".".join(ionic_SMILES_list)
                    ionic_SMILES_set = set(ionic_SMILES_list)
                    if not len(ionic_SMILES_set) > 1:
                        raise ValueError("SMILES does not represent a salt")

                charge_sum = 0
                for ionic_SMILES in ionic_SMILES_set:
                    stoichometric_coefficients.append(SMILES.count(ionic_SMILES))

                    ionic_mol = MoleculeResolver.get_from_SMILES(ionic_SMILES)
                    ionic_charge = Chem.rdmolops.GetFormalCharge(ionic_mol)

                    ionic_info = self.find_single_molecule_cross_checked(
                        ionic_SMILES,
                        modes="SMILES",
                        required_charge=ionic_charge,
                        required_structure_type="ion",
                        services_to_use=services_to_use,
                        search_iupac_name=search_iupac_name,
                        minimum_number_of_cross_checks=minimum_number_of_cross_checks,
                        ignore_exceptions=ignore_exceptions,
                    )

                    if ionic_info[0] is None:
                        ionic_info = (
                            ionic_SMILES,
                            [],
                            [],
                            f"given identifiers: {ionic_SMILES}",
                            "smiles",
                            "manual",
                            1,
                        )
                        if interactive:
                            ionic_info = self.find_single_molecule(
                                ionic_SMILES,
                                modes="SMILES",
                                required_charge=ionic_charge,
                                required_structure_type="ion",
                                services_to_use=[],
                                search_iupac_name=search_iupac_name,
                                interactive=interactive,
                                ignore_exceptions=ignore_exceptions,
                            )

                    if ionic_info[0] is not None:
                        all_molecules.append(ionic_info)

                        charge_sum += stoichometric_coefficients[-1] * ionic_charge

                if charge_sum != 0:
                    all_molecules = [salt_info]

        # not all molecules have been found
        if len(all_molecules) < 3:
            cation_info = None
            anion_info = None
            all_molecules = []
            if salt_info[0] is not None:
                all_molecules.append(salt_info)
                stoichometric_coefficients = [-1]

            for synonym in synonyms:
                synonym_parts = synonym.split(" ")
                possible_cation_name = synonym_parts[0]

                cation_info = self.find_single_molecule_cross_checked(
                    possible_cation_name,
                    modes=["name"],
                    required_charge="positive",
                    required_structure_type="ion",
                    services_to_use=services_to_use,
                    search_iupac_name=search_iupac_name,
                    minimum_number_of_cross_checks=minimum_number_of_cross_checks,
                    ignore_exceptions=ignore_exceptions,
                )

                if cation_info[0] is None:
                    if interactive:
                        cation_info = self.find_single_molecule(
                            possible_cation_name,
                            modes=["name"],
                            required_charge="positive",
                            required_structure_type="ion",
                            services_to_use=[],
                            search_iupac_name=search_iupac_name,
                            interactive=interactive,
                            ignore_exceptions=ignore_exceptions,
                        )

                possible_anion_name = synonym_parts[1:]

                anion_info = self.find_single_molecule_cross_checked(
                    possible_anion_name,
                    modes=["name"],
                    required_charge="negative",
                    required_structure_type="ion",
                    services_to_use=services_to_use,
                    search_iupac_name=search_iupac_name,
                    minimum_number_of_cross_checks=minimum_number_of_cross_checks,
                    ignore_exceptions=ignore_exceptions,
                )

                if anion_info[0] is None:
                    if interactive:
                        anion_info = self.find_single_molecule(
                            possible_anion_name,
                            modes=["name"],
                            required_charge="negative",
                            required_structure_type="ion",
                            services_to_use=[],
                            search_iupac_name=search_iupac_name,
                            interactive=interactive,
                            ignore_exceptions=ignore_exceptions,
                        )

                if cation_info[0] is not None and anion_info[0] is not None:
                    all_molecules.append(cation_info)
                    all_molecules.append(anion_info)

                if len(all_molecules) > 2:
                    break

            if len(all_molecules) > 2:
                if len(all_molecules) != 3:
                    raise NotImplementedError(
                        "Functionality for salts with more than 2 ions has not been implemented."
                    )

                cation_charge = Chem.rdmolops.GetFormalCharge(
                    MoleculeResolver.get_from_SMILES(cation_info[0])
                )
                anion_charge = Chem.rdmolops.GetFormalCharge(
                    MoleculeResolver.get_from_SMILES(anion_info[0])
                )

                if cation_charge >= abs(anion_charge):
                    if cation_charge % abs(anion_charge) == 0:
                        stoichometric_coefficients.append(1)
                        stoichometric_coefficients.append(
                            abs(int(cation_charge / anion_charge))
                        )
                    else:
                        stoichometric_coefficients.append(abs(anion_charge))
                        stoichometric_coefficients.append(abs(cation_charge))
                else:
                    if abs(anion_charge) % cation_charge == 0:
                        stoichometric_coefficients.append(
                            abs(int(anion_charge / cation_charge))
                        )
                        stoichometric_coefficients.append(1)
                    else:
                        stoichometric_coefficients.append(abs(anion_charge))
                        stoichometric_coefficients.append(abs(cation_charge))

                if stoichometric_coefficients[0] != -1:
                    temp_smiles = stoichometric_coefficients[1] * [cation_info[0]]
                    temp_smiles.extend(stoichometric_coefficients[2] * [anion_info[0]])
                    salt_info = (
                        MoleculeResolver.standardize_SMILES(
                            ".".join(temp_smiles), True
                        ),
                        synonyms,
                        list(CAS),
                        "from consisting ions",
                        "smiles",
                        "automatic",
                        1,
                    )
                    all_molecules.insert(0, salt_info)
                    stoichometric_coefficients.insert(0, -1)

        if len(all_molecules) < 3:
            stoichometric_coefficients = []

        return all_molecules, stoichometric_coefficients

    def find_single_molecule(
        self,
        identifiers: list[str],
        modes: list[str] = ["name"],
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        services_to_use: Optional[list[str]] = None,
        standardize: bool = True,
        search_iupac_name: bool = False,
        interactive: bool = False,
        ignore_exceptions: bool = False,
    ) -> Optional[Molecule]:
        if services_to_use is None:
            services_to_use = MoleculeResolver._available_services

        (
            flattened_identifiers,
            flattened_modes,
            synonyms,
            CAS,
            given_SMILES,
        ) = MoleculeResolver._check_and_flatten_identifiers_and_modes(
            identifiers, modes
        )
        MoleculeResolver._check_parameters(
            services=services_to_use,
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            context="find_single",
        )

        if required_formula is None and modes.count("formula") == 1:
            required_formula = identifiers[modes.index("formula")]

        SMILES = None
        additional_information = None
        mode_used = None
        identifier_used = None
        current_service = None
        try:
            for service in services_to_use:
                current_service = service
                if service == "cas_registry":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_CAS_registry(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS = set(
                                    cmp.CAS
                                )  # overwrite CAS with data from the CAS registry
                                additional_information = cmp.service
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break
                elif service == "pubchem":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_pubchem(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS.update(cmp.CAS)
                                additional_information = (
                                    f"{cmp.service} id: {cmp.additional_information}"
                                )
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break
                elif service == "cir":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_CIR(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                additional_information = cmp.service
                                mode_used = mode
                                identifier_used = cmp.identifier
                                break
                elif service == "opsin":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_OPSIN(
                                identifier,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                additional_information = cmp.service
                                mode_used = mode
                                identifier_used = cmp.identifier
                                break
                elif service == "chebi":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_ChEBI(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS.update(cmp.CAS)
                                additional_information = (
                                    f"{cmp.service} id: {cmp.additional_information}"
                                )
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break
                elif service == "srs":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_SRS(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS.update(cmp.CAS)
                                additional_information = (
                                    f"{cmp.service} id: {cmp.additional_information}"
                                )
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break
                elif service == "comptox":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_CompTox(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS.update(cmp.CAS)
                                additional_information = (
                                    f"{cmp.service} id: {cmp.additional_information}"
                                )
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break
                elif service == "chemeo":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_Chemeo(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS.update(cmp.CAS)
                                additional_information = (
                                    f"{cmp.service} id: {cmp.additional_information}"
                                )
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break
                elif service == "cts":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_CTS(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS.update(cmp.CAS)
                                additional_information = "cts"
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break
                elif service == "nist":
                    for identifier, mode in zip(flattened_identifiers, flattened_modes):
                        if (
                            mode
                            in MoleculeResolver._supported_modes_by_services[service]
                        ):
                            cmp = self.get_molecule_from_NIST(
                                identifier,
                                mode,
                                required_formula,
                                required_charge,
                                required_structure_type,
                                standardize,
                            )
                            if cmp is not None:
                                SMILES = cmp.SMILES
                                synonyms.extend(cmp.synonyms)
                                CAS.update(cmp.CAS)
                                additional_information = (
                                    f"{cmp.service} id: {cmp.additional_information}"
                                )
                                mode_used = cmp.mode
                                identifier_used = cmp.identifier
                                break

                if SMILES is not None:
                    break

            if SMILES is None:
                if given_SMILES is not None:
                    if MoleculeResolver.check_SMILES(
                        given_SMILES,
                        required_formula,
                        required_charge,
                        required_structure_type,
                    ):
                        SMILES = given_SMILES
                        additional_information = "given SMILES"
                        mode_used = "smiles"
                        if synonyms is None or len(synonyms) == 0:
                            search_iupac_name = True

            if SMILES is None:
                if len(synonyms) > 0 and required_charge is not None:
                    # if searching for an ion
                    # search for salts in pubchem and extract the single ions, take the one found most often
                    if required_charge != "zero" and required_charge != 0:
                        for identifier, mode in zip(
                            flattened_identifiers, flattened_modes
                        ):
                            molecules = (
                                self.get_molecule_for_ion_from_partial_pubchem_search(
                                    identifier,
                                    required_formula,
                                    required_charge,
                                    standardize,
                                )
                            )
                            if molecules is not None:
                                if len(molecules) > 0:
                                    mode_used = mode
                                    identifier_used = identifier
                                    current_service = "pubchem"
                                    additional_information = (
                                        "get_SMILES_for_ion_from_partial_pubchem_search"
                                    )
                                    SMILES = molecules[0][
                                        0
                                    ].SMILES  # get most likely candidate
                                    break

            if interactive and SMILES is None:
                return self.find_single_molecule_interactively(
                    identifiers,
                    modes,
                    required_formula=required_formula,
                    required_charge=required_charge,
                    standardize=standardize,
                )

            if SMILES is not None:
                # keep first synonym
                kept_synonyms = []
                if len(synonyms) > 0:
                    kept_synonyms = [synonyms.pop(0)]

                if search_iupac_name:
                    # add iupac name if available
                    iupac_name = self.get_iupac_name_from_CIR(SMILES)

                    if iupac_name is not None:
                        if not isinstance(iupac_name, list):
                            iupac_name = [iupac_name]
                        kept_synonyms.extend(iupac_name)

                # add rest of names
                if len(synonyms) > 0:
                    kept_synonyms.extend(synonyms)

                kept_synonyms = MoleculeResolver.filter_and_sort_synonyms(kept_synonyms)

                if len(kept_synonyms) == 0:
                    if interactive:
                        manual_name = ""
                        print("Could not infere name for following SMILES: " + SMILES)
                        while len(manual_name) < 3:
                            manual_name = input("please input manually:")
                            manual_name = manual_name.strip()

                        kept_synonyms = [manual_name]

                synonyms = kept_synonyms
        except Exception:
            if ignore_exceptions:
                print()
                print("Error with searching for", identifiers, "on", current_service)
                print(traceback.format_exc())
                print()
                return None
            else:
                raise

        if len(synonyms) == 0 or SMILES is None:
            return None

        SMILES = MoleculeResolver.standardize_SMILES(SMILES, standardize)

        return Molecule(
            SMILES,
            synonyms,
            list(CAS),
            additional_information,
            mode_used,
            current_service,
            1,
            identifier_used,
        )

    def find_single_molecule_interactively(
        self,
        identifiers: list[str],
        modes: list[str] = ["name"],
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        standardize: bool = True,
    ) -> Optional[Molecule]:
        (
            flattened_identifiers,
            flattened_modes,
            synonyms,
            CAS,
            given_SMILES,
        ) = MoleculeResolver._check_and_flatten_identifiers_and_modes(
            identifiers, modes
        )
        MoleculeResolver._check_parameters(
            required_formulas=required_formula,
            required_charges=required_charge,
            required_structure_types=required_structure_type,
            context="find_single",
        )

        SMILES = None
        start_geometry_comment = None

        # print header
        print("------------------------------------")
        formula = ""
        if required_formula is not None:
            formula = " | formula: " + required_formula
        charge = ""
        if required_charge is not None:
            charge = " | charge: " + str(required_charge)
        print(f"Searching for molecule {formula} {charge}")
        print("with the following identifiers:")
        for identifier, mode in zip(flattened_identifiers, flattened_modes):
            print(mode + ":", identifier)
        print("")

        mode = ""
        return_value = ""
        mode_used = "interactive"
        prompt_session = PromptSession()
        # regex for inchi, inchikey, SMILES, etc.
        # https://gist.github.com/lsauer/1312860/264ae813c2bd2c27a769d261c8c6b38da34e22fb
        while SMILES is None:
            if given_SMILES is not None:
                return_value = given_SMILES
            else:
                return_value = prompt_session.prompt(
                    "Please input the pubchem cid, name, CAS, SMILES, InChI or none to skip:",
                    default=str(return_value),
                )

            if return_value.lower() == "none":
                return None

            identifier = return_value
            try:
                return_value = int(return_value)
                mode = "cid"
                identifier = str(return_value)
            except Exception:
                if MoleculeResolver.is_valid_SMILES(return_value):
                    mode = "smiles"
                elif MoleculeResolver.is_valid_CAS(return_value):
                    mode = "cas"
                elif MoleculeResolver.is_valid_InChI(return_value):
                    mode = "inchi"
                else:
                    mode = "name"

            try:
                if mode == "smiles":
                    (
                        temptative_SMILES,
                        temptative_synonyms,
                        temptative_CAS,
                        temptative_start_geometry_comment,
                        mode_used,
                    ) = (identifier, synonyms, CAS, "manual geometry", "interactive")
                else:
                    cmp = self.find_single_molecule(
                        [identifier],
                        [mode],
                        required_formula,
                        required_charge,
                        standardize=standardize,
                    )
                    (
                        temptative_SMILES,
                        temptative_synonyms,
                        temptative_CAS,
                        temptative_start_geometry_comment,
                        mode_used,
                    ) = (
                        (
                            cmp.SMILES,
                            cmp.synonyms,
                            cmp.CAS,
                            cmp.additional_information,
                            cmp.mode,
                        )
                        if cmp is not None
                        else (None, [], [], "", "")
                    )
            except Exception:
                if mode == "cid":
                    print(
                        "There was a problem retreiving the molecule from pubchem. Does it exist?"
                    )

            if MoleculeResolver.check_SMILES(
                temptative_SMILES,
                required_formula,
                required_charge,
                required_structure_type,
            ):
                temptative_synonyms = MoleculeResolver.filter_and_sort_synonyms(
                    synonyms + temptative_synonyms
                )
                if len(temptative_synonyms) == 0:
                    possible_name = ""
                    while True:
                        possible_name = prompt_session.prompt(
                            "Input one name for the molecule", default=possible_name
                        )
                        if yes_no_dialog("Is this name correct?").run():
                            break
                    temptative_synonyms = [possible_name]

                temp_mol = MoleculeResolver.get_from_SMILES(temptative_SMILES)

                if temp_mol is not None:
                    MoleculeResolver.show_and_pause(temp_mol, temptative_synonyms[0])
                else:
                    print(
                        "Error: Image could not be gernerated from SMILES, this however does not mean always that the SMILES ist wrong."
                    )

                if yes_no_dialog(
                    title="Confirmation", text="Is the molecule structure correct?"
                ).run():
                    SMILES = temptative_SMILES
                    synonyms = temptative_synonyms
                    CAS.update(temptative_CAS)
                    temptative_CAS = CAS
                    start_geometry_comment = temptative_start_geometry_comment

                    while True:
                        temptative_CAS = prompt_session.prompt(
                            "Input valid CAS or empty:",
                            default=",".join(temptative_CAS),
                        )

                        if len(temptative_CAS) == 0 or all(
                            [
                                MoleculeResolver.is_valid_CAS(x.strip())
                                for x in temptative_CAS.split(",")
                            ]
                        ):
                            CAS = temptative_CAS.split(",")
                            break
        SMILES = MoleculeResolver.standardize_SMILES(SMILES, standardize)

        return Molecule(
            SMILES,
            synonyms,
            list(CAS),
            start_geometry_comment,
            mode_used,
            "interactive",
            1,
        )

    def find_single_molecule_cross_checked(
        self,
        identifiers: list[str],
        modes: list[str] = ["name"],
        required_formula: Optional[str] = None,
        required_charge: Optional[int] = None,
        required_structure_type: Optional[str] = None,
        services_to_use: Optional[list[str]] = None,
        standardize: bool = True,
        search_iupac_name: bool = False,
        minimum_number_of_cross_checks: Optional[int] = 1,
        try_to_choose_best_structure: bool = True,
        ignore_exceptions: bool = False,
    ) -> Union[Optional[Molecule], list[Optional[Molecule]]]:
        if services_to_use is None:
            services_to_use = MoleculeResolver._available_services

        if minimum_number_of_cross_checks is None:
            minimum_number_of_cross_checks = 1
        if not minimum_number_of_cross_checks <= len(services_to_use):
            raise ValueError(
                "The minimum_number_of_cross_checks exceeds the number of services that are used."
            )

        molecules = []

        for service in services_to_use:
            molecule = self.find_single_molecule(
                identifiers=identifiers,
                modes=modes,
                required_formula=required_formula,
                required_charge=required_charge,
                required_structure_type=required_structure_type,
                services_to_use=[service],
                standardize=standardize,
                search_iupac_name=search_iupac_name,
                ignore_exceptions=ignore_exceptions,
            )

            molecules.append(molecule)

        filtered_molecules = MoleculeResolver.filter_molecules(
            molecules,
            required_formula,
            required_charge,
            required_structure_type,
            standardize,
        )

        if not filtered_molecules:
            return None

        grouped_molecules = MoleculeResolver.group_molecules_by_structure(
            filtered_molecules, False
        )

        maximum_number_of_crosschecks_found = max(
            [len(v) for v in grouped_molecules.values()]
        )
        SMILES_with_highest_number_of_crosschecks = []
        for group_SMILES, group_molecules in grouped_molecules.items():
            if len(group_molecules) >= minimum_number_of_cross_checks:
                if len(group_molecules) == maximum_number_of_crosschecks_found:
                    SMILES_with_highest_number_of_crosschecks.append(group_SMILES)

        if try_to_choose_best_structure:
            SMILES_preferred = sorted(SMILES_with_highest_number_of_crosschecks)[0]
            if len(SMILES_with_highest_number_of_crosschecks) > 1:
                # if SMILES are the same ignoring isomeric info, use the more specific one:
                unique_non_isomeric_SMILES = set(
                    [
                        self.standardize_SMILES(smi, standardize, isomeric_SMILES=False)
                        for smi in SMILES_with_highest_number_of_crosschecks
                    ]
                )
                if len(unique_non_isomeric_SMILES) == 1:
                    SMILES_preferred = sorted(
                        SMILES_with_highest_number_of_crosschecks, key=len
                    )[-1]
                else:
                    # trust opsin algorithm: if not sure and opsin available
                    SMILES_preferred_by_opsin = None
                    for SMILES in SMILES_with_highest_number_of_crosschecks:
                        for molecule in grouped_molecules[SMILES]:
                            if molecule.mode == "name" and molecule.service == "opsin":
                                SMILES_preferred_by_opsin = SMILES

                    # if opsin result not available, or searched by another mode
                    # try getting all structures from the names and see if they agree
                    # with the SMILES found
                    if not SMILES_preferred_by_opsin:
                        SMILES_map = []
                        names_map = []
                        for SMILES in SMILES_with_highest_number_of_crosschecks:
                            for molecule in grouped_molecules[SMILES]:
                                for name in molecule.synonyms:
                                    SMILES_map.append(SMILES)
                                    names_map.append(name)

                        use_opsin_batch = False  # self._OPSIN_tempfolder is not None
                        if use_opsin_batch:
                            opsin_results = self.get_molecule_from_OPSIN_batchmode(
                                names_map
                            )
                        else:
                            opsin_results = [
                                self.get_molecule_from_OPSIN(name) for name in names_map
                            ]

                        SMILES_preferred_by_opsin = []
                        for (
                            original_SMILES_found,
                            molecule_found_by_opsin_from_synonym,
                        ) in zip(SMILES_map, opsin_results):
                            if molecule_found_by_opsin_from_synonym:
                                if (
                                    original_SMILES_found
                                    == molecule_found_by_opsin_from_synonym.SMILES
                                ):
                                    SMILES_preferred_by_opsin.append(
                                        original_SMILES_found
                                    )

                        SMILES_preferred_by_opsin = set(SMILES_preferred_by_opsin)
                        if len(SMILES_preferred_by_opsin) == 1:
                            SMILES_preferred_by_opsin = SMILES_preferred_by_opsin.pop()
                        else:
                            SMILES_preferred_by_opsin = None

                    if SMILES_preferred_by_opsin:
                        SMILES_preferred = SMILES_preferred_by_opsin
                    else:
                        # usually when chebi does not agree with others chebi is wrong
                        # this is used to our advantage here
                        SMILES_not_from_chebi = []
                        for (
                            temptative_SMILES,
                            temptative_molecules,
                        ) in grouped_molecules.items():
                            molecules_from_chebi = [
                                t_mol
                                for t_mol in temptative_molecules
                                if "chebi" in t_mol.service
                            ]
                            if len(molecules_from_chebi) == 0:
                                SMILES_not_from_chebi.extend(
                                    [temptative_SMILES] * len(temptative_molecules)
                                )

                        c = collections.Counter(SMILES_not_from_chebi).most_common()
                        if len(c) == 1 or (len(c) > 1 and c[0][1] > c[1][1]):
                            SMILES_preferred = c[0][0]
                        else:
                            temp = len(SMILES_with_highest_number_of_crosschecks)
                            warnings.warn(
                                f"\n\n{temp} molecules were found equally as often. First one sorted by SMILES was taken: \n{grouped_molecules}\n"
                            )
            molec = MoleculeResolver.combine_molecules(
                grouped_molecules[SMILES_preferred]
            )
            molec.found_molecules.append(grouped_molecules)
            return molec
        else:
            return [
                MoleculeResolver.combine_molecules(grouped_molecules[SMILES])
                for SMILES in SMILES_with_highest_number_of_crosschecks
            ]

    def find_multiple_molecules_parallelized(
        self,
        identifiers: list[str],
        modes: list[str],
        required_formulas: Optional[list[str]] = None,
        required_charges: Optional[list[int]] = None,
        required_structure_types: Optional[list[str]] = None,
        services_to_use: Optional[list[str]] = None,
        standardize: bool = True,
        search_iupac_name: bool = False,
        minimum_number_of_cross_checks: Optional[int] = 1,
        try_to_choose_best_structure: bool = True,
        progressbar: bool = True,
        max_workers: int = 5,
        ignore_exceptions: bool = True,
    ) -> list[Optional[Molecule]]:
        # reinitialize session
        self._session = None
        self._init_session(pool_maxsize=max_workers * 2)

        if services_to_use is None:
            services_to_use = MoleculeResolver._available_services

        services_to_use = [service.lower() for service in services_to_use]
        if required_formulas is None:
            required_formulas = [None] * len(identifiers)

        if required_structure_types is not None:
            if isinstance(required_structure_types, str):
                required_structure_types = [required_structure_types] * len(identifiers)
        else:
            required_structure_types = [None] * len(identifiers)
        if required_charges is not None:
            if isinstance(required_charges, int):
                required_charges = [required_charges] * len(identifiers)
        else:
            required_charges = [None] * len(identifiers)
        MoleculeResolver._check_parameters(
            modes=modes,
            services=services_to_use,
            required_charges=required_charges,
            required_formulas=required_formulas,
            required_structure_types=required_structure_types,
            context="find_multiple",
        )

        # when lots of molecules, use batch capabilities if possible
        message_getting_from_services = "Getting data from services"
        if len(identifiers) > 100:
            message_getting_from_services = "Getting data from all other services"

            services_asked_for_with_batch_capabilities = set(
                services_to_use
            ).intersection(set(self._available_services_with_batch_capabilities))
            services_asked_for_with_batch_capabilities = sorted(
                list(services_asked_for_with_batch_capabilities)
            )

            # sometimes there are issues with opsin, so always try to get from opsin first
            # in the case of an error you don't loose the rest
            if "opsin" in services_asked_for_with_batch_capabilities:
                services_asked_for_with_batch_capabilities.remove("opsin")
                services_asked_for_with_batch_capabilities.insert(0, "opsin")

            for service in services_asked_for_with_batch_capabilities:
                print("Getting data in batchmode for the following service:", service)

                if not MoleculeResolver._is_list_of_list_of_str(
                    identifiers
                ):  # Convert list[str] to list[list[str]] for usage in batchmode
                    identifiers = [[idf] for idf in identifiers]
                if not MoleculeResolver._is_list_of_list_of_str(modes):
                    modes = [[md] for md in modes]

                self.get_molecules_using_batchmode_from(
                    identifiers,
                    modes,
                    service,
                    progressbar=progressbar,
                    ignore_exceptions=ignore_exceptions,
                )

        def _find(generator):
            return list(
                tqdm(generator, total=len(identifiers), disable=not progressbar)
            )

        args = []
        for (
            identifier,
            mode,
            required_formula,
            required_charge,
            required_structure_type,
        ) in zip(
            identifiers,
            modes,
            required_formulas,
            required_charges,
            required_structure_types,
        ):
            if minimum_number_of_cross_checks is None:
                args.append(
                    (
                        identifier,
                        mode,
                        required_formula,
                        required_charge,
                        required_structure_type,
                        services_to_use,
                        standardize,
                        search_iupac_name,
                        False,
                        ignore_exceptions,
                    )
                )
            else:
                args.append(
                    (
                        identifier,
                        mode,
                        required_formula,
                        required_charge,
                        required_structure_type,
                        services_to_use,
                        standardize,
                        search_iupac_name,
                        minimum_number_of_cross_checks,
                        try_to_choose_best_structure,
                        ignore_exceptions,
                    )
                )

        print(message_getting_from_services)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            if not minimum_number_of_cross_checks:
                results = _find(executor.map(self.find_single_molecule, *zip(*args)))
            else:
                results = _find(
                    executor.map(self.find_single_molecule_cross_checked, *zip(*args))
                )

        return results
