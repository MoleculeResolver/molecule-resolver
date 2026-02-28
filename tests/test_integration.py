from datetime import datetime
import pytest
from moleculeresolver import Molecule, MoleculeResolver
import json
import os
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Union
from tqdm import tqdm
from dataclasses import dataclass


# IUPAC names
dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
with open(dir_path / "benchmark_component_molecules_iupac.json", "r") as f:
    benchmark = json.load(f)

RESPONSES_PATH = dir_path / "responses.json"
SMILES = "SMILES"


class MoleculeResolverPatched(MoleculeResolver):
    def __init__(
        self,
        json_path: Optional[str] = None,
        available_service_API_keys: Optional[dict[str, Optional[str]]] = None,
        molecule_cache_db_path: Optional[str] = None,
        molecule_cache_expiration: Optional[datetime] = None,
        differentiate_isotopes: bool = False,
    ):
        super().__init__(
            available_service_API_keys,
            molecule_cache_db_path,
            molecule_cache_expiration,
            differentiate_isotopes,
        )
        self.json_path = json_path
        if self.json_path is not None:
            self.json_path = Path(self.json_path)

        if self.json_path is not None and self.json_path.exists():
            with open(self.json_path, "r") as f:
                self.json_data = json.load(f)
        else:
            self.json_data = {}
        
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
    ) -> Optional[str]:
        if self.json_path is not None and url in self.json_data:
            res =  self.json_data[url]
        else:
            res =  super()._resilient_request(
                url,
                kwargs,
                request_type,
                accepted_status_codes,
                rejected_status_codes,
                max_retries,
                sleep_time,
                allow_redirects,
                json,
                return_response,
            )
            if self.json_path is not None and url not in self.json_data:
                self.json_data[url] = res
        return res


    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.json_path is not None:
            with open(self.json_path, "w") as f:
                json.dump(self.json_data, f)
        return super().__exit__(exc_type, exc_value, exc_traceback)
    



@pytest.mark.parametrize("data", benchmark.values())
class TestServices:
    json_path: Optional[str] = RESPONSES_PATH

    @staticmethod
    def _test_service(
        call_method: Callable,
        input_identifier: str,
        output_identifier_type: str,
        output_identifier,
        kwargs: Optional[Dict] = None,
    ):
        """
        Test a service by calling it with an input identifier and checking that the output identifier matches the expected value.

        Parameters
        ----------
        call_method : Callable
            The method to call
        input_identifier : str
            The input identifier
        output_identifier_type : str
            The type of the output identifier
        output_identifier : str
            The expected output identifier
        kwargs : Optional[Dict], optional
            Additional keyword arguments to pass to the call method, by default None


        """
        if kwargs is None:
            kwargs = {}
        res = call_method(input_identifier, **kwargs)
        if res is None:
            raise ValueError(f"No molecule found for {input_identifier}")

        res_txt = res.__dict__[output_identifier_type]
        if res_txt == output_identifier:
            return
        else:
            raise ValueError(f"Expected {output_identifier} but got {res_txt}")

    def test_opsin(self, data):
        with MoleculeResolverPatched(json_path=self.json_path) as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_OPSIN,
                iupac_name,
                SMILES,
                data["SMILES"],
            )

    # def test_pubchem(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_pubchem,
    #             iupac_name,
    #             SMILES,
    #             data["SMILES"],
    #             {"mode": "name"},
    #         )

    # def test_comptox(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_CompTox,
    #             iupac_name,
    #             SMILES,
    #             data["SMILES"],
    #             {"mode": "name"},
    #         )

    # def test_cts(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_CTS,
    #             iupac_name,
    #             SMILES,
    #             data["SMILES"],
    #             {"mode": "name"},
    #         )

    # # Need API key
    # # def test_chemeo(self, data):
    # #     with MoleculeResolver() as mr:
    # #         iupac_name = data["iupac_name"]
    # #         self._test_service(
    # #             mr.get_molecule_from_Chemeo,
    # #             iupac_name,
    # #             SMILES,
    # #             data["SMILES"],
    # #             {"mode": "name"},
    # #         )

    # def test_cas(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_CAS_registryx,
    #             iupac_name,
    #             SMILES,
    #             data["SMILES"],
    #             {"mode": "name"},
    #         )

    # def test_cir(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_CIR,
    #             iupac_name,
    #             SMILES,
    #             data["SMILES"],
    #             {"mode": "name"},
    #         )

    # def test_nist(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_NIST,
    #             iupac_name,
    #             SMILES,
    #             data["SMILES"],
    #             {"mode": "name"},
    #         )

    # def test_chebi(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_ChEBI, iupac_name, "chebi_id", data["chebi_id"]
    #         )


# def test_opsin_batchmode():
#     names = [d["iupac_name"] for d in benchmark.values()]
#     smiles = [d["SMILES"] for d in benchmark.values()]
#     with MoleculeResolver() as mr:
#         res = mr.get_molecule_from_OPSIN_batchmode(names)
#     for i, r in enumerate(res):
#         if r[0].SMILES == smiles[i]:
#             continue
#         else:
#             raise ValueError("Expected " + smiles[i] + " but got " + r.SMILES)


def generate_data(json_path, overwrite_json: bool = False):
    """ Generate data for unit tests
    
    Parameters
    ----------
    json_path : str
        Path to JSON file to save data to
    overwrite_json : bool, optional
        Whether to overwrite the JSON file, by default False
    
    """

    # Remove JSON path
    if overwrite_json:
        if os.path.exists(json_path):
            os.remove(json_path)

    # Get all test methods
    test_services = TestServices(json_path=json_path)
    methods  = dir(TestServices)
    test_methods = [m for m in methods if m.startswith("test_")]

    # Run tests
    bar = tqdm(test_methods, desc="Running tests")
    for m in bar:
        method = getattr(test_services, m)
        bar.set_description(f"Running {m}")
        for data in benchmark.values():
            method(data)

def test_multi_name_exhaustive_benchmark_case(monkeypatch):
    mr = MoleculeResolver()
    mr.supported_modes_by_services = {"svc1": ["name"], "svc2": ["name"]}
    mr._available_services = ["svc1", "svc2"]

    monkeypatch.setattr(mr, "_check_parameters", lambda **_: None)
    monkeypatch.setattr(
        mr,
        "_check_and_flatten_identifiers_and_modes",
        lambda identifiers, modes: (
            ["isopropanol", "2-propanol"],
            ["name", "name"],
            ["isopropanol", "2-propanol"],
            set(),
            None,
        ),
    )
    monkeypatch.setattr(mr, "standardize_SMILES", lambda smiles: smiles)

    result_map = {
        ("svc1", "isopropanol"): {
            "SMILES": "CC(O)C",
            "synonyms": ["isopropanol"],
            "CAS": {"67-63-0"},
            "additional_information": "svc1",
            "mode_used": "name",
            "identifier_used": "isopropanol",
            "service": "svc1",
            "cas_is_authoritative": False,
        },
        ("svc1", "2-propanol"): {
            "SMILES": "CC(C)O",
            "synonyms": ["2-propanol"],
            "CAS": {"67-63-0"},
            "additional_information": "svc1",
            "mode_used": "name",
            "identifier_used": "2-propanol",
            "service": "svc1",
            "cas_is_authoritative": False,
        },
        ("svc2", "isopropanol"): {
            "SMILES": "CC(C)O",
            "synonyms": ["isopropanol"],
            "CAS": {"67-63-0"},
            "additional_information": "svc2",
            "mode_used": "name",
            "identifier_used": "isopropanol",
            "service": "svc2",
            "cas_is_authoritative": False,
        },
    }

    monkeypatch.setattr(
        mr,
        "_resolve_single_service_candidate",
        lambda service, identifier, mode, *_: result_map.get((service, identifier)),
    )

    result = mr.find_single_molecule(
        identifiers=["isopropanol", "2-propanol"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        search_strategy="exhaustive",
    )

    assert result is not None
    assert result.SMILES == "CC(C)O"


def test_consensus_vs_legacy_mode_comparison(monkeypatch):
    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    legacy_winner = Molecule(
        "CCO", ["ethanol"], [], "svc1", "name", "svc1", 1, "ethanol"
    )
    consensus_winner = Molecule(
        "CCO", ["ethyl alcohol"], [], "svc2", "name", "svc2", 1, "ethyl alcohol"
    )

    def fake_find_single_molecule(*args, **kwargs):
        service = kwargs["services_to_use"][0]
        if service == "svc1":
            return legacy_winner
        return consensus_winner

    monkeypatch.setattr(mr, "find_single_molecule", fake_find_single_molecule)
    monkeypatch.setattr(
        mr,
        "filter_molecules",
        lambda molecules, *_: [molecule for molecule in molecules if molecule is not None],
    )
    monkeypatch.setattr(
        mr,
        "group_molecules_by_structure",
        lambda molecules, *_: {"CCO": molecules},
    )

    legacy_result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol"],
        modes=["name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="legacy",
    )
    consensus_result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol"],
        modes=["name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="consensus",
    )

    assert legacy_result is not None
    assert consensus_result is not None
    assert legacy_result.SMILES == consensus_result.SMILES


def test_strict_isomer_acceptance_case(monkeypatch):
    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    mol_a = Molecule("CCO", ["ethanol"], [], "svc1", "name", "svc1", 1, "ethanol")
    mol_b = Molecule(
        "C[C@H](O)C", ["(S)-2-butanol"], [], "svc2", "name", "svc2", 1, "(S)-2-butanol"
    )

    def fake_find_single_molecule(*args, **kwargs):
        service = kwargs["services_to_use"][0]
        return mol_a if service == "svc1" else mol_b

    monkeypatch.setattr(mr, "find_single_molecule", fake_find_single_molecule)
    monkeypatch.setattr(
        mr,
        "filter_molecules",
        lambda molecules, *_: [molecule for molecule in molecules if molecule is not None],
    )
    monkeypatch.setattr(
        mr,
        "group_molecules_by_structure",
        lambda molecules, *_: {"CCO": [mol_a], "C[C@H](O)C": [mol_b]},
    )
    monkeypatch.setattr(
        mr,
        "_collect_opsin_isomer_matches",
        lambda grouped, smiles: {smiles[0]: False, smiles[1]: True},
    )

    result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol", "(S)-2-butanol"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="strict_isomer",
    )

    assert result is not None
    assert result.SMILES == "C[C@H](O)C"


if __name__ == "__main__":
    generate_data(RESPONSES_PATH, overwrite_json=True)
