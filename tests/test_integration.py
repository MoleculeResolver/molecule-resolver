from datetime import datetime
import pytest
from moleculeresolver import MoleculeResolver
import json
import os
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Union


# IUPAC names
dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
with open(dir_path / "benchmark_component_molecules_iupac.json", "r") as f:
    benchmark = json.load(f)

SMILES = "SMILES"

# PATCH_STATE = "SAVE"


# class PatchResilientRequest:
#     def __init__(self, json_data, patch_state):
#         self.json_data = json_data
#         self.patch_state = patch_state

#     def __call__(self, url: str, **kwargs) -> str:
#         if self.patch_state == "SAVE":
#             self.json_data[url] = kwargs["json"]
#         elif self.patch_state == "LOAD":
#             return self.json_data[url]

# Add ability to call TestServices with a patch of resilient request that saves response
# Use Mock that wraps _resilient_request and has a side effect that saves the json response.


class MoleculeResolverSaver(MoleculeResolver):
    def __init__(
        self,
        save_path: Optional[str] = None,
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
        self.save_path = save_path
        if self.save_path is not None:
            self.save_path = Path(self.save_path)

        if self.save_path is not None and self.save_path.exists():
            with open(self.save_path, "r") as f:
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
        if self.save_path is not None:
            self.json_data[url] = res
        return res


    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.save_path is not None:
            with open(self.save_path, "w") as f:
                json.dump(self.json_data, f)
        return super().__exit__(exc_type, exc_value, exc_traceback)
    




@pytest.mark.parametrize("data", benchmark.values())
class TestServices:
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
        with MoleculeResolverSaver(save_path="test.json") as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_OPSIN,
                iupac_name,
                SMILES,
                data["SMILES"],
            )

    def test_pubchem(self, data):
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_pubchem,
                iupac_name,
                SMILES,
                data["SMILES"],
                {"mode": "name"},
            )

    def test_comptox(self, data):
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_CompTox,
                iupac_name,
                SMILES,
                data["SMILES"],
                {"mode": "name"},
            )

    def test_cts(self, data):
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_CTS,
                iupac_name,
                SMILES,
                data["SMILES"],
                {"mode": "name"},
            )

    # Need API key
    # def test_chemeo(self, data):
    #     with MoleculeResolver() as mr:
    #         iupac_name = data["iupac_name"]
    #         self._test_service(
    #             mr.get_molecule_from_Chemeo,
    #             iupac_name,
    #             SMILES,
    #             data["SMILES"],
    #             {"mode": "name"},
    #         )

    def test_cas(self, data):
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_CAS_registryx,
                iupac_name,
                SMILES,
                data["SMILES"],
                {"mode": "name"},
            )

    def test_cir(self, data):
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_CIR,
                iupac_name,
                SMILES,
                data["SMILES"],
                {"mode": "name"},
            )

    def test_nist(self, data):
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_NIST,
                iupac_name,
                SMILES,
                data["SMILES"],
                {"mode": "name"},
            )

    def test_chebi(self, data):
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_ChEBI, iupac_name, "chebi_id", data["chebi_id"]
            )


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


def generate_data():
    # Run each test with a patch of resilient request that saves response
    pass


if __name__ == "__main__":
    generate_data()
