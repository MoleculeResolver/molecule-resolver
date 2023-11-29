import pytest
from moleculeresolver import MoleculeResolver
import json
import os
from pathlib import Path
from typing import Any, Callable, Dict, Optional


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
        with MoleculeResolver() as mr:
            iupac_name = data["iupac_name"]
            self._test_service(
                mr.get_molecule_from_OPSIN,
                iupac_name,
                SMILES,
                data["SMILES"],
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
