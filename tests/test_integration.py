import pytest
from moleculeresolver import MoleculeResolver
import json
import os 
from pathlib import Path


# IUPAC names
dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
with open(dir_path / "benchmark_component_molecules_iupac.json", "r") as f:
    benchmark = json.load(f)


@pytest.mark.parametrize("data", benchmark.values())
def test_opsin(data):
    with MoleculeResolver() as mr:
        iupac_name = data["iupac_name"]
        res = mr.get_molecule_from_OPSIN(iupac_name)
        if res is None:
            raise ValueError("No molecule found for " + iupac_name)
        if res.SMILES == data["SMILES"]:
            return
        else:
            raise ValueError("Expected " + data["SMILES"] + " but got " + res.SMILES)


if __name__ == "__main__":
    for d in benchmark.values():
        test_opsin(d)