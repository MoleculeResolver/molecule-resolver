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


def test_opsin_batchmode():
    names = [d["iupac_name"] for d in benchmark.values()]
    smiles = [d["SMILES"] for d in benchmark.values()]
    with MoleculeResolver() as mr:
        res = mr.get_molecule_from_OPSIN_batchmode(names)
    for i, r in enumerate(res):
        if r[0].SMILES == smiles[i]:
            continue
        else:
            raise ValueError("Expected " + smiles[i] + " but got " + r.SMILES)
            
@pytest.mark.parametrize("data", benchmark.values())
def test_pubchem(data):
    with MoleculeResolver() as mr:
        iupac_name = data["iupac_name"]
        res = mr.get_molecule_from_pubchem(iupac_name, mode="name")
        if res is None:
            raise ValueError("No molecule found for " + iupac_name)
        if res.SMILES == data["SMILES"]:
            return
        else:
            raise ValueError("Expected " + data["SMILES"] + " but got " + res.SMILES)

if __name__ == "__main__":
    test_opsin_batchmode()