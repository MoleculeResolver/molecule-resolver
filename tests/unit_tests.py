import pytest
from moleculeresolver import MoleculeResolver
import json

with open("benchmark_component_molecules.json", "r") as f:
    benchmark = json.load(f)
    breakpoint()

@pytest.mark.parametrize("name, data", benchmark.items())
def test_opsin(name, data):
    with MoleculeResolver() as mr:
        res = mr.get_molecule_from_OPSIN(name)
        if res is None:
            assert data["SMILES"] == "None"
        assert res.SMILES == data["SMILES"]

if __name__ == "__main__":
    test_opsin()
