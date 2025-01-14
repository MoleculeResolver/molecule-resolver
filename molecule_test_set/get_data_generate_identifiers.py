
import json
import os
from moleculeresolver import MoleculeResolver
from rdkit import Chem

with open('benchmark_molecule_names.json', 'r') as f:
    benchmark_molecule_names = json.load(f)

if not os.path.exists('benchmark_component_molecules.json'):
        
    benchmark_component_molecules = {}
    with MoleculeResolver(available_service_API_keys={"chemeo": 'YOUR_KEY'}, molecule_cache_db_path='molecule_cache.db') as mr:
        mr._available_services.remove('cts')
        temp = mr.find_multiple_molecules_parallelized(benchmark_molecule_names, [['name']] * len(benchmark_molecule_names))

        for name, molecule in zip(benchmark_molecule_names, temp, strict=True):
            molecule.found_molecules = []
            molecule = molecule.__dict__
            mol = Chem.MolFromSmiles(molecule['SMILES'])
            pubchem_cid = [v.strip() for v in molecule['additional_information'].split(';') if 'pubchem' in v]
            if pubchem_cid:
                pubchem_cid = int(pubchem_cid[0].split(':')[-1])
            else:
                pubchem_cid = None

            molecule['pubchem_cid'] = pubchem_cid
            molecule['formula'] = Chem.rdMolDescriptors.CalcMolFormula(mol)
            molecule['hill_formula'] = mr.to_hill_formula(mol)
            molecule['inchi'] = Chem.MolToInchi(mol)
            molecule['inchikey'] = Chem.InchiToInchiKey(molecule['inchi'])
            benchmark_component_molecules[name] = molecule
    
    with open('benchmark_component_molecules.json', 'w') as f:
        json.dump(benchmark_component_molecules, f, indent=4)
