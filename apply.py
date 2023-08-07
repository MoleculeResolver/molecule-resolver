from moleculeresolver import MoleculeResolver

with MoleculeResolver(available_service_API_keys={"chemeo": 'YOUR_CHEMEO_API_KEY'}) as cf:

    molecules = []
    for name in ['2-bromobutane', 'ethanol', 'methanol', 'propane', 'butane']:
        names = cf.expand_name_heuristically(name)
        molecule = cf.find_single_molecule_cross_checked(names, 'name') # optional
        molecules.append(molecule)