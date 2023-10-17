from moleculeresolver import MoleculeResolver

with MoleculeResolver(available_service_API_keys={"chemeo": 'YOUR_CHEMEO_API_KEY'}, molecule_cache_db_path='test.db') as mr:

    names_to_find = ['2-bromobutane', 'ethanol', 'methanol', 'propane', 'butane']
    
    # search for the names in parallel
    all_names = []
    all_modes = []
    for name in names_to_find:
        names = mr.expand_name_heuristically(name)
        all_names.append(names)
        all_modes.append(['name'])

    molecules_found_in_parallel = mr.find_multiple_molecules_parallelized(all_names, all_modes)
    print('all_found_in_parallel:', all(molecules_found_in_parallel))

    # search for the names sequentially
    molecules = []
    for name in names_to_find:
        names = mr.expand_name_heuristically(name)
        molecule = mr.find_single_molecule_cross_checked(names, 'name')
        molecules.append(molecule)
    print('all_found:', all(molecules))

    # search for CAS numbers
    molecules_found_by_CAS = []
    CAS_numbers = ['7732-18-5', '78-76-2', '64-17-5', '67-56-1', '74-98-6', '106-97-8']
    for CAS in CAS_numbers:
        molecule = mr.find_single_molecule_cross_checked(CAS, 'CAS')
        molecules_found_by_CAS.append(molecule)

    print('all_found_by_CAS:', all(molecules_found_by_CAS))