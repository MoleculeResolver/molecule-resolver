
# MoleculeResolver

The **moleculeresolver** was born out of the need to annotate large datasets with accurate structural information fast and to crosscheck whether given metadata (name, SMILES) agrees with each other. It also allows to efficiently compare whether structures are available in two large datasets. 

In short it's a Python module that allows you to retrieve molecular structures from multiple chemical databases, perform crosschecks to ensure data reliability, and standardize the best representation of molecules. It also provides functions for comparing molecules and sets of molecules based on specific configurations. This makes it a useful tool for researchers, chemists, or anyone working in computational chemistry / cheminformatics who needs to ensure they are working with the best available data for a molecule.

## Installation

The package is available on [pypi](https://pypi.org/project/molecule-resolver/):

```sh
pip install molecule-resolver
```
While the source code is available here: [https://github.com/MoleculeResolver/molecule-resolver](https://github.com/MoleculeResolver/molecule-resolver)

## Features

- **🔍 Retrieve Molecular Structures**: Fetch molecular structures from different chemical databases, including PubChem, Comptox, Chemo, and others.
- **🆔 Support for Different Identifier Types**: Retrieve molecular structures using a variety of identifier types, including CAS numbers, SMILES, InChI, InChIkey and common names.
- **✅ Cross-check Capabilities**: Use data from multiple sources to verify molecular structures and identify the best representation.
- **🔄 Molecule Comparison**: Compare molecules or sets of molecules based on their structure, properties, and specified ⚙️ configurations.
- **⚙️ Standardization**: Standardize molecular structures, including handling isomers, tautomers, and isotopes.
- **💾 Caching Mechanism**: Use local caching to store molecules and reduce the number of repeated requests to external services, improving performance and reducing latency.

## Services used
At this moment, the following services are used to get the best structure for a given identifier. In the future, this list might be reviewed to improve perfomance, adding new services or removing some.
In case you want to add an additional service, open an issue or a pull request.

The MoleculeResolver does not offer all options/configurations for each service available with the specific related repos as it focusses on getting the structure based on the identifiers and doing so as accurate as possible while still being fast using parallelization under the hood.
| Service                                                                 | Name | CAS | Formula | SMILES | InChI | InChIKey | CID | Batch search | Repos                                                                 |
|-------------------------------------------------------------------------|------|-----|---------|--------|-------|----------|-----|--------------------|------------------------------------------------------------------------------|
| [cas_registry](https://commonchemistry.cas.org/)                        | ✅   | ✅  | ❌      | ✅     | ✅    | ❌       | ❌  | ❌                 |  |
| [chebi](https://www.ebi.ac.uk/chebi/)                                   | ✅   | ✅  | ✅      | ✅     | ✅    | ✅       | ❌  | ❌                 |                               |
| [chemeo](https://www.chemeo.com/)                                       | ✅   | ✅  | ❌      | ✅     | ✅    | ✅       | ❌  | ❌                 |                                                        |
| [cir](https://cactus.nci.nih.gov/chemical/structure)                    | ✅   | ✅  | ✅      | ✅     | ✅    | ✅       | ❌  | ❌                 | - [CIRpy](https://github.com/mcs07/CIRpy "wrapper for the CIR. FYI, CIR uses OPSIN under the hood, unless specified otherwise.")                                      |
| [comptox](https://comptox.epa.gov/dashboard)                            | ✅   | ✅  | ❌      | ❌     | ❌    | ✅       | ❌  | ✅                 |                                                         |
| [cts](https://cts.fiehnlab.ucdavis.edu/)                                | (✅)   | ✅  | ❌      | ✅     | ❌    | ❌       | ❌  | ❌                 |                                                        |
| [nist](https://webbook.nist.gov/chemistry/)                             | ✅   | ✅  | ✅      | ✅     | ❌    | ❌       | ❌  | ❌                 | - [NistChemPy](https://github.com/IvanChernyshov/NistChemPy "unofficial wrapper for search and data extraction of the NIST Webbook.")                                                       |
| [opsin](https://opsin.ch.cam.ac.uk/)                                    | ✅   | ❌  | ❌      | ❌     | ❌    | ❌       | ❌  | ✅                 | - [py2opsin](https://github.com/JacksonBurns/py2opsin "lightweight OPSIN wrapper only depending on having Java installed.") <br> - [pyopsin](https://github.com/Dingyun-Huang/pyopsin "lightweight OPSIN wrapper depending on having Java installed + additional dependencies.")           |
| [pubchem](https://pubchem.ncbi.nlm.nih.gov/)</li></ul>                           | ✅   | ✅  | ✅      | ✅     | ✅    | ✅       | ✅  | ✅                 | - [PubChemPy](https://github.com/mcs07/PubChemPy "wrapper for the pubchem PUG API")                              |
| [srs](https://cdxapps.epa.gov/oms-substance-registry-services/search)   | ✅   | ✅  | ❌      | ❌     | ❌    | ❌       | ❌  | ✅                 |                                                         |

ChemSpider was not used as it is already included in CIR [[1]](https://matt-swain.com/blog/2012-03-20-cirpy-python-nci-chemical-identifier-resolver) [[2]](https://cactus.nci.nih.gov/blog/?p=1456) [[3]](https://github.com/mcs07/ChemSpiPy). ChemIDplus and the Drug Information Portal were retired in 2022 [[4]](https://www.nlm.nih.gov/pubs/techbull/ja22/ja22_pubchem.html).

## 🚀 Usage

### Initialization

To use **Molecule Resolver**, first import and initialize the `MoleculeResolver` class. it is supposed to be used as a context manager:

```python
from moleculeresolver import MoleculeResolver

with MoleculeResolver(available_service_API_keys={"chemeo": "YOUR_API_KEY"}) as mr:
    ...
```

### Retrieve and Compare Molecules by Name and CAS

Retrieve a molecule using both its common name and CAS number, then compare the two to ensure they represent the same structure:

```python
from rdkit import Chem
from moleculeresolver import MoleculeResolver

with MoleculeResolver(available_service_API_keys={"chemeo": "YOUR_API_KEY"}) as mr:
    molecule_name = mr.find_single_molecule(["aspirin"], ["name"])
    molecule_cas = mr.find_single_molecule(["50-78-2"], ["cas"])
    
    are_same = mr.are_equal(Chem.MolFromSmiles(molecule_name.SMILES), 
                            Chem.MolFromSmiles(molecule_cas.SMILES))
    print(f"Are the molecules the same? {are_same}")
```

### Parallelized Molecule Retrieval and Saving to JSON

Use the parallelized version to retrieve multiple molecules. If a large number of molecules is searched, moleculeresolver will try to use batch download capabilities whenever the database supports this.

```python
import json
from moleculeresolver import MoleculeResolver

molecule_names = ["aspirin", "propanol", "ibuprofen", "non-exixtent-name"]
not_found_molecules = []
molecules_dicts = {}

with MoleculeResolver(available_service_API_keys={"chemeo": "YOUR_API_KEY"}) as mr:
    molecules = mr.find_multiple_molecules_parallelized(molecule_names, [["name"]] * len(molecule_names))
    for name, molecule in zip(molecule_names, molecules):
        if molecule:
            molecules_dicts[name] = molecule.to_dict(found_molecules='remove')
        else:
            not_found_molecules.append(name)

with open("molecules.json", "w") as json_file:
    json.dump(molecules_dicts, json_file, indent=4)

print(f"Molecules not found: {not_found_molecules}")
```

## ⚙️ Configuration

The `MoleculeResolver` class allows users to configure various options like:

- **API Keys**: Set API keys for accessing different molecular databases. Currently only chemeo needs one.
- **Standardization Options**: Choose how to handle molecular standardization (e.g., normalizing functional groups, disconnecting metals, handling isomers, etc.).
- **Differentiation Settings**: Options for distinguishing between isomers, tautomers, and isotopes.

## ⚠️ Warning

**Inchi** is included in the set of valid identifiers for various [services](#services-used). You should be aware that using Inchi to get SMILES using RDKit is not the most robust approach. You can read more about it [here](https://github.com/rdkit/rdkit/issues/542).

## 🤝 Contributing

Contributions are welcome! If you have suggestions for improving the Molecule Resolver or want to add new features, feel free to submit an issue or a pull request on GitHub.

