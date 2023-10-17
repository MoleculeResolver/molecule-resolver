# MoleculeResolver

A python package allowing to use several web services to find molecule structures/names/CAS based on their identifiers such as name, CAS, SMILES, InChI, etc.


## Installation

The package is available on [pypi](https://pypi.org/project/molecule-resolver/).

```
pip intall molecule-resolver
```

## Usage

Other examples available in `apply.py`, it is supposed to be called as context manager:

```python
with MoleculeResolver(available_service_API_keys={'chemeo': 'YOUR_API_KEY'}) as cf:
    molecule = cf.find_single_molecule_cross_checked(['ethanol'], ['name'], minimum_number_of_cross_checks=1)
```

If you call it with as a context manager it will automatically silence all output from rdkit, as this can be substantial when parsing the data from the different web services.

If you want to have more control over what is muted from the rdkit output, you can also use the `MoleculeResolver` class like so:

```python
from moleculeresolver.rdkitmods import disabling_rdkit_logger

with disabling_rdkit_logger(mute_errors = True, mute_warning = True, mute_info = True, mute_debug = True):
    cf = MoleculeResolver(available_service_API_keys={'chemeo': 'YOUR_API_KEY'})
    molecule = cf.find_single_molecule_cross_checked(['ethanol'], ['name'], minimum_number_of_cross_checks=1)
```

Under `moleculeresolver.rdkitmods` you will find the `disabling_rdkit_logger` class which can be used as a context manager or as a decorator that you can use if you want to mute input on one function alone.

```python
from moleculeresolver.rdkitmods import disabling_rdkit_logger

@disabling_rdkit_logger(mute_errors = True, mute_warning = True, mute_info = True, mute_debug = True)
def yourfunction():
    pass
```