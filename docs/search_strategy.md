# Search Strategy

`find_single_molecule`, `find_single_molecule_crosschecked`, and
`find_multiple_molecules_parallelized` now support `search_strategy`.

## Options

- `first_hit` (default): preserves the previous behavior and returns as soon as
  the first matching candidate is found in service order.
- `exhaustive`: evaluates all provided identifier/mode combinations across the
  selected services, deduplicates candidates, and chooses the best supported
  structure group.

## Example

```python
from moleculeresolver import MoleculeResolver

with MoleculeResolver() as mr:
    molecule = mr.find_single_molecule(
        identifiers=["ethanol", "ethyl alcohol"],
        modes=["name", "name"],
        search_strategy="exhaustive",
    )
```
