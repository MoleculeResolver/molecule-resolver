# Runtime Context Lifecycle

`MoleculeResolver` context management is intentionally limited to runtime resources.

## Enter

The `with MoleculeResolver() as mr:` enter step performs only:

1. RDKit log suppression context setup.
2. Molecule cache context setup.
3. OPSIN temp folder creation (when OPSIN batch support is enabled).

## Exit

The exit step performs only:

1. Molecule cache context teardown.
2. RDKit log suppression context teardown.
3. OPSIN temp folder cleanup (skipped when an exception is bubbling up).

This separation keeps search and resolution logic outside lifecycle handling, so runtime setup/teardown stays predictable and easier to maintain.
