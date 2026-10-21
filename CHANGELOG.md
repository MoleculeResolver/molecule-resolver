# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- The `cts` service now uses [CTS-Lite](https://cts-lite.metabolomics.us/), because the
  old CTS REST API is gone (HTTP 410). Supported modes are now `smiles`, `inchi`,
  `inchikey`, `formula` and `cid`; `cas` and `name` are no longer supported.
  CTS-Lite guesses the identifier type itself, so results whose guessed type differs
  from the requested mode are discarded. See `docs/cts.md` for the limitations.

### Added

- Batch mode for `cts` (`get_molecules_from_CTS_batchmode`), 1000 identifiers per
  request by default.
- Integration tests for all services (#4, #34). Benchmark lookups are replayed from
  recorded responses in `tests/responses.json`; set `MOLECULERESOLVER_RECORD_RESPONSES=1`
  to record missing ones. Services known to fail for a benchmark entry are listed in
  its `expected_failures`. Tests marked `online` look up ethanol in every live service;
  for Chemeo and CAS Registry they need `CHEMEO_API_KEY` / `CAS_REGISTRY_API_KEY`,
  which CI takes from the repository secrets of the same name.

### Fixed

- Bare `raise` statements outside `except` blocks (molblock/url check, PubChem
  download) raised `RuntimeError: No active exception to reraise` instead of a
  `ValueError` with a message.
- The CAS Registry API key was validated with the Chemeo token regex.
- `cas_registry` rejected every API key (HTTP 401): requests now go to the
  `direct-api` endpoint with the `x-origin` header that CAS Common Chemistry requires.
  A rejected key now raises a `ValueError` right away instead of retrying.
- Error messages of failed requests printed the request headers including API keys;
  `X-API-KEY` and `Authorization` headers are now shown as `***`.
- `cas_registry` found no structure for molecules whose detail record has an empty
  `smile` field (e.g. urea); it now falls back to `canonicalSmile`.

## [0.6.0] - 2026-10-21

### Added

- This changelog.

- New structure type `"ionic_liquid"`: a salt whose cations are all organic (contain
  carbon and no metal). Anions are not checked. It can be used as
  `required_structure_type` in all search methods and in `check_structure_type`.

### Changed

- `get_structure_type_from_SMILES` now returns `"ionic_liquid"` instead of `"salt"` for
  salts with only organic cations. Code comparing its result to `"salt"` directly
  should use `check_structure_type(SMILES, "salt")` instead, which matches both
  `"salt"` and `"ionic_liquid"`.
- `"ionic_liquid"` is a stricter type of `"salt"`: searches with
  `required_structure_type="salt"` (including `find_salt_molecule`) still accept
  ionic liquids, and resonance structures are compared for ionic liquids as for salts.
- Bumped dependency `urllib3` from 2.7.0 to 2.8.0.
- Bumped dev dependency `anyio` from 4.9.0 to 4.14.2.

### Fixed

- `get_structure_type_from_SMILES` returns `None` instead of the string `"None"` if a
  part of the SMILES cannot be parsed.

## [0.5.4] - 2026-08-21

### Fixed

- Support for Python 3.14.

## [0.5.3] - 2026-07-26

### Added

- Service adapter layer: a `ServiceAdapter` base class and registry, with adapters for
  OPSIN, PubChem, CIR, ChEBI, CompTox, SRS, Chemeo, CTS, NIST and CAS Registry.
- `search_strategy` API with `"first_hit"` and `"exhaustive"` strategies, query
  expansion and deduplication of results.
- Weighted consensus scoring of candidates with a deterministic tie-break pipeline.
- `resolution_mode` API, including a `"strict_isomer"` mode that rejects stereo
  mismatches using an OPSIN isomer verification step.
- `include_evidence` option returning `CandidateEvidence` and `ResolutionResult`
  objects with service agreement, identifier concordance, synonym overlap and score
  breakdown.
- `differentiate_deuterium` option.

### Changed

- Scoring uses RDKit stereo/chirality metrics instead of SMILES length.
- Resolution and strategy parameters are typed with `Literal`s.
- Lifecycle (context manager) handling of cache, logger and temporary directory was
  extracted into helpers.
- Updated the bundled OPSIN version.
- Bumped dependencies `pillow` to 12.3.0, `idna` to 3.15, `urllib3` to 2.7.0 and
  `requests` to 2.33.0, and dev dependencies `starlette` to 1.3.1, `pygments` to
  2.20.0 and `pytest` to 9.0.3.

## [0.5.2] - 2026-02-24

### Changed

- Improved CompTox results.
- RDKit log output is disabled.
- Dependency version constraints were relaxed to allow other versions.

### Fixed

- Regular expressions use raw strings.

## [0.5.1] - 2026-02-22

### Fixed

- ChEBI and CompTox lookups.

## [0.5.0] - 2026-02-20

### Changed

- Much improved CIR lookups and better name filtering.

## [0.4.4] - 2026-02-18

### Changed

- Faster SMILES retrieval from CIR and faster search for the Java executable.

## [0.4.3] - 2026-02-16

### Fixed

- Issue on macOS (#33).

## [0.4.2] - 2026-02-12

### Added

- Ability to save results.
- Clear error message if Java is not found.

### Changed

- Bumped `pillow` from 11.1.0 to 12.1.1.

## [0.4.1] - 2026-01-23

### Fixed

- Rare CompTox issue.

## [0.4.0] - 2026-01-09

### Changed

- Bumped dependencies `urllib3` to 2.6.3, `requests` to 2.32.4, `starlette` to 0.49.1
  and `h11` to 0.16.0.

## [0.3.9] - 2025-09-12

### Changed

- Improved searching for salts.
- Better name cleaning.

## [0.3.8] - 2025-07-11

### Fixed

- PubChem SMILES retrieval.

## [0.3.7] - 2025-07-07

### Added

- `MoleculeResolver` can be instantiated without a `with` block.

## [0.3.6] - 2025-07-07

### Added

- Proper user agent for web requests.

### Fixed

- Getting salts.

## [0.3.5] - 2025-05-15

### Changed

- New OPSIN URL.
- Sulfinyl groups are normalized using SMARTS.

### Fixed

- Extra normalization issues from RDKit in SMILES normalization.
- Default value in `to_dict`.

## [0.3.4] - 2025-03-28

### Added

- Support for newer Python versions.

### Changed

- Small improvements.

## [0.3.3] - 2025-01-22

### Added

- Documentation website.

### Changed

- Simpler and better algorithm.
- Updated dependencies.

## [0.3.2] - 2024-12-31

### Added

- Bundled offline OPSIN executable.

## [0.3.1] - 2024-12-28

### Added

- GitHub link on PyPI.

## [0.3.0] - 2024-12-28

### Added

- `expand_CAS_heuristically`.

### Changed

- Much improved ChEBI support.
- The SQL molecule cache is consistent with web results.

## [0.2.4] - 2024-11-15

### Changed

- Improved chemical name cleaning.

## [0.2.3] - 2024-11-09

### Fixed

- Argument renaming.

## [0.2.2] - 2024-11-03

### Added

- Published to PyPI.

## [0.2.1] - 2024-10-25

[Unreleased]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.6.0...HEAD
[0.6.0]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.5.4...0.6.0
[0.5.4]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.5.3...0.5.4
[0.5.3]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.5.2...0.5.3
[0.5.2]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.5.1...0.5.2
[0.5.1]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.5.0...0.5.1
[0.5.0]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.4.4...0.5.0
[0.4.4]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.4.3...0.4.4
[0.4.3]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.4.2...0.4.3
[0.4.2]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.4.1...0.4.2
[0.4.1]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.4.0...0.4.1
[0.4.0]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.9...0.4.0
[0.3.9]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.8...0.3.9
[0.3.8]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.7...0.3.8
[0.3.7]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.6...0.3.7
[0.3.6]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.5...0.3.6
[0.3.5]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.4...0.3.5
[0.3.4]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.3...0.3.4
[0.3.3]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.2...0.3.3
[0.3.2]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.1...0.3.2
[0.3.1]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.2.4...0.3.0
[0.2.4]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.2.3...0.2.4
[0.2.3]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.2.2...0.2.3
[0.2.2]: https://github.com/MoleculeResolver/molecule-resolver/compare/0.2.1...0.2.2
[0.2.1]: https://github.com/MoleculeResolver/molecule-resolver/releases/tag/0.2.1
