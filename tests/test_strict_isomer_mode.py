from __future__ import annotations

from moleculeresolver import MoleculeResolver
from moleculeresolver.molecule import Molecule


def _build_stubbed_crosscheck_resolver(monkeypatch):
    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    mol_a = Molecule("CCO", ["alpha"], [], "svc1", "name", "svc1", 1, "alpha")
    mol_b = Molecule("C[C@H](O)C", ["beta"], [], "svc2", "name", "svc2", 1, "beta")

    def fake_find_single_molecule(*args, **kwargs):
        service = kwargs["services_to_use"][0]
        return mol_a if service == "svc1" else mol_b

    monkeypatch.setattr(mr, "find_single_molecule", fake_find_single_molecule)
    monkeypatch.setattr(
        mr,
        "filter_molecules",
        lambda molecules, *_: [molecule for molecule in molecules if molecule is not None],
    )
    monkeypatch.setattr(
        mr,
        "group_molecules_by_structure",
        lambda molecules, *_: {"CCO": [mol_a], "C[C@H](O)C": [mol_b]},
    )
    return mr


def test_strict_isomer_rejects_when_no_opsin_match(monkeypatch):
    mr = _build_stubbed_crosscheck_resolver(monkeypatch)
    monkeypatch.setattr(
        mr,
        "_collect_opsin_isomer_matches",
        lambda grouped, smiles: {smiles[0]: False, smiles[1]: False},
    )

    result = mr.find_single_molecule_crosschecked(
        identifiers=["alpha", "beta"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="strict_isomer",
    )

    assert result is None


def test_strict_isomer_returns_verified_candidate(monkeypatch):
    mr = _build_stubbed_crosscheck_resolver(monkeypatch)
    monkeypatch.setattr(
        mr,
        "_collect_opsin_isomer_matches",
        lambda grouped, smiles: {smiles[0]: False, smiles[1]: True},
    )

    result = mr.find_single_molecule_crosschecked(
        identifiers=["alpha", "beta"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="strict_isomer",
    )

    assert result is not None
    assert result.SMILES == "C[C@H](O)C"


def test_legacy_mode_skips_strict_isomer_filter(monkeypatch):
    mr = _build_stubbed_crosscheck_resolver(monkeypatch)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("strict OPSIN filtering should not run in legacy mode")

    monkeypatch.setattr(mr, "_collect_opsin_isomer_matches", fail_if_called)

    result = mr.find_single_molecule_crosschecked(
        identifiers=["alpha", "beta"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="legacy",
    )

    assert result is not None
