from __future__ import annotations

import inspect

import pytest

from moleculeresolver import MoleculeResolver


def test_multi_name_exhaustive_benchmark_case(monkeypatch):
    if "search_strategy" not in inspect.signature(
        MoleculeResolver.find_single_molecule
    ).parameters:
        pytest.skip("search_strategy is not available in this branch state")

    mr = MoleculeResolver()
    mr.supported_modes_by_services = {"svc1": ["name"], "svc2": ["name"]}
    mr._available_services = ["svc1", "svc2"]

    monkeypatch.setattr(mr, "_check_parameters", lambda **_: None)
    monkeypatch.setattr(
        mr,
        "_check_and_flatten_identifiers_and_modes",
        lambda identifiers, modes: (
            ["isopropanol", "2-propanol"],
            ["name", "name"],
            ["isopropanol", "2-propanol"],
            set(),
            None,
        ),
    )
    monkeypatch.setattr(mr, "standardize_SMILES", lambda smiles: smiles)

    result_map = {
        ("svc1", "isopropanol"): {
            "SMILES": "CC(O)C",
            "synonyms": ["isopropanol"],
            "CAS": {"67-63-0"},
            "additional_information": "svc1",
            "mode_used": "name",
            "identifier_used": "isopropanol",
            "service": "svc1",
            "cas_is_authoritative": False,
        },
        ("svc1", "2-propanol"): {
            "SMILES": "CC(C)O",
            "synonyms": ["2-propanol"],
            "CAS": {"67-63-0"},
            "additional_information": "svc1",
            "mode_used": "name",
            "identifier_used": "2-propanol",
            "service": "svc1",
            "cas_is_authoritative": False,
        },
        ("svc2", "isopropanol"): {
            "SMILES": "CC(C)O",
            "synonyms": ["isopropanol"],
            "CAS": {"67-63-0"},
            "additional_information": "svc2",
            "mode_used": "name",
            "identifier_used": "isopropanol",
            "service": "svc2",
            "cas_is_authoritative": False,
        },
    }

    monkeypatch.setattr(
        mr,
        "_resolve_single_service_candidate",
        lambda service, identifier, mode, *_: result_map.get((service, identifier)),
    )

    result = mr.find_single_molecule(
        identifiers=["isopropanol", "2-propanol"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        search_strategy="exhaustive",
    )

    assert result is not None
    assert result.SMILES == "CC(C)O"


def test_consensus_vs_legacy_mode_comparison(monkeypatch):
    if "resolution_mode" not in inspect.signature(
        MoleculeResolver.find_single_molecule_crosschecked
    ).parameters:
        pytest.skip("resolution_mode is not available in this branch state")

    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    from moleculeresolver.molecule import Molecule

    legacy_winner = Molecule("CCO", ["ethanol"], [], "svc1", "name", "svc1", 1, "ethanol")
    consensus_winner = Molecule("CCO", ["ethyl alcohol"], [], "svc2", "name", "svc2", 1, "ethyl alcohol")

    def fake_find_single_molecule(*args, **kwargs):
        service = kwargs["services_to_use"][0]
        if service == "svc1":
            return legacy_winner
        return consensus_winner

    monkeypatch.setattr(mr, "find_single_molecule", fake_find_single_molecule)
    monkeypatch.setattr(
        mr,
        "filter_molecules",
        lambda molecules, *_: [molecule for molecule in molecules if molecule is not None],
    )
    monkeypatch.setattr(
        mr,
        "group_molecules_by_structure",
        lambda molecules, *_: {"CCO": molecules},
    )

    legacy_result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol"],
        modes=["name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="legacy",
    )
    consensus_result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol"],
        modes=["name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="consensus",
    )

    assert legacy_result is not None
    assert consensus_result is not None
    assert legacy_result.SMILES == consensus_result.SMILES


def test_strict_isomer_acceptance_case(monkeypatch):
    if "resolution_mode" not in inspect.signature(
        MoleculeResolver.find_single_molecule_crosschecked
    ).parameters:
        pytest.skip("resolution_mode is not available in this branch state")
    if not hasattr(MoleculeResolver, "_collect_opsin_isomer_matches"):
        pytest.skip("strict isomer verification helpers are not available in this branch state")

    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    from moleculeresolver.molecule import Molecule

    mol_a = Molecule("CCO", ["ethanol"], [], "svc1", "name", "svc1", 1, "ethanol")
    mol_b = Molecule("C[C@H](O)C", ["(S)-2-butanol"], [], "svc2", "name", "svc2", 1, "(S)-2-butanol")

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
    monkeypatch.setattr(
        mr,
        "_collect_opsin_isomer_matches",
        lambda grouped, smiles: {smiles[0]: False, smiles[1]: True},
    )

    result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol", "(S)-2-butanol"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="strict_isomer",
    )

    assert result is not None
    assert result.SMILES == "C[C@H](O)C"
