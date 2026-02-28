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
