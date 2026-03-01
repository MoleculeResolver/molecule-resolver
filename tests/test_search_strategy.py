from __future__ import annotations

from types import SimpleNamespace

from moleculeresolver import MoleculeResolver


def _make_exhaustive_adapter_result(
    smiles: str,
    synonym: str,
    *,
    identifier: str,
    mode: str = "name",
):
    return SimpleNamespace(
        molecule=SimpleNamespace(SMILES=smiles),
        synonyms=[synonym],
        cas=set(),
        additional_information=f"{identifier}:{mode}",
        mode_used=mode,
        identifier_used=identifier,
    )


def _make_strategy_test_resolver(monkeypatch) -> MoleculeResolver:
    mr = MoleculeResolver()
    mr.supported_modes_by_services = {
        "svc1": ["name"],
        "svc2": ["name"],
    }
    mr._available_services = ["svc1", "svc2"]

    monkeypatch.setattr(mr, "_check_parameters", lambda **_: None)
    monkeypatch.setattr(
        mr,
        "_check_and_flatten_identifiers_and_modes",
        lambda identifiers, modes: (
            ["alpha", "beta"],
            ["name", "name"],
            ["alpha", "beta"],
            set(),
            None,
        ),
    )
    monkeypatch.setattr(mr, "standardize_SMILES", lambda smi: smi)
    return mr


def test_first_hit_stops_after_first_match(monkeypatch):
    mr = _make_strategy_test_resolver(monkeypatch)
    calls = []

    def fake_resolve_with_adapter(service, identifiers, modes, *_):
        calls.append((service, tuple(identifiers), tuple(modes)))
        if service == "svc1":
            return SimpleNamespace(
                molecule=SimpleNamespace(SMILES="CCO"),
                synonyms=["alpha"],
                additional_information="svc1:alpha",
                mode_used="name",
                identifier_used="alpha",
                cas=set(),
            )
        return None

    monkeypatch.setattr(mr, "_resolve_service_with_adapter", fake_resolve_with_adapter)

    result = mr.find_single_molecule(
        identifiers=["alpha", "beta"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        search_strategy="first_hit",
    )

    assert result is not None
    assert result.SMILES == "CCO"
    assert calls == [("svc1", ("alpha", "beta"), ("name", "name"))]


def test_exhaustive_search_checks_all_pairs_and_uses_consensus(monkeypatch):
    mr = _make_strategy_test_resolver(monkeypatch)
    calls = []
    response_map = {
        ("svc1", "alpha"): _make_exhaustive_adapter_result("CCO", "alpha", identifier="alpha"),
        ("svc1", "beta"): _make_exhaustive_adapter_result("CCC", "beta", identifier="beta"),
        ("svc2", "alpha"): _make_exhaustive_adapter_result("CCO", "alpha", identifier="alpha"),
        ("svc2", "beta"): _make_exhaustive_adapter_result("CCO", "beta", identifier="beta"),
    }

    def fake_resolve(service, identifier, mode, *_):
        calls.append((service, identifier, mode))
        return response_map.get((service, identifier))

    monkeypatch.setattr(mr, "_resolve_identifier_with_adapter", fake_resolve)

    result = mr.find_single_molecule(
        identifiers=["alpha", "beta"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        search_strategy="exhaustive",
    )

    assert result is not None
    assert result.SMILES == "CCO"
    assert len(calls) == 4
    assert set(result.synonyms) == {"alpha", "beta"}
