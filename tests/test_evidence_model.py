from __future__ import annotations

from moleculeresolver import MoleculeResolver
from moleculeresolver.molecule import Molecule
from moleculeresolver.resolution import ResolutionResult


def test_include_evidence_returns_resolution_result(monkeypatch):
    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    mol_a = Molecule("CCO", ["Ethanol", "EtOH"], ["64-17-5"], "svc1", "name", "svc1", 1, "a")
    mol_b = Molecule("CCO", ["ethanol"], [], "svc2", "name", "svc2", 1, "b")

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
        mr, "group_molecules_by_structure", lambda molecules, *_: {"CCO": [mol_a, mol_b]}
    )

    result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol"],
        modes=["name"],
        services_to_use=["svc1", "svc2"],
        include_evidence=True,
    )

    assert isinstance(result, ResolutionResult)
    assert result.best_molecule is not None
    assert result.selected_smiles == "CCO"
    assert len(result.ranked_candidates) == 1
    evidence = result.ranked_candidates[0]
    assert evidence.service_agreement_count == 2
    assert evidence.identifier_concordance_count == 2
    assert evidence.synonym_overlap_count >= 1
    assert "service_agreement" in evidence.score_breakdown
    assert evidence.total_score == sum(evidence.score_breakdown.values())


def test_include_evidence_with_tied_structures_has_no_best_molecule(monkeypatch):
    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    mol_a = Molecule("CCO", ["ethanol"], [], "svc1", "name", "svc1", 1, "a")
    mol_b = Molecule("CCC", ["propane"], [], "svc2", "name", "svc2", 1, "b")

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
        lambda molecules, *_: {"CCO": [mol_a], "CCC": [mol_b]},
    )

    result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol", "propane"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        include_evidence=True,
        try_to_choose_best_structure=False,
    )

    assert isinstance(result, ResolutionResult)
    assert result.best_molecule is None
    assert result.selection_reason == "multiple_structures_tied"
