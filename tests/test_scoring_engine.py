from __future__ import annotations

from dataclasses import dataclass
from rdkit import Chem

from moleculeresolver.resolution import (
    build_structure_group_candidates,
    score_structure_groups,
    select_best_scored_structure,
)


@dataclass
class _FakeMolecule:
    mode: str
    service: str


class _FakeResolver:
    @staticmethod
    def get_from_SMILES(smiles):
        return Chem.MolFromSmiles(smiles)

    @staticmethod
    def to_SMILES(smiles, isomeric=False):
        return Chem.MolToSmiles(smiles, isomericSmiles=isomeric)


def test_weighted_scoring_prefers_higher_crosscheck_count():
    resolver = _FakeResolver()
    grouped = {
        "CCO": [_FakeMolecule(mode="name", service="opsin"), _FakeMolecule(mode="name", service="pubchem")],
        "CCC": [_FakeMolecule(mode="name", service="pubchem")],
    }

    candidates = build_structure_group_candidates(resolver, grouped)
    scored = score_structure_groups(candidates)
    best = select_best_scored_structure(scored)

    assert best.smiles == "CCO"
    assert best.crosscheck_count == 2


def test_tie_breaker_is_deterministic_for_equal_scores():
    # Same score components. Lexicographical SMILES decides.
    scored = [
        type(
            "Score",
            (),
            {
                "smiles": "CCN",
                "total_score": 200,
                "crosscheck_count": 2,
                "opsin_bonus": 0,
                "stereo_specificity_bonus": 10,
                "defined_chirality_bonus": 10,
                "bond_stereo_bonus": 0,
            },
        )(),
        type(
            "Score",
            (),
            {
                "smiles": "CCO",
                "total_score": 200,
                "crosscheck_count": 2,
                "opsin_bonus": 0,
                "stereo_specificity_bonus": 10,
                "defined_chirality_bonus": 10,
                "bond_stereo_bonus": 0,
            },
        )(),
    ]

    best = select_best_scored_structure(scored)

    assert best.smiles == "CCN"


def test_stereo_specific_candidate_wins_when_crosscheck_is_equal():
    resolver = _FakeResolver()
    grouped = {
        "CCO": [_FakeMolecule(mode="name", service="pubchem")],
        "C[C@H](O)F": [_FakeMolecule(mode="name", service="pubchem")],
    }

    candidates = build_structure_group_candidates(resolver, grouped)
    scored = score_structure_groups(candidates)
    best = select_best_scored_structure(scored)

    assert best.smiles == "C[C@H](O)F"
    assert best.stereo_specificity_bonus > 0
