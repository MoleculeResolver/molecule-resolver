from __future__ import annotations

from dataclasses import dataclass

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
        return smiles

    @staticmethod
    def to_SMILES(smiles, isomeric=False):
        return smiles.replace("@", "") if not isomeric else smiles


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
    # Same molecule count and bonus signals. Lexicographical SMILES decides.
    scored = [
        type(
            "Score",
            (),
            {
                "smiles": "CCN",
                "total_score": 200,
                "crosscheck_count": 2,
                "opsin_bonus": 0,
                "isomer_specificity_bonus": 3,
                "smiles_length_bonus": 3,
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
                "isomer_specificity_bonus": 3,
                "smiles_length_bonus": 3,
            },
        )(),
    ]

    best = select_best_scored_structure(scored)

    assert best.smiles == "CCN"
