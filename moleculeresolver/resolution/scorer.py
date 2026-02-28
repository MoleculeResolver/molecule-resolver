from __future__ import annotations

from moleculeresolver.resolution.models import (
    StructureGroupCandidate,
    WeightedStructureScore,
)


def score_structure_groups(
    candidates: list[StructureGroupCandidate],
) -> list[WeightedStructureScore]:
    """Score grouped structures with explicit weighted components."""
    scored = []
    for candidate in candidates:
        crosscheck_points = candidate.crosscheck_count * 100
        opsin_bonus = 30 if candidate.supports_opsin_name_match else 0
        # Prefer richer isomer-specific structures for matching non-isomeric groups.
        isomer_specificity_bonus = len(candidate.smiles)
        smiles_length_bonus = len(candidate.smiles)

        scored.append(
            WeightedStructureScore(
                smiles=candidate.smiles,
                total_score=(
                    crosscheck_points
                    + opsin_bonus
                    + isomer_specificity_bonus
                    + smiles_length_bonus
                ),
                crosscheck_count=candidate.crosscheck_count,
                opsin_bonus=opsin_bonus,
                isomer_specificity_bonus=isomer_specificity_bonus,
                smiles_length_bonus=smiles_length_bonus,
            )
        )

    return scored
