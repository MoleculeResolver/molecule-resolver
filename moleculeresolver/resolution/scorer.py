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
        unresolved_chiral_centers = max(
            0, candidate.chiral_center_count - candidate.defined_chiral_center_count
        )
        defined_chirality_bonus = candidate.defined_chiral_center_count * 20
        bond_stereo_bonus = candidate.bond_stereo_count * 15
        stereo_specificity_bonus = (
            defined_chirality_bonus + bond_stereo_bonus + unresolved_chiral_centers * 5
        )

        scored.append(
            WeightedStructureScore(
                smiles=candidate.smiles,
                total_score=(
                    crosscheck_points
                    + opsin_bonus
                    + stereo_specificity_bonus
                ),
                crosscheck_count=candidate.crosscheck_count,
                opsin_bonus=opsin_bonus,
                stereo_specificity_bonus=stereo_specificity_bonus,
                defined_chirality_bonus=defined_chirality_bonus,
                bond_stereo_bonus=bond_stereo_bonus,
            )
        )

    return scored
