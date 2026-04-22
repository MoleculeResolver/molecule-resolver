from __future__ import annotations

from moleculeresolver.resolution.models import WeightedStructureScore


def select_best_scored_structure(
    scored_structures: list[WeightedStructureScore],
) -> WeightedStructureScore:
    """Deterministically select the strongest structure candidate."""
    if not scored_structures:
        raise ValueError("select_best_scored_structure requires a non-empty list.")

    return sorted(
        scored_structures,
        key=lambda score: (
            -score.total_score,
            -score.crosscheck_count,
            -score.opsin_bonus,
            -score.stereo_specificity_bonus,
            -score.defined_chirality_bonus,
            -score.bond_stereo_bonus,
            score.smiles,
        ),
    )[0]
