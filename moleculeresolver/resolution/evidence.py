from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from moleculeresolver.molecule import Molecule


@dataclass
class CandidateEvidence:
    """Evidence and scoring payload for one candidate structure group."""

    smiles: str
    service_agreement_count: int
    service_names: list[str]
    identifiers: list[str]
    identifier_concordance_count: int
    synonym_overlap_count: int = 0
    score_breakdown: dict[str, int] = field(default_factory=dict)
    total_score: int = 0


@dataclass
class ResolutionResult:
    """Extended result payload for include_evidence=True requests."""

    best_molecule: "Molecule | None"
    ranked_candidates: list[CandidateEvidence]
    grouped_by_structure: dict[str, list["Molecule"]]
    selected_smiles: str | None
    selection_reason: str
