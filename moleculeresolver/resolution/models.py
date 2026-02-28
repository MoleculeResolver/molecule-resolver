from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from moleculeresolver.molecule import Molecule
    from moleculeresolver.moleculeresolver import MoleculeResolver


@dataclass
class StructureGroupCandidate:
    """Single grouped-structure candidate used in consensus scoring."""

    smiles: str
    molecules: list["Molecule"]
    crosscheck_count: int
    non_isomeric_smiles: str
    supports_opsin_name_match: bool


@dataclass
class WeightedStructureScore:
    """Explicit weighted score and breakdown for one structure group."""

    smiles: str
    total_score: int
    crosscheck_count: int
    opsin_bonus: int
    isomer_specificity_bonus: int
    smiles_length_bonus: int


def build_structure_group_candidates(
    resolver: "MoleculeResolver",
    grouped_molecules: dict[str, list["Molecule"]],
) -> list[StructureGroupCandidate]:
    """Build normalized candidates from grouped molecules for scoring."""
    candidates = []
    for smiles, molecules in grouped_molecules.items():
        non_isomeric = resolver.to_SMILES(resolver.get_from_SMILES(smiles), isomeric=False)
        candidates.append(
            StructureGroupCandidate(
                smiles=smiles,
                molecules=molecules,
                crosscheck_count=len(molecules),
                non_isomeric_smiles=non_isomeric,
                supports_opsin_name_match=any(
                    molecule.mode == "name" and molecule.service == "opsin"
                    for molecule in molecules
                ),
            )
        )
    return candidates
