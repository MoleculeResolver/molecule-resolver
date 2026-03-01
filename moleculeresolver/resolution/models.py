from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
from rdkit import Chem
from rdkit.Chem.rdchem import BondStereo

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
    chiral_center_count: int
    defined_chiral_center_count: int
    bond_stereo_count: int
    stereo_signal_count: int


@dataclass
class WeightedStructureScore:
    """Explicit weighted score and breakdown for one structure group."""

    smiles: str
    total_score: int
    crosscheck_count: int
    opsin_bonus: int
    stereo_specificity_bonus: int
    defined_chirality_bonus: int
    bond_stereo_bonus: int

    # Backward-compatible aliases for earlier score fields.
    @property
    def isomer_specificity_bonus(self) -> int:
        return self.stereo_specificity_bonus

    @property
    def smiles_length_bonus(self) -> int:
        return 0


def build_structure_group_candidates(
    resolver: "MoleculeResolver",
    grouped_molecules: dict[str, list["Molecule"]],
) -> list[StructureGroupCandidate]:
    """Build normalized candidates from grouped molecules for scoring."""
    candidates = []
    for smiles, molecules in grouped_molecules.items():
        mol = resolver.get_from_SMILES(smiles)
        non_isomeric = (
            resolver.to_SMILES(mol, isomeric=False) if mol is not None else smiles
        )
        chiral_center_count = 0
        defined_chiral_center_count = 0
        bond_stereo_count = 0
        if mol is not None:
            chiral_centers = Chem.FindMolChiralCenters(
                mol, includeUnassigned=True, useLegacyImplementation=False
            )
            chiral_center_count = len(chiral_centers)
            defined_chiral_center_count = sum(
                1 for _, label in chiral_centers if label != "?"
            )
            bond_stereo_count = sum(
                1
                for bond in mol.GetBonds()
                if bond.GetStereo() not in {BondStereo.STEREONONE, BondStereo.STEREOANY}
            )
        stereo_signal_count = defined_chiral_center_count + bond_stereo_count
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
                chiral_center_count=chiral_center_count,
                defined_chiral_center_count=defined_chiral_center_count,
                bond_stereo_count=bond_stereo_count,
                stereo_signal_count=stereo_signal_count,
            )
        )
    return candidates
