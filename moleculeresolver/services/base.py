from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from moleculeresolver.molecule import Molecule
    from moleculeresolver.moleculeresolver import MoleculeResolver


@dataclass
class ServiceSearchResult:
    """Normalized single-service search result consumed by MoleculeResolver."""

    molecule: "Molecule"
    mode_used: str
    identifier_used: str
    additional_information: Optional[str]
    current_service: str
    synonyms: list[str]
    cas: set[str]


class ServiceAdapter(ABC):
    """Contract for a resolver service adapter."""

    name: str

    @abstractmethod
    def resolve(
        self,
        resolver: "MoleculeResolver",
        flattened_identifiers: list[str],
        flattened_modes: list[str],
        required_formula: Optional[str],
        required_charge: Optional[int],
        required_structure_type: Optional[str],
    ) -> Optional[ServiceSearchResult]:
        """Resolve a molecule for this adapter or return None."""
        raise NotImplementedError
