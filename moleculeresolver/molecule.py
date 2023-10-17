from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Molecule:
    SMILES: Optional[str] = None
    synonyms: list[str] = field(default_factory=list)
    CAS: list[str] = field(default_factory=list)
    additional_information: str = ""
    mode: str = ""
    service: str = ""
    number_of_crosschecks: int = 1
    identifier: str = ""
    found_molecules: list = field(default_factory=list)
