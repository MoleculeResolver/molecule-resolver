from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import copy


@dataclass
class Molecule:
    """
    Represents a molecule with various properties and identifiers.

    Attributes:
        SMILES (Optional[str]): The SMILES (Simplified Molecular Input Line Entry System) representation of the molecule.

        synonyms (Optional[list[str]]): A list of alternative names or synonyms for the molecule.

        CAS (Optional[list[str]]): A list of CAS (Chemical Abstracts Service) registry numbers for the molecule.

        additional_information (Optional[str]): Any additional information about the molecule.

        mode (Optional[str]): The mode associated with the molecule.
        
        service (Optional[str]): The service associated with the molecule.

        number_of_crosschecks (Optional[int]): The number of cross-checks performed on the molecule.

        identifier (Optional[str]): A unique identifier for the molecule.

        found_molecules (Optional[list]): A list of related molecules found during processing.
    """

    SMILES: Optional[str] = None
    synonyms: Optional[List[str]] = field(default_factory=list)
    CAS: Optional[List[str]] = field(default_factory=list)
    additional_information: Optional[str] = ""
    mode: Optional[str] = ""
    service: Optional[str] = ""
    number_of_crosschecks: Optional[int] = 1
    identifier: Optional[str] = ""
    found_molecules: Optional[list] = field(default_factory=list)

    def to_dict(self, found_molecules: Optional[str] = None) -> Dict[str, Any]:
        """
        Convert the Molecule object to a dictionary.

        Args:
            found_molecules (Optional[str]): Determines how 'found_molecules' are handled.
                - If 'remove', the 'found_molecules' field will be excluded.
                - If 'recursive', 'found_molecules' will be recursively converted to dictionaries.
                - If None, 'found_molecules' will be included as is.

        Returns:
            Dict[str, Any]: A dictionary representation of the Molecule object.

        Note:
            This method creates a deep copy of the object's `__dict__` attribute.
            Depending on the `found_molecules` parameter, it may exclude or recursively convert
            the 'found_molecules' field before returning the dictionary.
        """
        d = copy.deepcopy(self.__dict__)
        if found_molecules == "remove":
            if "found_molecules" in d:
                d.pop("found_molecules")
        elif found_molecules == "recursive":
            if "found_molecules" in d:
                new_found_molecules = []
                for grouped_item in d["found_molecules"]:
                    key = list(grouped_item.keys())[0]
                    value = list(grouped_item.values())[0]
                    new_found_molecules.append(
                        {key: [m.to_dict("recursive") for m in value]}
                    )
                d["found_molecules"] = new_found_molecules
        return d
