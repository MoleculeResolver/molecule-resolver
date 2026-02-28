from __future__ import annotations

from typing import Optional

from moleculeresolver.services.base import ServiceAdapter, ServiceSearchResult


class OPSINServiceAdapter(ServiceAdapter):
    name = "opsin"

    def resolve(
        self,
        resolver,
        flattened_identifiers: list[str],
        flattened_modes: list[str],
        required_formula: Optional[str],
        required_charge: Optional[int],
        required_structure_type: Optional[str],
    ) -> Optional[ServiceSearchResult]:
        for identifier, mode in zip(flattened_identifiers, flattened_modes, strict=True):
            if mode not in resolver.supported_modes_by_services[self.name]:
                continue
            cmp = resolver.get_molecule_from_OPSIN(
                identifier,
                required_formula,
                required_charge,
                required_structure_type,
            )
            if cmp is not None:
                return ServiceSearchResult(
                    molecule=cmp,
                    mode_used=mode,
                    identifier_used=cmp.identifier,
                    additional_information=cmp.additional_information,
                    current_service=self.name,
                    synonyms=list(cmp.synonyms),
                    cas=set(),
                )
        return None
