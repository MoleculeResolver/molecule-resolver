from __future__ import annotations

from typing import Optional

from moleculeresolver.services.base import ServiceAdapter, ServiceSearchResult


class PubChemServiceAdapter(ServiceAdapter):
    name = "pubchem"

    def resolve_one(
        self,
        resolver,
        identifier: str,
        mode: str,
        required_formula: Optional[str],
        required_charge: Optional[int],
        required_structure_type: Optional[str],
    ) -> Optional[ServiceSearchResult]:
        if mode not in resolver.supported_modes_by_services[self.name]:
            return None
        cmp = resolver.get_molecule_from_pubchem(
            identifier,
            mode,
            required_formula,
            required_charge,
            required_structure_type,
        )
        if cmp is None:
            return None
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information=f"{cmp.service} id: {cmp.additional_information}",
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )
