from __future__ import annotations

from typing import Optional

from moleculeresolver.services.base import ServiceAdapter, ServiceSearchResult


class _ModeDrivenAdapter(ServiceAdapter):
    fetch_method_name: str

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        raise NotImplementedError

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
        fetch_method = getattr(resolver, self.fetch_method_name)
        cmp = fetch_method(
            identifier,
            mode,
            required_formula,
            required_charge,
            required_structure_type,
        )
        if cmp is not None:
            return self._build_result(cmp, mode)
        return None


class CASRegistryServiceAdapter(_ModeDrivenAdapter):
    name = "cas_registry"
    fetch_method_name = "get_molecule_from_CAS_registry"

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information=cmp.service,
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )


class ChEBIServiceAdapter(_ModeDrivenAdapter):
    name = "chebi"
    fetch_method_name = "get_molecule_from_ChEBI"

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information=f"{cmp.service} id: {cmp.additional_information}",
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )


class SRSServiceAdapter(_ModeDrivenAdapter):
    name = "srs"
    fetch_method_name = "get_molecule_from_SRS"

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information=f"{cmp.service} id: {cmp.additional_information}",
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )


class CompToxServiceAdapter(_ModeDrivenAdapter):
    name = "comptox"
    fetch_method_name = "get_molecule_from_CompTox"

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information=f"{cmp.service} id: {cmp.additional_information}",
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )


class ChemeoServiceAdapter(_ModeDrivenAdapter):
    name = "chemeo"
    fetch_method_name = "get_molecule_from_Chemeo"

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information=f"{cmp.service} id: {cmp.additional_information}",
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )


class CTSServiceAdapter(_ModeDrivenAdapter):
    name = "cts"
    fetch_method_name = "get_molecule_from_CTS"

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information="cts",
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )


class NISTServiceAdapter(_ModeDrivenAdapter):
    name = "nist"
    fetch_method_name = "get_molecule_from_NIST"

    def _build_result(self, cmp, mode: str) -> ServiceSearchResult:
        return ServiceSearchResult(
            molecule=cmp,
            mode_used=cmp.mode,
            identifier_used=cmp.identifier,
            additional_information=f"{cmp.service} id: {cmp.additional_information}",
            current_service=self.name,
            synonyms=list(cmp.synonyms),
            cas=set(cmp.CAS),
        )
