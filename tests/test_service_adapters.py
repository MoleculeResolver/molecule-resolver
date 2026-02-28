from __future__ import annotations

from dataclasses import dataclass

from moleculeresolver.services import (
    CASRegistryServiceAdapter,
    OPSINServiceAdapter,
    PubChemServiceAdapter,
    ServiceAdapterRegistry,
)


@dataclass
class _FakeMolecule:
    SMILES: str
    synonyms: list[str]
    CAS: list[str]
    service: str
    mode: str
    identifier: str
    additional_information: str


class _FakeResolver:
    def __init__(self):
        self.supported_modes_by_services = {
            "opsin": ["name"],
            "pubchem": ["name", "cas"],
            "cas_registry": ["name", "cas"],
        }

    def get_molecule_from_OPSIN(
        self,
        identifier,
        required_formula,
        required_charge,
        required_structure_type,
    ):
        if identifier == "ethanol":
            return _FakeMolecule(
                "CCO",
                ["ethanol"],
                [],
                "opsin",
                "name",
                identifier,
                "opsin resolution",
            )
        return None

    def get_molecule_from_pubchem(
        self,
        identifier,
        mode,
        required_formula,
        required_charge,
        required_structure_type,
    ):
        if identifier == "ethanol":
            return _FakeMolecule(
                "CCO",
                ["ethyl alcohol"],
                ["64-17-5"],
                "pubchem",
                mode,
                identifier,
                "702",
            )
        return None

    def get_molecule_from_CAS_registry(
        self,
        identifier,
        mode,
        required_formula,
        required_charge,
        required_structure_type,
    ):
        if identifier == "64-17-5":
            return _FakeMolecule(
                "CCO",
                ["ethanol"],
                ["64-17-5"],
                "cas_registry",
                mode,
                identifier,
                "cas_registry",
            )
        return None


def test_service_adapter_registry_roundtrip():
    registry = ServiceAdapterRegistry()
    adapter = OPSINServiceAdapter()
    registry.register(adapter)

    assert registry.get("opsin") is adapter
    assert registry.names() == ["opsin"]


def test_opsin_adapter_resolves_name_only():
    resolver = _FakeResolver()
    adapter = OPSINServiceAdapter()

    result = adapter.resolve(
        resolver,
        flattened_identifiers=["ethanol"],
        flattened_modes=["name"],
        required_formula=None,
        required_charge=None,
        required_structure_type=None,
    )
    assert result is not None
    assert result.molecule.SMILES == "CCO"
    assert result.mode_used == "name"
    assert result.current_service == "opsin"


def test_pubchem_adapter_formats_additional_information():
    resolver = _FakeResolver()
    adapter = PubChemServiceAdapter()

    result = adapter.resolve(
        resolver,
        flattened_identifiers=["ethanol"],
        flattened_modes=["name"],
        required_formula=None,
        required_charge=None,
        required_structure_type=None,
    )
    assert result is not None
    assert result.additional_information == "pubchem id: 702"
    assert "64-17-5" in result.cas


def test_cas_registry_adapter_provides_authoritative_cas_set():
    resolver = _FakeResolver()
    adapter = CASRegistryServiceAdapter()

    result = adapter.resolve(
        resolver,
        flattened_identifiers=["64-17-5"],
        flattened_modes=["cas"],
        required_formula=None,
        required_charge=None,
        required_structure_type=None,
    )
    assert result is not None
    assert result.current_service == "cas_registry"
    assert result.cas == {"64-17-5"}
