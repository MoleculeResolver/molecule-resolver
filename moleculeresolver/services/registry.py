from __future__ import annotations

from typing import Optional

from moleculeresolver.services.base import ServiceAdapter


class ServiceAdapterRegistry:
    """Simple runtime registry for all service adapters."""

    def __init__(self):
        self._adapters: dict[str, ServiceAdapter] = {}

    def register(self, adapter: ServiceAdapter) -> None:
        self._adapters[adapter.name] = adapter

    def get(self, service_name: str) -> Optional[ServiceAdapter]:
        return self._adapters.get(service_name)

    def names(self) -> list[str]:
        return sorted(self._adapters)
