from datetime import datetime
import pytest
from moleculeresolver import Molecule, MoleculeResolver
import json
import os
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Union
from tqdm import tqdm
from dataclasses import dataclass
from types import SimpleNamespace


# IUPAC names
dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
with open(dir_path / "benchmark_component_molecules_iupac.json", "r") as f:
    benchmark = json.load(f)

RESPONSES_PATH = dir_path / "responses.json"
SMILES = "SMILES"

# Set to "1" to query services that have no recorded response and save their answers
# to responses.json. Otherwise tests that would need a live request are skipped.
RECORD_RESPONSES = os.environ.get("MOLECULERESOLVER_RECORD_RESPONSES") == "1"

API_KEYS = {
    "cas_registry": os.environ.get("CAS_REGISTRY_API_KEY"),
    "chemeo": os.environ.get("CHEMEO_API_KEY"),
}


def _response_key(
    url: str,
    kwargs: Optional[dict[str, Any]],
    request_type: Optional[str],
    json_data: Any,
) -> str:
    # GET requests are keyed by URL only, POST requests also by their body, as e.g.
    # PubChem sends different queries to the same URL.
    if request_type != "post":
        return url

    body = json_data if json_data is not None else (kwargs or {}).get("data")
    if isinstance(body, bytes):
        body = body.decode("utf-8")
    elif not isinstance(body, str):
        body = json.dumps(body, sort_keys=True)
    return f"POST {url} {body}"


class MoleculeResolverPatched(MoleculeResolver):
    def __init__(
        self,
        json_path: Optional[str] = None,
        record: Optional[bool] = None,
        available_service_API_keys: Optional[dict[str, Optional[str]]] = None,
        **kwargs,
    ):
        # An in-memory molecule cache makes sure every lookup reaches _resilient_request.
        kwargs.setdefault("molecule_cache_db_path", ":memory:")
        super().__init__(available_service_API_keys=available_service_API_keys, **kwargs)
        self.json_path = json_path
        self.record = RECORD_RESPONSES if record is None else record
        self._json_data_changed = False
        if self.json_path is not None:
            self.json_path = Path(self.json_path)

        if self.json_path is not None and self.json_path.exists():
            with open(self.json_path, "r") as f:
                self.json_data = json.load(f)
        else:
            self.json_data = {}

    def _resilient_request(
        self,
        url: str,
        kwargs: Optional[dict[str, Any]] = None,
        request_type: Optional[str] = "get",
        accepted_status_codes: list[int] = [200],
        rejected_status_codes: list[int] = [404],
        offline_status_codes: list[int] = [],
        max_retries: Optional[int] = 10,
        sleep_time: Union[int, float] = 2,
        allow_redirects: Optional[bool] = False,
        json: Optional[str] = None,
        return_response: Optional[bool] = False,
    ) -> Optional[str]:
        key = _response_key(url, kwargs, request_type, json)
        if self.json_path is not None and not return_response:
            if key in self.json_data:
                return self.json_data[key]
            if not self.record:
                pytest.skip(f"No recorded response for {key[:200]}")

        res = super()._resilient_request(
            url,
            kwargs=kwargs,
            request_type=request_type,
            accepted_status_codes=accepted_status_codes,
            rejected_status_codes=rejected_status_codes,
            offline_status_codes=offline_status_codes,
            max_retries=max_retries,
            sleep_time=sleep_time,
            allow_redirects=allow_redirects,
            json=json,
            return_response=return_response,
        )
        # Raw response objects cannot be stored as JSON.
        if self.json_path is not None and not return_response:
            self.json_data[key] = res
            self._json_data_changed = True
        return res

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.json_path is not None and self._json_data_changed:
            with open(self.json_path, "w") as f:
                json.dump(self.json_data, f)
        return super().__exit__(exc_type, exc_value, exc_traceback)


@pytest.mark.parametrize("data", benchmark.values())
class TestServices:
    """Benchmark of each service, replayed from recorded responses in responses.json.

    A benchmark entry can list services that are known to fail for it in
    "expected_failures", e.g. ["cir", "nist"]; these are marked as xfail.
    """

    json_path: Optional[str] = RESPONSES_PATH

    @staticmethod
    def _test_service(
        call_method: Callable,
        input_identifier: str,
        output_identifier_type: str,
        output_identifier,
        kwargs: Optional[Dict] = None,
    ):
        """
        Test a service by calling it with an input identifier and checking that the output identifier matches the expected value.

        Parameters
        ----------
        call_method : Callable
            The method to call
        input_identifier : str
            The input identifier
        output_identifier_type : str
            The type of the output identifier
        output_identifier : str
            The expected output identifier
        kwargs : Optional[Dict], optional
            Additional keyword arguments to pass to the call method, by default None


        """
        if kwargs is None:
            kwargs = {}
        res = call_method(input_identifier, **kwargs)
        if res is None:
            raise ValueError(f"No molecule found for {input_identifier}")

        res_txt = res.__dict__[output_identifier_type]
        if res_txt == output_identifier:
            return
        else:
            raise ValueError(f"Expected {output_identifier} but got {res_txt}")

    def _test_service_by_name(self, service: str, method_name: str, data: dict, **kwargs):
        if service in data.get("expected_failures", []):
            pytest.xfail(f"{service} is expected to fail for {data['iupac_name']}")
        with MoleculeResolverPatched(json_path=self.json_path, available_service_API_keys=API_KEYS) as mr:
            self._test_service(
                getattr(mr, method_name),
                data["iupac_name"],
                SMILES,
                data["SMILES"],
                kwargs,
            )

    def test_opsin(self, data):
        self._test_service_by_name("opsin", "get_molecule_from_OPSIN", data)

    def test_pubchem(self, data):
        self._test_service_by_name("pubchem", "get_molecule_from_pubchem", data, mode="name")

    def test_cts(self, data):
        # CTS-Lite does not support searching by name or CAS
        if not data.get("inchikey"):
            pytest.skip("No InChIKey in benchmark entry")
        if "cts" in data.get("expected_failures", []):
            pytest.xfail(f"cts is expected to fail for {data['iupac_name']}")
        with MoleculeResolverPatched(json_path=self.json_path) as mr:
            self._test_service(
                mr.get_molecule_from_CTS,
                data["inchikey"],
                SMILES,
                data["SMILES"],
                {"mode": "inchikey"},
            )

    def test_chemeo(self, data):
        if API_KEYS["chemeo"] is None:
            pytest.skip("CHEMEO_API_KEY not set")
        self._test_service_by_name("chemeo", "get_molecule_from_Chemeo", data, mode="name")

    def test_cas(self, data):
        if API_KEYS["cas_registry"] is None:
            pytest.skip("CAS_REGISTRY_API_KEY not set")
        self._test_service_by_name("cas_registry", "get_molecule_from_CAS_registry", data, mode="name")

    def test_chebi(self, data):
        self._test_service_by_name("chebi", "get_molecule_from_ChEBI", data, mode="name")

    def test_cir(self, data):
        self._test_service_by_name("cir", "get_molecule_from_CIR", data, mode="name")

    def test_nist(self, data):
        self._test_service_by_name("nist", "get_molecule_from_NIST", data, mode="name")

    def test_srs(self, data):
        self._test_service_by_name("srs", "get_molecule_from_SRS", data, mode="name")

    # CompTox puts a unique request id into every URL, so its responses cannot be
    # replayed; it is only covered by test_service_online.


@pytest.mark.online
@pytest.mark.parametrize(
    "service, method_name, kwargs",
    [
        ("opsin", "get_molecule_from_OPSIN", {}),
        ("pubchem", "get_molecule_from_pubchem", {"mode": "name"}),
        ("comptox", "get_molecule_from_CompTox", {"mode": "name"}),
        ("chemeo", "get_molecule_from_Chemeo", {"mode": "name"}),
        ("cas_registry", "get_molecule_from_CAS_registry", {"mode": "name"}),
        ("chebi", "get_molecule_from_ChEBI", {"mode": "name"}),
        ("cir", "get_molecule_from_CIR", {"mode": "name"}),
        ("nist", "get_molecule_from_NIST", {"mode": "name"}),
        ("srs", "get_molecule_from_SRS", {"mode": "name"}),
    ],
)
def test_service_online(service, method_name, kwargs):
    """Check that each service is reachable and still finds a well-known compound."""
    if service in API_KEYS and API_KEYS[service] is None:
        pytest.skip(f"No API key for {service}")
    with MoleculeResolver(
        available_service_API_keys=API_KEYS, molecule_cache_db_path=":memory:"
    ) as mr:
        res = getattr(mr, method_name)("ethanol", **kwargs)
    assert res is not None, f"{service} did not find ethanol"
    assert res.SMILES == "CCO"


@pytest.mark.online
def test_cts_online():
    with MoleculeResolver(molecule_cache_db_path=":memory:") as mr:
        res = mr.get_molecule_from_CTS("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", mode="inchikey")
    assert res is not None, "cts did not find ethanol"
    assert res.SMILES == "CCO"


def test_opsin_batchmode():
    names = [d["iupac_name"] for d in benchmark.values()]
    smiles = [d["SMILES"] for d in benchmark.values()]
    with MoleculeResolver(molecule_cache_db_path=":memory:") as mr:
        if not mr._java_path:
            pytest.skip("Java is needed for OPSIN batch mode")
        res = mr.get_molecule_from_OPSIN_batchmode(names)
        for i, r in enumerate(res):
            if r is None:
                raise ValueError(f"No molecule found for {names[i]}")
            # Batch mode returns OPSIN's SMILES without standardization.
            res_SMILES = mr.standardize_SMILES(r[0].SMILES)
            if res_SMILES == smiles[i]:
                continue
            else:
                raise ValueError("Expected " + smiles[i] + " but got " + res_SMILES)


def generate_data(json_path, overwrite_json: bool = False):
    """ Record responses of all services for the benchmark tests

    Parameters
    ----------
    json_path : str
        Path to JSON file to save data to
    overwrite_json : bool, optional
        Whether to overwrite the JSON file, by default False

    """

    # Remove JSON path
    if overwrite_json:
        if os.path.exists(json_path):
            os.remove(json_path)

    # Get all test methods
    test_services = TestServices()
    test_services.json_path = json_path
    methods = dir(TestServices)
    test_methods = [m for m in methods if m.startswith("test_")]

    # Run tests, recording responses; failures are printed so they can be added
    # to "expected_failures" in the benchmark file.
    global RECORD_RESPONSES
    RECORD_RESPONSES = True
    bar = tqdm(test_methods, desc="Running tests")
    for m in bar:
        method = getattr(test_services, m)
        bar.set_description(f"Running {m}")
        for data in benchmark.values():
            try:
                method(data)
            except (
                Exception,
                pytest.fail.Exception,
                pytest.skip.Exception,
                pytest.xfail.Exception,
            ) as e:
                print(f"{m} failed for {data['iupac_name']}: {e}")

def test_multi_name_exhaustive_benchmark_case(monkeypatch):
    mr = MoleculeResolver()
    mr.supported_modes_by_services = {"svc1": ["name"], "svc2": ["name"]}
    mr._available_services = ["svc1", "svc2"]

    monkeypatch.setattr(mr, "_check_parameters", lambda **_: None)
    monkeypatch.setattr(
        mr,
        "_check_and_flatten_identifiers_and_modes",
        lambda identifiers, modes: (
            ["isopropanol", "2-propanol"],
            ["name", "name"],
            ["isopropanol", "2-propanol"],
            set(),
            None,
        ),
    )
    monkeypatch.setattr(mr, "standardize_SMILES", lambda smiles: smiles)

    def _adapter_result(smiles, synonym, identifier):
        return SimpleNamespace(
            molecule=SimpleNamespace(SMILES=smiles),
            synonyms=[synonym],
            cas={"67-63-0"},
            additional_information="mock",
            mode_used="name",
            identifier_used=identifier,
        )

    result_map = {
        ("svc1", "isopropanol"): _adapter_result("CC(O)C", "isopropanol", "isopropanol"),
        ("svc1", "2-propanol"): _adapter_result("CC(C)O", "2-propanol", "2-propanol"),
        ("svc2", "isopropanol"): _adapter_result("CC(C)O", "isopropanol", "isopropanol"),
    }

    monkeypatch.setattr(
        mr,
        "_resolve_identifier_with_adapter",
        lambda service, identifier, mode, *_: result_map.get((service, identifier)),
    )

    result = mr.find_single_molecule(
        identifiers=["isopropanol", "2-propanol"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        search_strategy="exhaustive",
    )

    assert result is not None
    assert result.SMILES == "CC(C)O"


def test_consensus_vs_legacy_mode_comparison(monkeypatch):
    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    legacy_winner = Molecule(
        "CCO", ["ethanol"], [], "svc1", "name", "svc1", 1, "ethanol"
    )
    consensus_winner = Molecule(
        "CCO", ["ethyl alcohol"], [], "svc2", "name", "svc2", 1, "ethyl alcohol"
    )

    def fake_find_single_molecule(*args, **kwargs):
        service = kwargs["services_to_use"][0]
        if service == "svc1":
            return legacy_winner
        return consensus_winner

    monkeypatch.setattr(mr, "find_single_molecule", fake_find_single_molecule)
    monkeypatch.setattr(
        mr,
        "filter_molecules",
        lambda molecules, *_: [molecule for molecule in molecules if molecule is not None],
    )
    monkeypatch.setattr(
        mr,
        "group_molecules_by_structure",
        lambda molecules, *_: {"CCO": molecules},
    )

    legacy_result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol"],
        modes=["name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="legacy",
    )
    consensus_result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol"],
        modes=["name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="consensus",
    )

    assert legacy_result is not None
    assert consensus_result is not None
    assert legacy_result.SMILES == consensus_result.SMILES


def test_strict_isomer_acceptance_case(monkeypatch):
    mr = MoleculeResolver()
    mr._available_services = ["svc1", "svc2"]

    mol_a = Molecule("CCO", ["ethanol"], [], "svc1", "name", "svc1", 1, "ethanol")
    mol_b = Molecule(
        "C[C@H](O)C", ["(S)-2-butanol"], [], "svc2", "name", "svc2", 1, "(S)-2-butanol"
    )

    def fake_find_single_molecule(*args, **kwargs):
        service = kwargs["services_to_use"][0]
        return mol_a if service == "svc1" else mol_b

    monkeypatch.setattr(mr, "find_single_molecule", fake_find_single_molecule)
    monkeypatch.setattr(
        mr,
        "filter_molecules",
        lambda molecules, *_: [molecule for molecule in molecules if molecule is not None],
    )
    monkeypatch.setattr(
        mr,
        "group_molecules_by_structure",
        lambda molecules, *_: {"CCO": [mol_a], "C[C@H](O)C": [mol_b]},
    )
    monkeypatch.setattr(
        mr,
        "_collect_opsin_isomer_matches",
        lambda grouped, smiles: {smiles[0]: False, smiles[1]: True},
    )

    result = mr.find_single_molecule_crosschecked(
        identifiers=["ethanol", "(S)-2-butanol"],
        modes=["name", "name"],
        services_to_use=["svc1", "svc2"],
        resolution_mode="strict_isomer",
    )

    assert result is not None
    assert result.SMILES == "C[C@H](O)C"


if __name__ == "__main__":
    generate_data(RESPONSES_PATH, overwrite_json=True)
