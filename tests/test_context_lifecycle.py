from __future__ import annotations

from dataclasses import dataclass

from moleculeresolver import MoleculeResolver


@dataclass
class _DummyContext:
    exited: bool = False

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.exited = True


@dataclass
class _DummyTempFolder:
    cleaned: bool = False

    def cleanup(self):
        self.cleaned = True


def test_enter_delegates_to_lifecycle_helpers(mocker):
    mr = MoleculeResolver()

    called = []
    mocker.patch.object(
        mr, "_enter_rdkit_log_context", side_effect=lambda: called.append("rdkit")
    )
    mocker.patch.object(
        mr, "_enter_molecule_cache_context", side_effect=lambda: called.append("cache")
    )
    mocker.patch.object(
        mr, "_create_opsin_tempfolder", side_effect=lambda: called.append("opsin")
    )

    result = mr.__enter__()

    assert result is mr
    assert called == ["rdkit", "cache", "opsin"]


def test_exit_delegates_to_teardown_helper(mocker):
    mr = MoleculeResolver()

    teardown_mock = mocker.patch.object(mr, "_teardown_runtime_contexts")

    mr.__exit__(None, None, None)
    teardown_mock.assert_called_once_with(error_ocurred=False)

    teardown_mock.reset_mock()
    mr.__exit__(RuntimeError, RuntimeError("boom"), object())
    teardown_mock.assert_called_once_with(error_ocurred=True)


def test_teardown_runtime_contexts_cleans_tempfolder_on_success():
    mr = MoleculeResolver()
    cache_ctx = _DummyContext()
    rdkit_ctx = _DummyContext()
    tempfolder = _DummyTempFolder()

    mr.molecule_cache = cache_ctx
    mr._disabling_rdkit_logger = rdkit_ctx
    mr._OPSIN_tempfolder = tempfolder

    mr._teardown_runtime_contexts(error_ocurred=False)

    assert cache_ctx.exited
    assert rdkit_ctx.exited
    assert tempfolder.cleaned


def test_teardown_runtime_contexts_keeps_tempfolder_on_error():
    mr = MoleculeResolver()
    cache_ctx = _DummyContext()
    rdkit_ctx = _DummyContext()
    tempfolder = _DummyTempFolder()

    mr.molecule_cache = cache_ctx
    mr._disabling_rdkit_logger = rdkit_ctx
    mr._OPSIN_tempfolder = tempfolder

    mr._teardown_runtime_contexts(error_ocurred=True)

    assert cache_ctx.exited
    assert rdkit_ctx.exited
    assert not tempfolder.cleaned
