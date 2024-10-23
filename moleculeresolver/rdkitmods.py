from contextlib import ContextDecorator
from rdkit import rdBase


class disabling_rdkit_logger(ContextDecorator):
    """
    A context manager for disabling RDKit logging
    based on https://github.com/rdkit/rdkit/issues/2320#issuecomment-731261149

    This class allows for temporary suppression of RDKit log messages at various levels
    (error, warning, info, debug) to reduce noise in the output during specific operations.

    Attributes:
        previous_status (dict): The log status before entering the context manager.
        desired_status (dict): The log status desired during the context manager's execution.
    """

    def __init__(
        self,
        mute_errors: bool = True,
        mute_warning: bool = True,
        mute_info: bool = True,
        mute_debug: bool = True,
    ) -> None:
        """
        Initializes the disabling_rdkit_logger context manager.

        Args:
            mute_errors (bool): If True, suppress error messages. Defaults to True.
            mute_warning (bool): If True, suppress warning messages. Defaults to True.
            mute_info (bool): If True, suppress info messages. Defaults to True.
            mute_debug (bool): If True, suppress debug messages. Defaults to True.
        """

        self.previous_status = self._get_log_status()

        self.desired_status = {
            "rdApp.error": not mute_errors,
            "rdApp.warning": not mute_warning,
            "rdApp.debug": not mute_debug,
            "rdApp.info": not mute_info,
        }

    def _get_log_status(self) -> dict[str, bool]:
        """
        Get the current log status of RDKit logs.

        Returns:
            dict[str, bool]: A dictionary indicating the log status (enabled/disabled) for each log level.
        """
        log_status = rdBase.LogStatus()
        log_status = {
            st.split(":")[0]: st.split(":")[1] for st in log_status.split("\n")
        }
        log_status = {k: v == "enabled" for k, v in log_status.items()}
        return log_status

    def _apply_log_status(self, log_status: dict[str, bool]) -> None:
        """
        Apply an RDKit log status.

        Args:
            log_status (dict[str, bool]): A dictionary with log levels as keys and their desired status (True/False) as values.
        """
        for k, v in log_status.items():
            if v:
                rdBase.EnableLog(k)
            else:
                rdBase.DisableLog(k)

    def __enter__(self) -> "disabling_rdkit_logger":
        """
        Enter the runtime context related to this object.

        Applies the desired log status when entering the context.

        Returns:
            disabling_rdkit_logger: The context manager itself.
        """
        self._apply_log_status(self.desired_status)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback) -> None:
        """
        Exit the runtime context related to this object.

        Restores the previous log status when exiting the context.
        """
        self._apply_log_status(self.previous_status)
