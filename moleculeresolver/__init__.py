from .moleculeresolver import Molecule
from .moleculeresolver import MoleculeResolver

from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("molecule-resolver")
except PackageNotFoundError:
    __version__ = "dev"