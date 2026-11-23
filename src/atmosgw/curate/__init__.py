"""AtmosGW Curate: loading, standardising and analysing atmospheric gravity-wave
observations.

Part of the AtmosGW family of packages (``atmosgw.curate``, later
``atmosgw.fuse``). Status: pre-alpha skeleton; no scientific code yet and not
validated against the MATLAB reference implementation.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("atmosgw-curate")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
