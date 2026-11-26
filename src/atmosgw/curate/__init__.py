"""AtmosGW Curate: atmospheric gravity-wave observations in common formats.

Loaders are in ``atmosgw.curate.obs`` (``load_limb_sounders``,
``load_nadir_sounders``, ``load_time_height``); the formats, data states and
``to_profiles`` are in ``atmosgw.curate.data``. Each module's docstring
describes its options, science and validation.

Part of the AtmosGW family of packages (``atmosgw.curate``, later
``atmosgw.fuse``). Status: pre-alpha (0.x). The ported loaders are checked
against the group's MATLAB code; what has been checked is listed in the
repository's ``README.md``.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("atmosgw-curate")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
