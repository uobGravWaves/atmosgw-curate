"""
Where this computer keeps its data.

Purpose
-------
Find the root folder under which instrument data live (e.g. ``<root>/ACE/raw``),
so that loaders work on any machine without hard-coded paths. Instrument
sub-folders are defined by the loaders, not here.

Inputs and outputs
------------------
``local_data_dir(path=None, must_exist=True)``

=============  ===========  =======  ================================================
Option         Type         Default  Description
=============  ===========  =======  ================================================
``path``       str or Path  None     Use this folder; skips the lookup below
``must_exist`` bool         True     Raise an error if the folder does not exist
=============  ===========  =======  ================================================

Returns a ``pathlib.Path`` (``~`` and environment variables expanded).

The folder is taken from the first of these that is set:

1. the ``path`` argument;
2. the environment variable ``ATMOSGW_DATA``;
3. ``data_dir`` in the user's settings file, ``~/.config/atmosgw/config.toml``
   (e.g. ``data_dir = "C:/Data"``), or in the file named by the environment
   variable ``ATMOSGW_CONFIG``.

If none is set, a ``DataDirNotSetError`` explains how to set one. Each machine
needs this set once; the settings file lives in the home folder, so it is not
shared between machines.

Science
-------
No scientific processing.

Origin
------
Replaces ``cjwtools.paths.local_data_dir`` (python-common, C. Wright), which
chose the folder from the computer's name and so worked only on the author's
machines. Lookup design agreed with C. Wright, 2026-11-23.

Validation
----------
``tests/test_paths.py``: each source on its own, their order of precedence,
``~`` expansion, and the errors for nothing set, a missing folder and a
settings file without ``data_dir``.

AI involvement
--------------
AI-originated (Claude), to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.paths import local_data_dir
>>> ace_dir = local_data_dir() / "ACE" / "raw"

Status
------
Preliminary.

Change log
----------
2026/11/23 [Claude] initial version

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import os
import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # Python 3.10
    import tomli as tomllib

DATA_DIR_ENV = "ATMOSGW_DATA"
CONFIG_ENV = "ATMOSGW_CONFIG"


class DataDirNotSetError(RuntimeError):
    """No data folder has been set for this computer."""


def local_data_dir(path=None, must_exist=True):
    """Return this computer's data root folder. See the module docstring."""
    if path is not None:
        folder, source = path, "the path argument"
    elif os.environ.get(DATA_DIR_ENV):
        folder, source = os.environ[DATA_DIR_ENV], f"${DATA_DIR_ENV}"
    else:
        config_file = config_path()
        folder = _data_dir_from_config(config_file)
        source = str(config_file)
        if folder is None:
            raise DataDirNotSetError(
                "No data folder set for this computer. Either set the "
                f"environment variable {DATA_DIR_ENV}, or create {config_file} "
                'containing a line like:  data_dir = "C:/Data"'
            )

    folder = Path(os.path.expandvars(str(folder))).expanduser()
    if must_exist and not folder.is_dir():
        raise FileNotFoundError(f"Data folder {folder} (from {source}) does not exist.")
    return folder


def config_path():
    """Path of the user's settings file (it may not exist)."""
    if os.environ.get(CONFIG_ENV):
        return Path(os.environ[CONFIG_ENV]).expanduser()
    return Path.home() / ".config" / "atmosgw" / "config.toml"


# supporting functions


def _data_dir_from_config(config_file):
    """``data_dir`` from the settings file, or None if the file or key is absent."""
    if not config_file.is_file():
        return None
    with open(config_file, "rb") as f:
        settings = tomllib.load(f)
    return settings.get("data_dir")
