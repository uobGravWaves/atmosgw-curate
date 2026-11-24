"""
Where this computer keeps its data.

Purpose
-------
Find the folders that hold instrument data and test reference files, so that
code works on any machine without hard-coded paths.

Inputs and outputs
------------------
``local_data_dir(path=None, must_exist=True)``: the data root folder.

``dataset_dir(name, kind="raw", must_exist=True)``: the folder of one dataset.

``reference_dir(must_exist=True)``: the folder of reference outputs used by
the tests (e.g. MATLAB results that ports are compared against).

=============  ===========  =======  ================================================
Option         Type         Default  Description
=============  ===========  =======  ================================================
``path``       str or Path  None     Use this folder; skips the lookup below
``name``       str          -        Dataset name, e.g. ``"SABER"`` (exact case)
``kind``       str          "raw"    Which copy of the dataset, e.g. ``"raw"``,
                                     ``"pws"`` (see ``datasets.toml``)
``must_exist`` bool         True     Raise an error if the folder does not exist
=============  ===========  =======  ================================================

Each returns a ``pathlib.Path`` (``~`` and environment variables expanded).

**Data root.** Taken from the first of these that is set:

1. the ``path`` argument;
2. the environment variable ``ATMOSGW_DATA``;
3. ``data_dir`` in the user's settings file, ``~/.config/atmosgw/config.toml``
   (e.g. ``data_dir = "C:/Data"``), or in the file named by the environment
   variable ``ATMOSGW_CONFIG``.

If none is set, a ``DataDirNotSetError`` explains how to set one. Each machine
needs this set once; the settings file lives in the home folder, so it is not
shared between machines.

**Dataset folders.** The package's ``datasets.toml`` gives every dataset's
default location, relative to the data root. Any entry can be overridden in
the settings file under the same heading::

    data_dir = "D:/Data"

    [datasets.SABER]
    raw = "E:/big_disk/saber"     # absolute: used as it is
    pws = "SABER_v2/pws"          # relative: under data_dir

A settings file may also add datasets that ``datasets.toml`` does not list.
An unknown dataset or kind raises ``KeyError``.

**Reference folder.** ``reference_dir`` in the settings file, e.g.
``reference_dir = "tests/references"`` (relative: under data_dir) or an absolute
path. If it is not set, ``DataDirNotSetError`` is raised; tests that need it
catch this and skip.

Science
-------
No scientific processing.

Origin
------
``local_data_dir`` replaces ``cjwtools.paths.local_data_dir`` (python-common,
C. Wright), which chose the folder from the computer's name and so worked only
on the author's machines. Lookup design agreed with C. Wright, 2026-11-23.
``dataset_dir`` and ``reference_dir`` replace the folder names written into
the group's MATLAB loaders (e.g. ``[LocalDataDir,'/SABER/raw/']`` in
``get_limbsounders``); design agreed with C. Wright, 2026-11-24.

Validation
----------
``tests/test_paths.py``: each source on its own, their order of precedence,
``~`` expansion, and the errors for nothing set, a missing folder and a
settings file without ``data_dir``; dataset defaults, absolute and relative
overrides, datasets added by the user, unknown names and kinds; the reference
folder set, relative and unset.

AI involvement
--------------
AI-originated (Claude), to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.paths import dataset_dir, local_data_dir
>>> saber_dir = dataset_dir("SABER")            # e.g. <data root>/SABER/raw
>>> saber_pws = dataset_dir("SABER", "pws")
>>> root = local_data_dir()

Status
------
Preliminary.

Change log
----------
2026/11/23 [Claude] initial version
2026/11/24 [Claude] added dataset_dir() with datasets.toml, and reference_dir()
2026/11/24 [Claude] tomllib only (Python >= 3.11; tomli fallback removed)

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import os
import tomllib
from importlib.resources import files
from pathlib import Path

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


def dataset_dir(name, kind="raw", must_exist=True):
    """Return the folder of one dataset. See the module docstring."""
    defaults = _read_toml(files("atmosgw.curate") / "datasets.toml")["datasets"]
    config_file = config_path()
    user = _read_settings(config_file).get("datasets", {})

    if name not in defaults and name not in user:
        known = ", ".join(sorted(set(defaults) | set(user)))
        raise KeyError(f"Unknown dataset {name!r}. Known datasets: {known}.")

    entry = {**defaults.get(name, {}), **user.get(name, {})}
    if kind not in entry:
        raise KeyError(
            f"Dataset {name!r} has no {kind!r} folder (it has: "
            f"{', '.join(sorted(entry))}). Add one under [datasets.{name}] "
            f"in {config_file}."
        )
    source = (
        str(config_file) if kind in user.get(name, {}) else "the package defaults"
    )

    folder = Path(os.path.expandvars(str(entry[kind]))).expanduser()
    if not folder.is_absolute():
        folder = local_data_dir(must_exist=must_exist) / folder
    if must_exist and not folder.is_dir():
        raise FileNotFoundError(
            f"Folder {folder} for dataset {name!r} ({kind}, from {source}) "
            "does not exist."
        )
    return folder


def reference_dir(must_exist=True):
    """Return the folder of test reference outputs. See the module docstring."""
    config_file = config_path()
    entry = _read_settings(config_file).get("reference_dir")
    if entry is None:
        raise DataDirNotSetError(
            f"No reference folder set. Add a line like  reference_dir = "
            f'"tests/references"  to {config_file} (relative paths are under '
            "data_dir)."
        )
    folder = Path(os.path.expandvars(str(entry))).expanduser()
    if not folder.is_absolute():
        folder = local_data_dir(must_exist=must_exist) / folder
    if must_exist and not folder.is_dir():
        raise FileNotFoundError(
            f"Reference folder {folder} (from {config_file}) does not exist."
        )
    return folder


def config_path():
    """Path of the user's settings file (it may not exist)."""
    if os.environ.get(CONFIG_ENV):
        return Path(os.environ[CONFIG_ENV]).expanduser()
    return Path.home() / ".config" / "atmosgw" / "config.toml"


# supporting functions


def _data_dir_from_config(config_file):
    """``data_dir`` from the settings file, or None if the file or key is absent."""
    return _read_settings(config_file).get("data_dir")


def _read_settings(config_file):
    """Contents of the settings file, or an empty dict if there is none."""
    if not config_file.is_file():
        return {}
    return _read_toml(config_file)


def _read_toml(file):
    """Read a TOML file (a Path or an importlib.resources file)."""
    with file.open("rb") as f:
        return tomllib.load(f)
