"""
Finding and reading data files.

Purpose
-------
Small, format-level helpers shared by the instrument loaders: find files by a
wildcard pattern anywhere below a folder, and read MATLAB ``.mat`` files.

Inputs and outputs
------------------
``find_files(root, pattern, case_sensitive=False)``: list of paths (str) of
files, and folders, below ``root`` whose path matches ``pattern``.

``read_mat(path)``: dict of the variables in a ``.mat`` file.

``read_hdf5(path, name)``: one HDF5 dataset as a float64 array (as MATLAB's
``H5D.read(..., 'H5T_NATIVE_DOUBLE', ...)``), in the file's (C) dimension
order, i.e. the transpose of what MATLAB returns for a 2-D dataset.

``read_netcdf(path, names, rows=None)``: dict of netCDF variables, processed
as the group's MATLAB ``nph_getnet`` (used through ``rCDF``) does: numeric
variables converted to float64; values equal to ``_FillValue`` and then to
``missing_value`` set to NaN; then multiplied by ``scale_factor`` and
``add_offset`` added, if those attributes exist. (netCDF4-python's own
masking and scaling is switched off: it also applies ``valid_range`` and
can return float32.) Text variables are returned as read. Arrays are in the
file's (C) dimension order: for 2-D variables this is the transpose of what
``rCDF`` returns (``rCDF`` reverses dimensions only for 3-D and higher). If
``rows`` is given (an index array), only those indices of the first dimension
are read, which saves memory for large files.

``find_files`` reproduces the group's MATLAB ``wildcardsearch`` (which calls
``regexpdir``):

- ``pattern`` is a wildcard pattern: ``*`` matches any characters, ``?`` any
  one character, ``.`` only a dot; other characters are used as regular
  expression characters, as in MATLAB's ``regexptranslate('wildcard', ...)``.
  The path separator character (``:`` on Linux, ``;`` on Windows) separates
  alternative patterns.
- It matches anywhere in the path relative to ``root`` (not only the file
  name), and ignores case unless ``case_sensitive``.
- Search order: at each level, sub-folders in name order (each one listed if
  it matches, then searched), then the files in name order. Loaders that take
  the first match therefore get the same file as MATLAB.
- Symbolic links to folders are followed.

Difference from MATLAB: ``regexpdir`` remembers the root folder of its first
call in a MATLAB session, so later calls match against the full path rather
than the path relative to their own root. ``find_files`` always uses the
relative path. Results differ only if the root path itself matches the
pattern (Q5 in ``docs/quirks.md``).

``read_mat`` reads every MATLAB file version (including v7.3) with ``mat-io``.
Arrays keep MATLAB's shapes (a row vector is 1 x n, a scalar 1 x 1). MATLAB
tables become ``pandas.DataFrame``. A 1 x 1 MATLAB struct becomes a dict of
its fields (recursively); struct arrays and cell arrays are returned as
``mat-io`` gives them.

Science
-------
No scientific processing.

Origin
------
``find_files``: new code written to reproduce the behaviour of
``wildcardsearch.m`` and ``regexpdir.m`` (B. C. Hamans, 2007; third-party
functions used by the group's MATLAB code, in corwin365/MatlabFunctions,
commit c943399). No code or text is copied from them. ``read_mat``: new.
``read_hdf5``: new. ``read_netcdf``: port of the processing in
``nph_getnet.m`` (N. Hindley, corwin365/MatlabFunctions, commit c943399), the
group's own code.

Validation
----------
``tests/test_readers.py``: wildcard translation, case, matching on the
relative path, alternatives, search order, folders that match; ``read_mat``
shapes, structs and v5 files. ``find_files`` is also exercised by every
loader test that reads files found by a pattern.

AI involvement
--------------
AI-originated (Claude), to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.readers import find_files, read_mat  # doctest: +SKIP
>>> files = find_files("/data/ACE/raw", "*ace_v52_2007d001*.mat")  # doctest: +SKIP
>>> data = read_mat(files[0])  # doctest: +SKIP

Status
------
Preliminary.

Change log
----------
2026/11/24 [Claude] initial version

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import os
import re

import h5py
import matio
import netCDF4
import numpy as np


def find_files(root, pattern, case_sensitive=False):
    """Files (and folders) below root whose relative path matches a wildcard
    pattern, in MATLAB wildcardsearch order. See the module docstring."""
    alternatives = pattern.split(os.pathsep)
    regex = "|".join(f"({_wildcard_to_regex(a)})" for a in alternatives)
    flags = 0 if case_sensitive else re.IGNORECASE
    compiled = re.compile(regex, flags)
    root = str(root)
    if not root.endswith(os.sep):
        root = root + os.sep
    return _search(root, root, compiled)


def read_mat(path):
    """Variables in a MATLAB .mat file, as a dict. See the module docstring."""
    contents = matio.load_from_mat(str(path))
    return {
        name: _struct_to_dict(value)
        for name, value in contents.items()
        if not name.startswith("__")
    }


def read_hdf5(path, name):
    """One HDF5 dataset as float64, in file (C) order. See the module docstring."""
    with h5py.File(path, "r") as f:
        return np.asarray(f[name][()], dtype=np.float64)


def read_netcdf(path, names, rows=None):
    """netCDF variables processed as nph_getnet. See the module docstring."""
    out = {}
    with netCDF4.Dataset(path) as f:
        f.set_auto_maskandscale(False)
        f.set_auto_chartostring(False)
        for name in names:
            var = f.variables[name]
            raw = var[rows, ...] if rows is not None else var[...]
            raw = np.asarray(raw)
            if raw.dtype.kind not in "biuf":
                out[name] = raw
                continue
            value = raw.astype(np.float64)
            attrs = {a.lower(): var.getncattr(a) for a in var.ncattrs()}
            for key in ("_fillvalue", "missing_value"):
                flag = np.asarray(attrs.get(key, ""))
                if flag.dtype.kind in "biuf":
                    value[value == np.float64(flag.ravel()[0])] = np.nan
            if "scale_factor" in attrs:
                value = value * np.float64(attrs["scale_factor"])
            if "add_offset" in attrs:
                value = value + np.float64(attrs["add_offset"])
            out[name] = value
    return out


# supporting functions


def _wildcard_to_regex(pattern):
    """MATLAB regexptranslate('wildcard', ...): only '.', '*' and '?' change."""
    return pattern.replace(".", r"\.").replace("*", ".*").replace("?", ".")


def _search(folder, root, regex):
    """Depth-first search in regexpdir order (sub-folders, then files)."""
    found = []
    try:
        names = sorted(os.listdir(folder))
    except FileNotFoundError:
        return found
    folders = [n for n in names if os.path.isdir(os.path.join(folder, n))]
    files = [n for n in names if not os.path.isdir(os.path.join(folder, n))]
    for name in folders:
        path = os.path.join(folder, name) + os.sep
        if regex.search(path[len(root):]):
            found.append(path)
        found.extend(_search(path, root, regex))
    for name in files:
        path = os.path.join(folder, name)
        if regex.search(path[len(root):]):
            found.append(path)
    return found


def _struct_to_dict(value):
    """A 1 x 1 MATLAB struct (structured array) as a dict of its fields."""
    if isinstance(value, np.ndarray) and value.dtype.names and value.size == 1:
        return {
            name: _struct_to_dict(value[name].reshape(-1)[0])
            for name in value.dtype.names
        }
    return value
