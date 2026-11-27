"""
Internal helpers for the limb-sounder loader (``load_limb_sounders``).

- ``conform_vars``: put user-requested extra variables on the [profiles x
  levels] grid (port of ``module_conform_vars.m``).
- ``loop_days``: the day loop shared by readers with one file per day.
- ``newest_mls_version``: of several MLS files for a day, the newest version
  (new, not in the MATLAB); ``newest_version`` and ``VERSION_IN_NAME``: the
  same for HIRDLS, SABER, SOFIE and MIPAS (Q19).
- ``record_if_empty``: record a file with no finite main variable.
- ``ESTIMATED_PRES``, ``note_estimated_pres``: the note recorded for pressure
  computed from height.

Helpers shared with the nadir loader (time rules, the load record, file
choice) are in ``obs/_common.py``.

Ported from the MATLAB ``get_limbsounders`` and its modules
(uobGravWaves/MultiInstrumentFunctions, commit 78eef78), the group's own code.
Data dicts use the MATLAB variable names (``Lat``, ``Temp`` ...) and datenum
times; see ``limb.py``.

Change log:
2026/11/25 [Claude] time rules, load record, choose_file moved to _common.py
2026/11/25 [Claude] note_estimated_pres: files with estimated pressure recorded
2026/11/25 [Claude] note_estimated_pres: which points is now the readers' PresEstimated
2026/11/27 [Claude] newest_version (Q19); record_if_empty
"""

import math
import os
import re
import warnings

import numpy as np
import pandas as pd

from atmosgw.curate.obs._common import day_label, note, record
from atmosgw.curate.profiles import concat_profiles


def conform_vars(store, settings):
    """
    Put settings["additional_vars"] onto the [profiles x levels] grid of
    store["Lat"] (port of module_conform_vars.m):

    - [profiles x levels (x ...)]: unchanged;
    - [levels x profiles (x ...)]: first two dimensions swapped;
    - vector with one value per profile: repeated over levels;
    - vector with one value per level: repeated over profiles;
    - table with one row per profile: unchanged;
    - anything else: ValueError, unless settings["ignore_non_modal"], when it
      is passed through.

    If there are as many profiles as levels, a matching square field is
    accepted with a warning, and a vector is an error (ambiguous).
    """
    n_prof, n_lev = np.shape(store["Lat"])[:2]
    for name in settings["additional_vars"]:
        if name not in store:
            raise KeyError(
                f"AdditionalVar '{name}' was not found in the "
                f"{settings['instrument']} data"
            )
        field = store[name]
        if isinstance(field, pd.DataFrame):
            if len(field) == n_prof:
                continue
        elif isinstance(field, np.ndarray) and field.dtype.kind in "biuf":
            field = np.atleast_2d(field)
            size = field.shape[:2]
            is_vector = field.ndim == 2 and min(size) == 1
            if size == (n_prof, n_lev):
                if n_prof == n_lev and n_prof > 1:
                    warnings.warn(
                        f"AdditionalVar '{name}': file has as many profiles as levels "
                        f"({n_prof}), so its orientation cannot be checked; assumed "
                        "[profiles x levels]",
                        stacklevel=2,
                    )
                continue
            if size == (n_lev, n_prof) and n_prof != n_lev:
                store[name] = np.swapaxes(field, 0, 1)
                continue
            if is_vector and field.size == n_prof and n_prof != n_lev:
                store[name] = np.tile(field.reshape(-1, 1, order="F"), (1, n_lev))
                continue
            if is_vector and field.size == n_lev and n_prof != n_lev:
                store[name] = np.tile(field.reshape(1, -1, order="F"), (n_prof, 1))
                continue
        if settings["ignore_non_modal"]:
            continue
        raise ValueError(
            f"AdditionalVar '{name}' from {settings['instrument']} has size "
            f"{list(np.shape(field))} ({type(field).__name__}), which is not "
            f"per-point [{n_prof} x {n_lev}], per-profile ({n_prof}) or per-level "
            f"({n_lev}). Set ignore_non_modal=True to pass it through without "
            "processing."
        )
    return store


def loop_days(settings, var_names, find_file, read_file, label, inst_info=None):
    """
    Day loop for readers with one file per day.

    For each day (an integer datenum) from floor(start) to floor(end):
    find_file(day) gives the file's path or None (skipped, and the date added
    to ``record(inst_info)["files_not_found"]`` if inst_info is given);
    read_file(path, file_count, day) gives that day's variables, already
    conformed. Returns (data, file_list): the days
    stacked with ``concat_profiles``, and the file names (without folders) in
    the order read. file_count is the 1-based number of the file, as stored
    in SourceFile.
    """
    first, last = (math.floor(t) for t in settings["time_range"])
    days, file_list = [], []
    for k, day in enumerate(range(first, last + 1), start=1):
        if settings["verbose"]:
            print(f"----> Loading {label} data: day {k} of {last - first + 1}")
        path = find_file(day)
        if path is None or not os.path.isfile(path):
            if inst_info is not None:
                record(inst_info)["files_not_found"].append(day_label(day))
            continue
        file_list.append(os.path.basename(path))
        day_data = read_file(path, len(file_list), day)
        for main in ("Temp", "Temp_Residual"):
            if main in day_data and np.size(day_data[main]) and not np.any(
                    np.isfinite(np.asarray(day_data[main], dtype=float))):
                if inst_info is not None:
                    record(inst_info).setdefault("files_without_data", []).append(
                        os.path.basename(path))
                break
        days.append(day_data)
    return concat_profiles(days, var_names, consume=True), file_list


def record_if_empty(inst_info, path, main):
    """Add path to ``files_without_data`` if main (the file's main variable,
    as read for the request) has values but none finite."""
    if main is not None and np.size(main) and not np.any(
            np.isfinite(np.asarray(main, dtype=float))):
        record(inst_info).setdefault("files_without_data", []).append(
            os.path.basename(path))


_MLS_VERSION = re.compile(r"_v(\d+)-(\d+)-c(\d+)_", re.IGNORECASE)


def newest_mls_version(files):
    """
    The file with the newest MLS version, from names such as
    ``MLS-Aura_L2GP-Temperature_v06-04-c01_2026d001.he5``: highest version
    (06-04 beats 06-03 beats 05-04), then highest cycle (c02 beats c01).
    Names without a version count as oldest; of equal versions, the first.
    """
    def version(path):
        found = _MLS_VERSION.search(os.path.basename(path))
        return tuple(int(g) for g in found.groups()) if found else (-1, -1, -1)

    return max(files, key=version)  # (max returns the first of equal keys)


# data versions in file names (Q19), for newest_version
VERSION_IN_NAME = {
    "HIRDLS": re.compile(r"_v(\d+)-(\d+)-(\d+)-c(\d+)_"),  # _v07-00-20-c01_
    "SABER": re.compile(r"_v(\d+)\.(\d+)\.nc$"),  # _v2.0.nc
    "MIPAS": re.compile(r"_v(\d+)\.(\d+)\.nc$"),  # _v1.40.nc
    "SOFIE": re.compile(r"_(\d+)\.(\d+)\.nc$"),  # _01.3.nc
}


def newest_version(pattern):
    """
    A ``choose_file`` ``prefer``: of several files for a day, the one whose
    name has the newest data version (pattern: a compiled regex whose groups
    are the version's numbers, compared as integers, most significant first).
    Names without a version count as oldest; of equal versions, the first in
    search order (as MATLAB, which always takes the first: Q19).
    """
    def version(path):
        found = pattern.search(os.path.basename(path))
        return tuple(int(g) for g in found.groups()) if found else (-1,)

    return lambda files: max(files, key=version)  # (first of equal keys)


# note for pressure computed from height
ESTIMATED_PRES = ("estimated from height with atmosphere.height_to_pressure "
                  "(U.S. Standard Atmosphere 1976 [R1])")


def note_estimated_pres(inst_info):
    """Record that a file had no pressure, so it was estimated from height.
    Which points: the reader's ``PresEstimated`` (1 estimated, 0 read), asked
    for by the core when it loads Pres from ACE/GNSS/Misc or PW files."""
    note(inst_info, "Pres", ESTIMATED_PRES + ", for files without Pres")
