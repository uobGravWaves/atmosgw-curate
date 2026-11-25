"""
Internal helpers for the limb-sounder loader (``load_limb_sounders``).

- ``expand_time_range``: the MATLAB ``TimeHandling`` rules, on datenums.
- ``conform_vars``: put user-requested extra variables on the [profiles x
  levels] grid (port of ``module_conform_vars.m``).
- ``loop_days``: the day loop shared by readers with one file per day.
- ``year_and_day``: calendar year and day of year of a datenum.
- ``choose_file``, ``newest_mls_version``: pick one file from those matching
  a day, and record when there was more than one (new, not in the MATLAB).
- ``record``, ``note_sources``, ``note``: what a load read, did not find, and
  how each variable was obtained, for the output attributes (new).

Ported from the MATLAB ``get_limbsounders`` and its modules
(uobGravWaves/MultiInstrumentFunctions, commit 78eef78), the group's own code.
Data dicts use the MATLAB variable names (``Lat``, ``Temp`` ...) and datenum
times; see ``limb.py``.
"""

import math
import os
import re
import warnings

import numpy as np
import pandas as pd

from atmosgw.curate.matlab_time import datenum_to_datetime
from atmosgw.curate.profiles import concat_profiles

TIME_HANDLING = ("strict", "generous", "fuzzy")


def expand_time_range(time_range, mode):
    """
    Apply the TimeHandling rules to one or two datenums.

    Returns ([start, end] as a float64 array, scalar_given). Modes, as in the
    MATLAB header (1, 2, 3 there):

    - "strict": exactly as given; one value t means the instant [t, t];
    - "generous": from the start of the first day to the end of the last day
      touched; an end exactly at midnight includes that whole day;
    - "fuzzy" (default): one non-integer value means its whole day; an end
      exactly at midnight is moved to just before the next midnight.

    "Just before" is 1e-8 day (about 1 ms) before midnight, as in MATLAB.
    """
    if mode not in TIME_HANDLING:
        raise ValueError(f"time_handling must be one of {TIME_HANDLING}, not {mode!r}.")
    tr = np.atleast_1d(np.asarray(time_range, dtype=float)).ravel()
    if tr.size > 2:
        raise ValueError(
            "time_range must be either one value (a day of interest) or two values "
            "(a time range)."
        )
    scalar = tr.size == 1
    if scalar:
        tr = np.array([tr[0], tr[0]])
    tr = np.sort(tr, kind="stable")
    if mode == "generous":
        tr[0] = math.floor(tr[0])
        if tr[1] % 1 == 0:
            tr[1] = tr[1] + 1 - 1e-8
        else:
            tr[1] = math.ceil(tr[1]) - 1e-8
    elif mode == "fuzzy":
        if scalar and tr[1] % 1 != 0:
            tr = np.array([math.floor(tr[0]), math.floor(tr[0]) + 1 - 1e-8])
        elif tr[1] % 1 == 0:
            tr[1] = tr[1] + 1 - 1e-8
    return tr, scalar


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
        days.append(read_file(path, len(file_list), day))
    return concat_profiles(days, var_names), file_list


def year_and_day(day_number):
    """Calendar year and day of year (1-based) of an integer datenum."""
    t = datenum_to_datetime(day_number)
    year = int(t.astype("datetime64[Y]").astype(int)) + 1970
    day = int((t.astype("datetime64[D]") - t.astype("datetime64[Y]")).astype(int)) + 1
    return year, day


def day_label(day_number):
    """ISO date (yyyy-mm-dd) of a datenum."""
    return str(datenum_to_datetime(math.floor(day_number)))[:10]


def record(inst_info):
    """
    The record of a load, kept in inst_info["record"] (created if absent):

    - ``files_not_found``: what was looked for and not found (dates for daily
      files, e.g. "2007-01-03"; month folders/files such as "January2007" for
      SABER; years for ALOMAR);
    - ``multiple_matches``: one line for each day with more than one matching
      file, naming them and the one used;
    - ``source_variables``: MATLAB name -> where it comes from in the file;
    - ``notes``: MATLAB name -> how it was obtained, where that is more than
      reading it (e.g. pressure estimated from height).
    """
    return inst_info.setdefault("record", {"files_not_found": [],
                                           "multiple_matches": [],
                                           "source_variables": {}, "notes": {}})


def note_sources(inst_info, sources):
    """Record where variables come from in the files (MATLAB name -> text)."""
    record(inst_info)["source_variables"].update(sources)


def note(inst_info, name, text):
    """Record how a variable was obtained (the first note for it is kept)."""
    record(inst_info)["notes"].setdefault(name, text)


def choose_file(found, inst_info, label, prefer=None):
    """
    One path from those found for a day (``readers.find_files`` output).

    Without ``prefer``: the first, as the MATLAB loaders do. With it:
    ``prefer(files)`` picks one of the files found (folders excluded). Either
    way, if more than one file matched, a line naming them and the one used
    is added to ``record(inst_info)["multiple_matches"]``. Returns the path,
    or None if nothing was found.
    """
    if not found:
        return None
    files = [path for path in found if os.path.isfile(path)]
    chosen = prefer(files) if prefer is not None and files else found[0]
    if len(files) > 1:
        names = ", ".join(os.path.basename(path) for path in files)
        record(inst_info)["multiple_matches"].append(
            f"{label}: {names} -> used {os.path.basename(chosen)}")
    return chosen


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


# note for pressure computed from height
ESTIMATED_PRES = ("estimated from height with atmosphere.height_to_pressure "
                  "(U.S. Standard Atmosphere 1976 [R1])")
