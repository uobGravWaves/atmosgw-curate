"""
Internal helpers for the limb-sounder loader (``load_limb_sounders``).

- ``expand_time_range``: the MATLAB ``TimeHandling`` rules, on datenums.
- ``conform_vars``: put user-requested extra variables on the [profiles x
  levels] grid (port of ``module_conform_vars.m``).
- ``loop_days``: the day loop shared by readers with one file per day.
- ``year_and_day``: calendar year and day of year of a datenum.

Ported from the MATLAB ``get_limbsounders`` and its modules
(uobGravWaves/MultiInstrumentFunctions, commit 78eef78), the group's own code.
Data dicts use the MATLAB variable names (``Lat``, ``Temp`` ...) and datenum
times; see ``limb.py``.
"""

import math
import os
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


def loop_days(settings, var_names, find_file, read_file, label):
    """
    Day loop for readers with one file per day.

    For each day (an integer datenum) from floor(start) to floor(end):
    find_file(day) gives the file's path or None (skipped);
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
