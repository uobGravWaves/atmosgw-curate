"""
Internal helpers shared by the limb and nadir loaders.

Purpose
-------
The parts of loading that are the same for every instrument family, so that
``load_limb_sounders`` and ``load_nadir_sounders`` treat times, options and
the record of a load the same way:

- times: ``to_datenums``, ``to_datetime64`` (user times to MATLAB datenums),
  ``TIME_HANDLING`` and ``expand_time_range`` (the MATLAB ``TimeHandling``
  rules), ``year_and_day``, ``day_label``;
- checks of shared options: ``check_flags``, ``check_time_handling``,
  ``check_additional_vars``, ``find_instrument``, ``warn_if_no_data``;
- the load record: ``record``, ``note_sources``, ``note`` (what a load read,
  did not find, and how each variable was obtained), ``choose_file`` (one
  file from those matching a day, recording when there was more than one);
- the physical range check of temperatures: ``outside_physical_range``;
- the output: ``COORD_ATTRS`` (attributes of latitude, longitude, time and
  altitude) and ``load_attrs`` (the dataset attributes recording the load),
  with ``options_json``.
- ``warn_multiple_matches``: one warning per load when a day had several files.
- ``record_attrs``: a variable's ``source_variable`` and ``comment`` from the record.

Data dicts use the MATLAB variable names (``Lat``, ``Temp`` ...) and datenum
times; see ``limb.py`` and ``nadir.py``.

Science
-------
No scientific processing, apart from the temperature limits of
``outside_physical_range`` (100 K and 400 K below 100 km, from the MATLAB
loaders).

Origin
------
Moved here unchanged from ``obs/_limb_tools.py``, ``obs/limb.py`` and
``obs/nadir.py`` (nadir design, checkpoint 6). The time rules, file choice
and range check are ports of the MATLAB ``get_limbsounders`` and
``prep_airs_3d`` (uobGravWaves/MultiInstrumentFunctions, corwin365/
MatlabFunctions), the group's own code; the load record is new.

Validation
----------
Through the loaders' tests (``tests/test_limb_*.py``, ``tests/test_nadir_*.py``);
the move left both loaders' outputs bit-identical.

AI involvement
--------------
AI-moved (Claude); to be reviewed by C. Wright.

Change log
----------
2026/11/25 [Claude] shared helpers moved here from _limb_tools.py, limb.py, nadir.py
2026/11/25 [Claude] warn_multiple_matches (now also for nadir loads)
2026/11/25 [Claude] record_attrs (from limb; now also used by nadir)
"""

import json
import math
import os
import re
import warnings

import numpy as np

from atmosgw.curate import __version__
from atmosgw.curate.matlab_time import (
    datenum_to_datetime,
    datetime_to_datenum,
    make_datetime,
)

TIME_HANDLING = ("strict", "generous", "fuzzy")

# attributes of the coordinates common to all formats
COORD_ATTRS = {
    "latitude": {"units": "degrees_north", "standard_name": "latitude"},
    "longitude": {"units": "degrees_east", "standard_name": "longitude"},
    "time": {"standard_name": "time"},
    "altitude": {"units": "km", "standard_name": "altitude"},
}


# times


def to_datenums(time_range):
    """Times as MATLAB datenums (float64 array of 1 or 2 values)."""
    values = np.atleast_1d(np.asarray(time_range, dtype=object)).ravel()
    if any(isinstance(v, (int, float, np.number)) and not isinstance(v, bool)
           for v in values):
        raise TypeError(
            "time_range must be dates or times (datetime64, ISO strings or "
            "datetime), not numbers. Convert MATLAB datenums with "
            "atmosgw.curate.matlab_time.datenum_to_datetime."
        )
    times = np.array([to_datetime64(v) for v in values])
    return np.atleast_1d(datetime_to_datenum(times))


def to_datetime64(value):
    """One time as datetime64[us]; a clear error for strings that are not real
    dates, e.g. MATLAB-style carry-over such as '2022-01-34'."""
    try:
        return np.datetime64(value, "us")
    except ValueError as error:
        message = f"time_range value {value!r} is not a valid date or time ({error})."
        parts = re.fullmatch(r"\s*(\d{4})-(\d{1,2})-(\d{1,3})\s*", str(value))
        if parts:
            year, month, day = (int(p) for p in parts.groups())
            if month >= 1 and day >= 1:
                meant = str(make_datetime(year, month, day))[:10]
                message += (
                    f" Days or months beyond the end of the month do not carry "
                    f"over as in MATLAB's datenum: for that, use {meant!r}, or "
                    f"atmosgw.curate.matlab_time.make_datetime({year}, {month}, "
                    f"{day})."
                )
        raise ValueError(message) from error


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


def year_and_day(day_number):
    """Calendar year and day of year (1-based) of an integer datenum."""
    t = datenum_to_datetime(day_number)
    year = int(t.astype("datetime64[Y]").astype(int)) + 1970
    day = int((t.astype("datetime64[D]") - t.astype("datetime64[Y]")).astype(int)) + 1
    return year, day


def day_label(day_number):
    """ISO date (yyyy-mm-dd) of a datenum."""
    return str(datenum_to_datetime(math.floor(day_number)))[:10]


# checks of shared options


def check_flags(options, names):
    """TypeError unless each named option is True or False."""
    for flag in names:
        if not isinstance(options[flag], (bool, np.bool_)):
            raise TypeError(f"{flag} must be True or False.")


def check_time_handling(options):
    """ValueError unless options["time_handling"] is a known mode."""
    if options["time_handling"] not in TIME_HANDLING:
        raise ValueError(f"time_handling must be one of {TIME_HANDLING}.")


def check_additional_vars(options):
    """TypeError unless options["additional_vars"] is a list of names."""
    additional = options["additional_vars"]
    if isinstance(additional, str) or not all(isinstance(a, str) for a in additional):
        raise TypeError("additional_vars must be a list of variable names.")


def find_instrument(instrument, instruments, note=""):
    """The registry name matching instrument (any case); ValueError listing the
    choices (and note) otherwise."""
    names = {name.lower(): name for name in instruments}
    if not isinstance(instrument, str) or instrument.lower() not in names:
        raise ValueError(f"Unknown instrument {instrument!r}. Choose from: "
                         f"{', '.join(instruments)}{note}.")
    return names[instrument.lower()]


def warn_if_no_data(raw, name, available, stacklevel):
    """Warn if datenums raw lie wholly outside the (first, last) datenums an
    instrument has data for."""
    if raw.min() > available[1] or raw.max() < available[0]:
        first, last = (str(datenum_to_datetime(t, "ms"))[:10] for t in available)
        warnings.warn(f"Data for {name} is only available from {first} to {last}",
                      stacklevel=stacklevel + 1)


# the load record


def record(inst_info):
    """
    The record of a load, kept in inst_info["record"] (created if absent):

    - ``files_not_found``: what was looked for and not found (dates for daily
      files, e.g. "2007-01-03"; month folders/files such as "January2007" for
      SABER; years for ALOMAR; granules for nadir sounders);
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


# range check


def outside_physical_range(temp, alt=None):
    """
    True where a temperature is physically implausible: below 100 K, or above
    400 K (only below 100 km if heights alt are given, as the MATLAB limb
    loader). NaN is not flagged. alt broadcasts against temp.
    """
    temp = np.asarray(temp)
    with np.errstate(invalid="ignore"):
        hot = temp > 400
        if alt is not None:
            hot = hot & (np.asarray(alt) < 100)
        return (temp < 100) | hot


# the output


def options_json(call):
    """The arguments of a loader call as a JSON string (arrays as lists,
    non-finite numbers and other objects as text)."""
    def plain(value):
        if isinstance(value, dict):
            return {k: plain(v) for k, v in value.items()}
        if isinstance(value, (list, tuple, np.ndarray)):
            return [plain(v) for v in np.asarray(value, dtype=object).ravel()]
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, (int, np.integer)):
            return int(value)
        if isinstance(value, (float, np.floating)):
            return float(value) if np.isfinite(value) else str(value)
        if value is None or isinstance(value, str):
            return value
        return str(value)

    return json.dumps(plain(call))


def load_attrs(settings, call, source):
    """
    The dataset attributes recording a load, in output order: time_range (as
    expanded), source (the public loader), atmosgw_curate_version, options
    (the call's arguments as JSON; "" without call), files_read,
    files_not_found and multiple_matches (from settings["record"]).
    """
    first, last = (str(datenum_to_datetime(t, "ms")) for t in settings["time_range"])
    return {
        "time_range": f"{first} to {last}",
        "source": source,
        "atmosgw_curate_version": __version__,
        "options": options_json(call) if call is not None else "",
        "files_read": ", ".join(settings["files_read"]),
        "files_not_found": ", ".join(settings["record"]["files_not_found"]),
        "multiple_matches": "; ".join(settings["record"]["multiple_matches"]),
    }


def warn_multiple_matches(settings, stacklevel):
    """Warn once per load if any day had more than one matching file (all are
    listed in settings["record"]["multiple_matches"] and the output attribute)."""
    multiple = settings["record"]["multiple_matches"]
    if multiple:
        warnings.warn(
            f"{len(multiple)} day(s) had more than one matching "
            f"{settings['instrument']} file; one was used for each (e.g. "
            f"{multiple[0]}). All are listed in "
            "the output attribute 'multiple_matches'.",
            stacklevel=stacklevel + 1,
        )


def record_attrs(name, record, interpolated=False):
    """
    ``source_variable`` and ``comment`` of one output variable (MATLAB name
    name) from the load record. Heights on a new grid (name "Alt" with
    interpolated True) have no source variable; their comment names the
    file heights the data were interpolated from.
    """
    attrs = {}
    if name in record["source_variables"]:
        attrs["source_variable"] = record["source_variables"][name]
    if name in record["notes"]:
        attrs["comment"] = record["notes"][name]
    if name == "Alt" and interpolated:
        attrs.pop("source_variable", None)
        origin = record["notes"].get(name) or record["source_variables"].get(name)
        attrs["comment"] = "the height scale; the data were interpolated to it" + (
            f" from the file's heights ({origin})" if origin else "")
    return attrs
