"""
Internal helpers shared by the limb and nadir loaders.

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
- ``note_version``, ``versions_text``, ``warn_mixed_versions``: the data
  version each file read states, the ``source_versions`` attribute, and one
  warning per load when the files read state different versions.
- ``record_attrs``: a variable's ``source_variable`` and ``comment`` from the record.

Data dicts use the MATLAB variable names (``Lat``, ``Temp`` ...) and datenum
times; see ``limb.py`` and ``nadir.py``.

Science
-------
The temperature limits of ``outside_physical_range``, from the MATLAB
loaders: below 100 K, or above 400 K (only below 100 km when heights are
passed, as by the limb and time-height loaders; the nadir loader passes
none, so at any height).

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

Status
------
AI-moved (Claude); to be reviewed by C. Wright.

Change log
----------
2026/11/25 [Claude] shared helpers moved here from _limb_tools.py, limb.py, nadir.py
2026/11/25 [Claude] warn_multiple_matches (moved from limb.py)
2026/11/25 [Claude] record_attrs (from limb; now also used by nadir)
2026/11/25 [Claude] outside_physical_range: floor_is_bad (IASI's 100 K lower bound)
2026/11/26 [Claude] nadir_fill_to_nan, note_fill_counts: AIRS/IASI fill, bad times
2026/11/26 [Claude] note_version, versions_text, warn_mixed_versions: source_versions
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
    if time_range is None or (np.ndim(time_range) and np.size(time_range) == 0):
        raise ValueError("time_range is empty: give one time (a day) or two "
                         "(start and end).")
    if np.asarray(time_range).dtype.kind == "M":  # (datetime64 of any unit)
        values = np.atleast_1d(np.asarray(time_range, dtype="datetime64[us]")).ravel()
    else:
        values = np.atleast_1d(np.asarray(time_range, dtype=object)).ravel()
    if any(isinstance(v, (int, float, np.number)) and not isinstance(v, bool)
           and not isinstance(v, np.datetime64) for v in values):
        raise TypeError(
            "time_range must be dates or times (datetime64, ISO strings or "
            "datetime), not numbers. Convert MATLAB datenums with "
            "atmosgw.curate.matlab_time.datenum_to_datetime."
        )
    times = np.array([to_datetime64(v) for v in values])
    if np.any(np.isnat(times)):
        raise ValueError("time_range contains NaT (not a time).")
    return np.atleast_1d(datetime_to_datenum(times))


def to_datetime64(value):
    """One time as datetime64[us]; a clear error for strings that are not real
    dates, e.g. MATLAB-style carry-over such as '2022-01-34'."""
    if isinstance(value, str) and re.fullmatch(r"\s*\d{8}\s*", value):
        v = value.strip()
        raise ValueError(f"time_range value {value!r} is read as the year {v}, not "
                         f"a date: write it as '{v[:4]}-{v[4:6]}-{v[6:]}'.")
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


def check_additional_vars(options, reserved=()):
    """TypeError unless options["additional_vars"] is a list of names;
    ValueError for a name in reserved (an output or working name, which the
    extra variable would replace). Returns the names, each once, in order."""
    additional = options["additional_vars"]
    if isinstance(additional, str) or not all(isinstance(a, str) for a in additional):
        raise TypeError("additional_vars must be a list of variable names.")
    clash = [a for a in additional if a in reserved or a.endswith("_flag")]
    if clash:
        raise ValueError(
            f"additional_vars {clash} would replace variables the loader makes "
            "itself (output names, flags, or the file variables those are made "
            "from); they cannot be loaded as extras.")
    return list(dict.fromkeys(additional))


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

    Folders whose names match are not files, so are left out (MATLAB took
    the first match even if it was a folder, then failed to read it).
    Without ``prefer``: the first file, as the MATLAB loaders do. With it:
    ``prefer(files)`` picks one. Either way, if more than one file matched, a
    line naming them and the one used is added to
    ``record(inst_info)["multiple_matches"]``. Returns the path, or None if
    no file was found.
    """
    files = [path for path in found if os.path.isfile(path)]
    if not files:
        return None
    chosen = prefer(files) if prefer is not None else files[0]
    if len(files) > 1:
        names = ", ".join(os.path.basename(path) for path in files)
        record(inst_info)["multiple_matches"].append(
            f"{label}: {names} -> used {os.path.basename(chosen)}")
    return chosen


# nadir fill values

# Fill values in the AIRS and IASI files, none declared as _FillValue (census of
# every AIRS daily and IASI file on eepc-0184, 2026-11-26): -999 in every
# variable of padded AIRS daily-file granules; -9999 in AIRS latitude and
# longitude of missing scans (times kept); latitude = longitude = 0 exactly in
# IASI scans without geolocation (some with time 0) and AIRS 2014-01-01 g228.
FILL_VALUES = (-999.0, -9999.0)
STRAY_SECONDS = 600.0  # a pixel this far from its scan's median time is stray
MISDATED_DAYS = 1.0  # a scan this far outside its file's span is misdated


def nadir_fill_to_nan(lat, lon, time, temp, span):
    """
    Fill values and impossible times in one nadir file's arrays to NaN, in
    place. lat, lon, time [along x cross] (time in seconds since 2000-01-01),
    temp [along x cross x level]; span: (first, last) second the file covers.

    - lat or lon -999 or -9999, or both exactly 0: lat and lon NaN;
    - time -999 or exactly 0 (2000-01-01T00:00, before either instrument): NaN;
    - a time more than STRAY_SECONDS from its scan's median: NaN (stray pixel);
    - a scan whose median time is more than MISDATED_DAYS outside span: every
      time of the scan NaN (misdated);
    - temp -999 or -9999: NaN.

    Returns {"position": n, "time": n, "stray": n, "misdated": n scans,
    "temp": n} (points set to NaN, except misdated).
    """
    fill_pos = (np.isin(lat, FILL_VALUES) | np.isin(lon, FILL_VALUES)
                | ((lat == 0) & (lon == 0)))
    lat[fill_pos] = np.nan
    lon[fill_pos] = np.nan
    fill_time = np.isin(time, FILL_VALUES) | (time == 0)
    time[fill_time] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN scans
        median = np.nanmedian(time, axis=1)
        stray = np.abs(time - median[:, None]) > STRAY_SECONDS
        time[stray] = np.nan
        median = np.nanmedian(time, axis=1)
    margin = MISDATED_DAYS * 86400
    misdated = (median < span[0] - margin) | (median > span[1] + margin)
    time[misdated, :] = np.nan
    fill_temp = np.isin(temp, FILL_VALUES)
    temp[fill_temp] = np.nan
    return {"position": int(fill_pos.sum()), "time": int(fill_time.sum()),
            "stray": int(stray.sum()), "misdated": int(misdated.sum()),
            "temp": int(fill_temp.sum())}


def note_fill_counts(inst_info, counts):
    """Add the totals of nadir_fill_to_nan to the load record's notes (AIRS:
    granules returned; IASI: files read, including any padding)."""
    text = {"Lat": ("position", "fill positions (-999, -9999 or 0/0) set to NaN"),
            "Lon": ("position", "fill positions (-999, -9999 or 0/0) set to NaN"),
            "Temp": ("temp", "fill values (-999, -9999) set to NaN")}
    notes = record(inst_info)["notes"]
    for name, (key, words) in text.items():
        if counts.get(key):
            _append_note(notes, name, f"{words}: {counts[key]} points")
    parts = [f"{counts[k]} {words}" for k, words in (
        ("time", "fill times (-999 or 0) set to NaN"),
        ("stray", f"stray times (> {STRAY_SECONDS:g} s from their scan) set to NaN"),
        ("misdated", f"misdated scans (> {MISDATED_DAYS:g} day outside their file) "
                     "set to NaN")) if counts.get(k)]
    if parts:
        _append_note(notes, "Time", "; ".join(parts))


def _append_note(notes, name, text):
    notes[name] = f"{notes[name]}; {text}" if notes.get(name) else text


# range check


def outside_physical_range(temp, alt=None, floor_is_bad=False):
    """
    True where a temperature is physically implausible: below 100 K, or above
    400 K (only below 100 km if heights alt are given, as the MATLAB limb
    loader). With floor_is_bad, exactly 100 K too (a retrieval's lower bound,
    e.g. IASI). NaN is not flagged. alt broadcasts against temp.
    """
    temp = np.asarray(temp)
    with np.errstate(invalid="ignore"):
        hot = temp > 400
        if alt is not None:
            hot = hot & (np.asarray(alt) < 100)
        cold = temp <= 100 if floor_is_bad else temp < 100
        return cold | hot


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
    source_versions, files_not_found and multiple_matches (from
    settings["record"]).
    """
    first, last = (str(datenum_to_datetime(t, "ms")) for t in settings["time_range"])
    return {
        "time_range": f"{first} to {last}",
        "source": source,
        "atmosgw_curate_version": __version__,
        "options": options_json(call) if call is not None else "",
        "files_read": ", ".join(settings["files_read"]),
        "source_versions": versions_text(settings["record"], settings["files_read"]),
        "files_not_found": ", ".join(settings["record"]["files_not_found"]),
        "files_without_data": ", ".join(settings["record"].get("files_without_data",
                                                               [])),
        "multiple_matches": "; ".join(settings["record"]["multiple_matches"]),
    }


def warn_missing_files(settings, stacklevel, what="files"):
    """Warn once per load if nothing was read, or if some of what was looked
    for was not found (settings["record"]["files_not_found"]) or was read but
    held no values ("files_without_data"); silenced by date_warning=False."""
    if not settings.get("date_warning", True):
        return
    rec = settings["record"]
    missing, empty = rec["files_not_found"], rec.get("files_without_data", [])
    if not settings["files_read"] and missing:
        warnings.warn(
            f"None of the {settings['instrument']} {what} looked for were found "
            f"(e.g. {missing[0]}; looked in "
            f"{settings.get('data_path', 'the data folder')}); the Dataset is empty. "
            "Check the data folder and the dates.", stacklevel=stacklevel + 1)
    elif missing:
        warnings.warn(
            f"{len(missing)} {what} looked for were not found (e.g. {missing[0]}); "
            "all are listed in the output attribute 'files_not_found'.",
            stacklevel=stacklevel + 1)
    if empty:
        warnings.warn(
            f"{len(empty)} {settings['instrument']} file(s) read hold no values of "
            f"the main variable (e.g. {empty[0]}); listed in the output attribute "
            "'files_without_data'.", stacklevel=stacklevel + 1)


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


# data versions of the files read


NOT_STATED = "not stated in the files"


def note_version(inst_info, path, text):
    """Record the data version stated in one file read: a text, or a list of
    texts for files holding several (GNSS: one per processing centre's
    retrieval version); "" or [] if it states none."""
    items = [text] if isinstance(text, str) else list(text)
    items = tuple(str(t).strip() for t in items if str(t).strip())
    record(inst_info).setdefault("versions", {})[os.path.basename(str(path))] = items


def file_versions(record, files):
    """The version text of each file in files (NOT_STATED where none)."""
    versions = record.get("versions", {})
    return [", ".join(versions.get(f, ())) or NOT_STATED for f in files]


def versions_text(record, files):
    """The ``source_versions`` attribute: each version stated by the files
    read, with how many files state it, e.g. 'Data_Product_Version 02.07 (30
    files); Data_Product_Version 02.08 (1 file)'; '' if no file was read."""
    versions = record.get("versions", {})
    counts = {}
    for f in files:
        for text in versions.get(f) or (NOT_STATED,):
            counts[text] = counts.get(text, 0) + 1
    return "; ".join(f"{text} ({n} file{'s' if n != 1 else ''})"
                     for text, n in counts.items())


def warn_mixed_versions(settings, stacklevel, exempt=()):
    """Warn once per load if the files read state more than one data version
    (instruments in exempt, whose files always mix versions, are not warned)."""
    if settings.get("instrument_name", settings.get("instrument")) in exempt:
        return
    versions = sorted(set(file_versions(settings["record"], settings["files_read"])))
    if len(versions) > 1:
        warnings.warn(
            f"The {settings['instrument']} files read state {len(versions)} different "
            f"data versions ({'; '.join(versions)}). Their counts are in the output "
            "attribute 'source_versions'.", stacklevel=stacklevel + 1)
