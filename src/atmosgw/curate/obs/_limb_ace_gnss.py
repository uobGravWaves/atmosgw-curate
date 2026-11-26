"""
Limb-sounder reader for ACE-FTS, GNSS radio occultation, and "Misc" files
(any files in the same MATLAB format). Port of ``module_load_ACE_GNSS.m``
(uobGravWaves/MultiInstrumentFunctions, commit 78eef78).

File format: one MATLAB .mat file per day, holding [profiles x levels] arrays
named as in the output (``Lat``, ``Lon``, ``Temp``, and for ACE and Misc
``Alt`` and ``Time``). GNSS files instead hold ``MSL_alt`` (one row of
heights) and a table ``MetaData`` with a ``time`` column (one per profile).
If there is no ``Pres``, it is estimated from ``Alt`` with
``height_to_pressure``.

File names: ACE ``<path>/<yyyy>/ace_v52_<yyyy>d<ddd>.mat``; GNSS
``<path>/<yyyy>/merged_ro_<yyyy>_<ddd>.mat``; Misc: the first file below
``misc_path`` matching ``*<identifier>*<yyyy>d<ddd>*.mat``.

Units and corrections (not in the MATLAB ``module_load_ACE_GNSS.m``, which
reads the files as they are; Q9-Q11 in ``docs/quirks.md``):

- ACE ``Pres`` is in atmospheres and GNSS ``Pres`` in pascals (median ratio
  to the U.S. Standard Atmosphere at 12-40 km: ACE 0.00088-0.0010 in every
  year 2004-2025, 1/1013.25 = 0.000987; GNSS 97-103 in every year
  2001-2025), so both are converted to hPa. Misc files are read as they are.
- GNSS (``gnss_fixes``): the netCDF default fill value 9.97e36 (in
  ``Temp``, ``Pres`` and ``Azim`` of most years) becomes NaN. Profiles with
  no position at any level take their ``MetaData`` latitude and longitude at
  every level (2-46 % of profiles on the days checked, 2002-2025, nearly all
  JPL retrievals). Where a profile's ``MetaData`` time is more than 10
  minutes from the time in its occultation id (``...-yyyymmddHHMM``), or is
  exactly on the hour and more than a minute from it, the id's time is used
  (JPL times in the files are rounded down to the hour: median 30 min, up
  to 62 min early; other centres agree with the id to within 6 min).
- GNSS: an occultation present more than once (the same id from several
  processing centres: 2-64 % of profiles on the days checked) is returned
  once, unless ``keep_duplicate_occultations``: of the copies with per-level
  positions and a good time (if any), the finest original vertical
  resolution (smallest median ``PTRes``: UCAR about 0.1, JPL 0.4, ROM SAF
  1.2, units not stated in the files), the first in the file if equal
  (decided by C. Wright, 2026-11-26). In practice UCAR or ROM SAF is kept
  over JPL, and UCAR over ROM SAF.
"""

import os
import re

import numpy as np

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.matlab_time import datetime_to_datenum
from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    note,
    note_sources,
    note_version,
    record,
    year_and_day,
)
from atmosgw.curate.obs._limb_tools import conform_vars, loop_days, note_estimated_pres
from atmosgw.curate.readers import find_files, read_mat


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    instrument = settings["instrument"]
    sources = {name: name for name in ("Lat", "Lon", "Temp", "Pres", "Alt", "Time")}
    if instrument == "GNSS":
        sources.update(Alt="MSL_alt", Time="MetaData.time")
    note_sources(inst_info, sources)

    def find_file(day):
        year, doy = year_and_day(day)
        if instrument == "GNSS":
            return os.path.join(
                inst_info["path"], f"{year:04d}", f"merged_ro_{year:04d}_{doy:03d}.mat"
            )
        if instrument == "ACE":
            return os.path.join(
                inst_info["path"], f"{year:04d}", f"ace_v52_{year:04d}d{doy:03d}.mat"
            )
        found = find_files(
            inst_info["path"],
            f"*{inst_info['identifier']}*{year:04d}d{doy:03d}*.mat",
        )
        return choose_file(found, inst_info, day_label(day))

    def read_file(path, file_count, day):
        inst_data = read_mat(path)
        rows = np.arange(np.shape(inst_data["Lat"])[0])  # file rows (0-based)
        if instrument == "GNSS":
            rows = gnss_fixes(inst_data, inst_info,
                              settings.get("keep_duplicate_occultations", False))
        note_version(inst_info, path, file_version(instrument, path, inst_data))
        needed = (("Lat", "Lon", "Temp", "MSL_alt", "MetaData")
                  if instrument == "GNSS" else ("Lat", "Lon", "Temp", "Alt", "Time"))
        for name in needed:
            if name not in inst_data and name in var_names + ["MSL_alt", "MetaData"]:
                raise KeyError(f"Variable {name} not found in {path}")
        if instrument in PRES_TO_HPA and "Pres" in inst_data:
            factor, text = PRES_TO_HPA[instrument]
            inst_data["Pres"] = np.asarray(inst_data["Pres"], dtype=float) * factor
            note(inst_info, "Pres", text)
        estimated = "Pres" not in inst_data
        if estimated:
            inst_data["Pres"] = height_to_pressure(inst_data["Alt"])
            note_estimated_pres(inst_info)
        # which points have estimated pressure (asked for by the core, for the
        # output's data states; not in the MATLAB)
        inst_data["PresEstimated"] = np.full(np.shape(inst_data["Pres"]),
                                             float(estimated))

        store = {}
        for name in var_names:
            computed = ("Alt", "Time", "SourceProf", "SourceFile")
            if name not in computed and name not in inst_data:
                raise KeyError(f"Variable {name} not found in {path}")
            if name == "Time":
                if instrument == "GNSS":
                    time = np.asarray(inst_data["MetaData"]["time"], dtype=float)
                    n_lev = np.shape(inst_data["Lat"])[1]
                    store["Time"] = np.tile(time.reshape(-1, 1), (1, n_lev))
                else:
                    store["Time"] = inst_data["Time"]
            elif name == "Alt":
                if instrument == "GNSS":
                    n_prof = np.shape(inst_data["Lat"])[0]
                    store["Alt"] = np.tile(inst_data["MSL_alt"], (n_prof, 1))
                else:
                    store["Alt"] = inst_data["Alt"]
            elif name == "SourceProf":
                n_lev = np.shape(store["Lat"])[1]
                profile = (rows + 1).astype(np.float32)  # 1-based row in the file
                store["SourceProf"] = np.tile(profile.reshape(-1, 1), (1, n_lev))
            elif name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            else:
                store[name] = inst_data[name]
        shape = np.shape(store["Lat"])
        for name in ("Alt", "Pres", "PresEstimated"):  # (one row: per level)
            if name in store and np.shape(store[name]) == (1,) + shape[1:]:
                store[name] = np.tile(store[name], (shape[0], 1))
        for name in ("Lon", "Temp", "Pres", "Alt", "Time", "PresEstimated"):
            if name in store and np.shape(store[name]) != shape:
                raise ValueError(f"{path}: {name} has shape {np.shape(store[name])}, "
                                 f"not that of Lat {shape} [profiles x levels].")
        return conform_vars(store, settings)

    out = loop_days(settings, var_names, find_file, read_file, instrument, inst_info)
    note_gnss_counts(inst_info)
    return out


# pressure in the files -> hPa (see the module docstring)
PRES_TO_HPA = {
    "ACE": (1013.25, "file Pres in atmospheres x 1013.25 (hPa)"),
    "GNSS": (0.01, "file Pres in Pa / 100 (hPa)"),
}

GNSS_FILL = 9.9e36  # values this large are the netCDF default fill (9.97e36)
ID_TIME_MINUTES = 10.0  # a MetaData time this far from its id's time is wrong
_ID_TIME = re.compile(r"(\d{12})$")  # yyyymmddHHMM at the end of an occid


def gnss_fixes(inst_data, inst_info, keep_duplicates=False):
    """
    The GNSS corrections of the module docstring, in place on one file's
    variables (read_mat output, [profiles x levels] arrays and the MetaData
    table): fill values to NaN, MetaData time from the occultation id where
    the file's is wrong (and the Time array, if the file has one, as the PW
    files), positions from MetaData for profiles with none, and one copy of
    each occultation. Returns the 0-based file rows kept, in file order.
    """
    n_prof = np.shape(inst_data["Lat"])[0]
    for name, value in list(inst_data.items()):
        value = np.asarray(value) if not hasattr(value, "columns") else value
        if isinstance(value, np.ndarray) and value.dtype.kind == "f":
            fill = np.abs(value) >= GNSS_FILL
            if fill.any():
                value = value.astype(float)
                value[fill] = np.nan
                inst_data[name] = value
    meta = inst_data.get("MetaData")
    if meta is None or len(meta) != n_prof:
        return np.arange(n_prof)
    meta = meta.copy()
    for column in ("latitude", "longitude", "time"):  # (fill values here too)
        if column in meta:
            value = np.asarray(meta[column], dtype=float)
            meta[column] = np.where(np.abs(value) >= GNSS_FILL, np.nan, value)
    meta["occid"] = meta["occid"].fillna("").astype(str)
    inst_data["MetaData"] = meta

    # times: the id's where the file's is more than ID_TIME_MINUTES away
    file_time = np.asarray(meta["time"], dtype=float)
    id_time = np.full(n_prof, np.nan)
    for k, occ in enumerate(meta["occid"]):
        found = _ID_TIME.search(occ)
        if found:
            stamp = found.group(1)
            try:
                id_time[k] = datetime_to_datenum(np.datetime64(
                    f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}:"
                    f"{stamp[10:12]}", "m"))
            except ValueError:  # (digits that are not a date: no id time)
                pass
    def on_hour(t):  # (to within 0.4 ms)
        return np.abs(t * 24 - np.round(t * 24)) < 1e-7

    with np.errstate(invalid="ignore"):
        minutes = np.abs(file_time - id_time) * 1440
        off = np.isfinite(id_time) & (
            (minutes > ID_TIME_MINUTES) | ~np.isfinite(file_time)
            | (on_hour(file_time) & ~on_hour(id_time) & (minutes > 1)))
    if off.any():
        meta["time"] = np.where(off, id_time, file_time)
        if "Time" in inst_data:  # (PW files: one time per point)
            time = np.asarray(inst_data["Time"], dtype=float).copy()
            fix = off[:, None] & np.isfinite(time)  # (levels with a time only)
            time[fix] = np.broadcast_to(id_time[:, None], time.shape)[fix]
            inst_data["Time"] = time
        note(inst_info, "Time", "GNSS profiles whose MetaData time is more than "
                                f"{ID_TIME_MINUTES:g} min from the time in their "
                                "occultation id, or on the hour and more than a "
                                "minute from it (JPL: rounded down to the hour in "
                                "the files), "
                                "take the id's time (to the minute)")

    # positions: MetaData's at every level, for profiles with none
    lat = np.asarray(inst_data["Lat"], dtype=float)
    lon = np.asarray(inst_data["Lon"], dtype=float)
    has_position = (np.isfinite(lat) & np.isfinite(lon)).any(axis=1)
    meta_lat = np.asarray(meta["latitude"], dtype=float)
    meta_lon = np.asarray(meta["longitude"], dtype=float)
    filled = ~has_position & np.isfinite(meta_lat) & np.isfinite(meta_lon)
    if filled.any():
        lat, lon = lat.copy(), lon.copy()
        lat[filled] = meta_lat[filled, None]
        lon[filled] = meta_lon[filled, None]
        inst_data["Lat"], inst_data["Lon"] = lat, lon

    # one copy of each occultation
    occid = meta["occid"].to_numpy()
    ids, counts = np.unique(occid[occid != ""], return_counts=True)
    if keep_duplicates or not np.any(counts > 1):
        _count_note(inst_info, "gnss_metadata_positions", int(filled.sum()))
        return np.arange(n_prof)
    resolution = np.full(n_prof, np.inf)
    if "PTRes" in inst_data:
        with np.errstate(invalid="ignore"):
            res = np.asarray(inst_data["PTRes"], dtype=float)
            finite = np.isfinite(res).any(axis=1)
            resolution[finite] = np.nanmedian(res[finite], axis=1)
    good = has_position & ~off
    keep = np.ones(n_prof, dtype=bool)
    for occ in ids[counts > 1]:
        copies = np.flatnonzero(occid == occ)
        pool = copies[good[copies]] if good[copies].any() else copies
        best = pool[np.argmin(resolution[pool])]  # (argmin: first of equal)
        keep[copies] = False
        keep[best] = True
    rows = np.flatnonzero(keep)
    for name, value in list(inst_data.items()):
        if name in ("MSL_alt", "Settings", "Note"):
            continue
        if hasattr(value, "iloc") and len(value) == n_prof:
            inst_data[name] = value.iloc[rows].reset_index(drop=True)
        elif isinstance(value, np.ndarray) and value.ndim >= 1 \
                and value.shape[0] == n_prof:
            inst_data[name] = value[rows]
        elif np.ndim(value) >= 1 and n_prof in np.shape(value)[1:]:
            raise ValueError(
                f"GNSS variable {name!r} has shape {np.shape(value)}, with profiles "
                "not along its first dimension, so repeated occultations cannot be "
                "removed from it; use keep_duplicate_occultations=True.")
    _count_note(inst_info, "gnss_metadata_positions", int(filled[rows].sum()))
    _count_note(inst_info, "gnss_duplicates_removed", int(n_prof - rows.size))
    return rows


def _count_note(inst_info, key, n):
    """Add n to a running count in the load record (summed over files)."""
    counts = record(inst_info).setdefault("counts", {})
    counts[key] = counts.get(key, 0) + n


def note_gnss_counts(inst_info):
    """The counts of gnss_fixes, summed over the files read, as notes."""
    counts = record(inst_info).get("counts", {})
    n = counts.get("gnss_metadata_positions", 0)
    if n:
        for name in ("Lat", "Lon"):
            note(inst_info, name, f"{n} profiles of the files read (before time and "
                                  "box selection) have no position in the file and "
                                  "take their MetaData position at every level")
    n = counts.get("gnss_duplicates_removed", 0)
    if n:
        for name in ("Temp", "Temp_Residual", "Temp_PW"):
            note(inst_info, name, f"{n} repeated copies of occultations (the same "
                                  "occultation id from several processing centres) "
                                  "left out of the files read: per occultation, the "
                                  "copy with per-level positions and a good time, "
                                  "then the finest original resolution (PTRes), is "
                                  "kept")


_RETRIEVAL = re.compile(r"atmosphericRetrieval_[^_]+_([A-Za-z]+)_([^_]+)_")


def file_version(instrument, path, inst_data):
    """The data version a file states (see note_version): ACE from its name;
    GNSS, a list: each processing centre's retrieval version in the
    source-file names of the profiles returned, and the merge version if the
    file records one; Misc none."""
    name = os.path.basename(str(path))
    if instrument == "ACE":
        found = re.search(r"ace_(v\d+)_", name)
        return f"{found.group(1)} (file name)" if found else ""
    if instrument != "GNSS" or "MetaData" not in inst_data:
        return ""
    meta = inst_data["MetaData"]
    pairs = set()
    if "AtmSource" in meta:
        for source in meta["AtmSource"].astype(str):
            found = _RETRIEVAL.search(source.replace("\\", "/").rsplit("/", 1)[-1])
            if found:
                pairs.add(f"{found.group(1)} {found.group(2)}")
    items = [f"{pair} retrieval" for pair in sorted(pairs)]
    settings = inst_data.get("Settings")
    if isinstance(settings, dict) and np.size(settings.get("Version", [])):
        items.append(f"merge {np.asarray(settings['Version']).ravel()[0]}")
    return items
