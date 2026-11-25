"""
Nadir-sounder reader for IASI 3-D temperature retrievals, in the group's
netCDF orbit files, with the IASI outlier screen.

Port of ``get_iasi_granule`` in ``prep_airs_3d.m`` (corwin365/MatlabFunctions,
commit f1a78c0, which includes the November 2026 fixes).

Files, below the IASI folder (``datasets.toml``, ``IASI`` ``raw``):
``YYYY/DDD/iasi_YYYYMMDDhhmmss.nc``, one per orbit dump, named by the start
of its first scan. Only names whose date matches the folder's day are used
(the 2007/149 folder also holds 2008-12-24 files). Files are back-to-back
(the next starts one 8-s scan after the last ends; checked at 67 boundaries,
2007-2016), so a file covers [its start, the next file's start). Files of the
day before and after are listed too, as orbits cross midnight.

Pseudo-granules (MATLAB): granule k of a day holds the scans whose first
cross-track position's time lies in [6(k-1), 6k) minutes after 00:00 UTC
(computed in seconds, so the bin edges are exact). A granule's time is the
mean of its times; it belongs to the time range if that lies in it, as for
AIRS. Granules can also be given explicitly.

Variables (reader contract, ``docs/design/nadir_loader.md`` section 5), as
for AIRS: ``Lat``, ``Lon``, ``Time`` (datenums), ``Temp`` [cross x along x
level], ``Alt``, and per row ``Granule``, ``Day``, ``SourceFile``,
``SourceRow`` (1-based row within its granule). 60 cross-track positions.

Differences from MATLAB, by design:

- files without temperature (e.g. all of 2007-05-29 and 2007-11-08 on
  eepc-0184) are skipped and recorded ("<file>: no ret_temp"); MATLAB stops
  with an error;
- the outlier screen (now the swappable ``obs._screen``, rolling by default)
  marks points but does not refill them: MATLAB
  refills screened points from the nearest along-track value; here they are
  NaN with their state, and filling is a separate, visible step (agreed);
- MATLAB's warnings when more than 5 % (``IasiFracBad``) or a connected 3 %
  (``IasiRegionFracBad``) of a granule is bad are replaced by the state
  counts in the output attributes.
"""

import os
import re
import warnings

import numpy as np

from atmosgw.curate.matlab_time import datenum_to_datetime, datetime_to_datenum
from atmosgw.curate.obs._common import (
    day_label,
    note,
    note_sources,
    record,
    year_and_day,
)
from atmosgw.curate.obs._nadir_airs import J2000
from atmosgw.curate.readers import read_netcdf

N_CROSS = 60
GRANULES_PER_DAY = 240
_NAME = re.compile(r"^iasi_(\d{14})\.nc$")
FILE_NAMES = {"Lat": "l1_lat", "Lon": "l1_lon", "Temp": "ret_temp", "Alt": "ret_z"}


def list_files(root, first_day, last_day):
    """(start datenum, path) of every file for the days, in time order."""
    out = []
    for day in range(int(first_day), int(last_day) + 1):
        t = datenum_to_datetime(day)
        year, doy = year_and_day(day)
        folder = os.path.join(root, f"{year:04d}", f"{doy:03d}")
        if not os.path.isdir(folder):
            continue
        stamp = str(t)[:10].replace("-", "")
        for name in sorted(os.listdir(folder)):
            match = _NAME.match(name)
            if match and match.group(1).startswith(stamp):
                s = match.group(1)
                start = datetime_to_datenum(np.datetime64(
                    f"{s[:4]}-{s[4:6]}-{s[6:8]}T{s[8:10]}:{s[10:12]}:{s[12:14]}"))
                out.append((float(start), os.path.join(folder, name)))
    return sorted(out)


def granule_of(time_seconds):
    """(day datenum, pseudo-granule number) of times in seconds since 2000."""
    t = np.asarray(time_seconds, dtype=float)
    day = J2000 + np.floor(t / 86400)
    return day, np.floor(np.mod(t, 86400) / 360) + 1


def load(settings, inst_info, var_names):
    """Read the granules; returns (data, file_list). See the module docstring."""
    note_sources(inst_info, {**FILE_NAMES, "Time": "l1_time"})
    note(inst_info, "Time", "l1_time: seconds since 2000-01-01T00:00Z")
    extra = [n for n in var_names if n not in ("Lat", "Lon", "Time", "Temp", "Alt")]
    explicit = settings.get("granules")
    if explicit is not None:
        wanted = {(float(d), float(g)) for d, g in explicit}
        spans = [(d + (g - 1) / 240, d + g / 240) for d, g in wanted]
    else:
        t0, t1 = settings["time_range"]
        spans = [(t0 - 6 / 1440, t1 + 6 / 1440)]
    low, high = min(s[0] for s in spans), max(s[1] for s in spans)

    files = list_files(inst_info["path"], np.floor(low) - 1, np.floor(high) + 1)
    starts = np.array([s for s, _ in files])
    finish = np.append(starts[1:], np.inf) if files else starts
    use = [k for k in range(len(files))
           if any(starts[k] < b and finish[k] > a for a, b in spans)]

    rows, file_list, alt = [], [], None
    rec = record(inst_info)
    for k in use:
        path = files[k][1]
        try:
            values = read_netcdf(path, ["l1_time", "l1_lat", "l1_lon", "ret_temp"]
                                 + extra)
        except KeyError as error:
            rec["files_not_found"].append(f"{os.path.basename(path)}: no {error}")
            continue
        z = read_netcdf(path, ["ret_z"])["ret_z"].ravel()
        if alt is None:
            alt = z
        elif not np.array_equal(alt, z):
            raise ValueError(f"{os.path.basename(path)}: heights differ from those "
                             "of the files before it.")
        day, gran = granule_of(values["l1_time"][:, 0])
        file_list.append(os.path.basename(path))
        rows.append((values, day, gran, len(file_list)))

    # assign rows to granules and choose granules
    if not rows:
        empty = np.empty((N_CROSS, 0))
        n_lev = 27 if alt is None else alt.size
        return ({"Lat": empty, "Lon": empty, "Time": empty,
                 "Temp": np.empty((N_CROSS, 0, n_lev)), "Granule": np.empty(0),
                 "Day": np.empty(0), "SourceFile": np.empty(0),
                 "SourceRow": np.empty(0),
                 "Alt": np.arange(0, 91, 3.0)[:27] if alt is None else alt},
                file_list)
    time = np.concatenate([J2000 + v["l1_time"] / 86400 for v, *_ in rows], axis=0)
    day = np.concatenate([d for _, d, _, _ in rows])
    gran = np.concatenate([g for _, _, g, _ in rows])
    source = np.concatenate([np.full(d.size, float(n)) for _, d, _, n in rows])
    key_id = day * 1000 + gran  # one number per (day, granule)
    ids, inverse = np.unique(key_id, return_inverse=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        row_mean = np.nanmean(time, axis=1)
        sums = np.bincount(inverse, weights=np.nan_to_num(row_mean))
        counts = np.bincount(inverse, weights=np.isfinite(row_mean).astype(float))
        mean = sums / counts
    if explicit is not None:
        wanted_ids = {d * 1000 + g for d, g in wanted}
        ok = np.isin(ids, list(wanted_ids))
        for d, g in sorted(wanted):
            if d * 1000 + g not in set(ids[ok].tolist()):
                rec["files_not_found"].append(f"{day_label(d)} granule {int(g)}")
    else:
        ok = (mean >= settings["time_range"][0]) & (mean <= settings["time_range"][1])
    order = np.flatnonzero(ok[inverse])
    order = order[np.argsort(time[order, 0], kind="stable")]

    def stack(name):
        return np.concatenate([v[name] for v, *_ in rows], axis=0)[order]

    data = {"Lat": stack("l1_lat").T, "Lon": stack("l1_lon").T, "Time": time[order].T,
            "Temp": np.transpose(stack("ret_temp"), (1, 0, 2)),
            "Granule": gran[order], "Day": day[order], "SourceFile": source[order],
            "Alt": alt}
    # row within its granule (1-based): rows are in time order
    inv = inverse[order]
    row = np.zeros(order.size)
    for i in np.unique(inv):
        sel = np.flatnonzero(inv == i)
        row[sel] = np.arange(1, sel.size + 1)
    data["SourceRow"] = row
    for name in extra:
        v = stack(name)
        data[name] = np.transpose(v, (1, 0) + tuple(range(2, v.ndim)))
    # keep only files that contributed rows
    used = sorted({int(n) for n in data["SourceFile"]})
    renumber = {old: new for new, old in enumerate(used, start=1)}
    data["SourceFile"] = np.array([renumber[int(n)] for n in data["SourceFile"]],
                                  dtype=float)
    return data, [file_list[n - 1] for n in used]

