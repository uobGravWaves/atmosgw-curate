"""
Reader for the ALOMAR lidar in the time_height format (see ``time_height.py``).

Files: as ``_limb_alomar.py`` (``<ALOMAR folder>/<yyyy>_<version>_file4juli.h5``,
one per year, one HDF5 group per measurement session, named by the date the
session started; ``T1_lidar``, ``T1_lidar_unc`` ... as [times x heights],
``time_lidar`` in POSIX seconds, ``z_lidar`` in km). The file name, product
and site position are taken from ``_limb_alomar`` so the two readers cannot
drift apart; ``_limb_alomar.py`` itself is not changed.

Differences from the limb reader (new code, not a port):

- sessions are named by their start date and can last several days (up to
  8.05 days in the files of 2003-2026: 2005-06-30 to 2005-07-07), so rows
  are selected by their own times in every session of the years touched and
  of the year before (a session that started in December can run into the
  range); the previous year's file is listed in the output only if it
  contributed rows. The limb reader reads only the years in the range;
- all rows are put on one height grid: every session's heights lie on a
  0.15 km lattice (5.1, 5.25 ... 120 km; checked on all files 2003-2026 to
  within 1e-11 km), so the union of the sessions' grids holds every value
  exactly, without interpolation. A height more than 1e-6 km off the lattice
  raises ValueError;
- times are kept to the microsecond (the limb loader rounds to the
  millisecond).

Returns ``(rows, file_list)``. rows is a dict:

- ``time_posix`` (n,), float64 POSIX seconds; ``time_datenum`` (n,);
- ``level_index`` (n_level,), int: lattice index k of each output level
  (height k * 0.15 km), sorted;
- ``T1_lidar``, ``T1_lidar_unc`` and each additional variable: (n, n_level)
  on that grid (NaN where a session has no such height) or (n,);
- ``present`` (n, n_level) bool: the session has that height;
- ``source_profile`` (n,): 1-based row number through the yearly file, in
  group order, counting every group (as ``SourceProf`` of the limb reader);
- ``source_file`` (n,): 1-based index into file_list;
- ``session`` (n,): 1-based session number in the output, in reading order.

Rows are in file order (not yet sorted by time) and only those inside
settings["time_range"] (datenums, inclusive, as the limb loader) are kept.
"""

import os

import h5py
import numpy as np

from atmosgw.curate.matlab_time import datenum_to_datetime
from atmosgw.curate.obs._common import note, note_sources, note_version, record
from atmosgw.curate.obs._limb_alomar import DATA_VERSION, RESOLUTION

LEVEL_STEP = 0.15  # km, the lattice of every session's heights
LATTICE_TOLERANCE = 1e-6  # km
_POSIX_EPOCH = 719529  # datenum of 1970-01-01


def load(settings, inst_info):
    """Read the sessions touching settings["time_range"]. See the module docstring."""
    temp_name, unc_name = f"{RESOLUTION}_lidar", f"{RESOLUTION}_lidar_unc"
    extra = list(dict.fromkeys(n for n in settings["additional_vars"]
                               if n not in (temp_name, unc_name)))
    for name in ("time_lidar", "z_lidar"):  # (already the time and altitude)
        if name in extra:
            raise ValueError(f"{name!r} is not an additional variable: it gives "
                             "the output's time or altitude.")
    start, end = (float(t) for t in settings["time_range"])
    years = range(_year(start) - 1, _year(end) + 1)  # (the year before: see above)

    note_sources(inst_info, {"Temp": temp_name, "TempUnc": unc_name,
                             "Time": "time_lidar", "Alt": "z_lidar"})
    note(inst_info, "Lat", "site latitude (fixed)")
    note(inst_info, "Lon", "site longitude (fixed)")

    sessions, file_list = [], []
    for year in years:
        file_name = os.path.join(inst_info["path"],
                                 f"{year}_{DATA_VERSION}_file4juli.h5")
        look_back_only = year < _year(start)
        if not os.path.isfile(file_name):
            if not look_back_only:
                record(inst_info)["files_not_found"].append(str(year))
            continue
        file_index = None
        note_version(inst_info, file_name, f"{DATA_VERSION} (file name)")
        if not look_back_only:  # listed, as the limb reader, even if no rows
            file_list.append(os.path.basename(file_name))
            file_index = len(file_list)
        rows_before = 0  # rows of the earlier groups in this file
        with h5py.File(file_name, "r") as f:
            for group in sorted(f.keys()):
                g = f[group]
                n_time = g["time_lidar"].shape[0]
                first_row = rows_before
                rows_before += n_time
                t = np.asarray(g["time_lidar"][()], dtype=np.float64).ravel()
                t_datenum = t / 86400 + _POSIX_EPOCH  # (as the limb reader)
                keep = np.flatnonzero((t_datenum >= start) & (t_datenum <= end))
                if keep.size == 0:
                    continue
                if settings["verbose"]:
                    print(f"----> Loading ALOMAR data: {group}")
                if file_index is None:  # the look-back file, now that it is used
                    file_list.append(os.path.basename(file_name))
                    file_index = len(file_list)
                k = _lattice_index(np.asarray(g["z_lidar"][()], dtype=np.float64),
                                   group)
                values = {}
                for name in [temp_name, unc_name] + extra:
                    if name not in g:
                        raise KeyError(f"Variable {name!r} is not in ALOMAR "
                                       f"session {group} "
                                       f"({os.path.basename(file_name)}).")
                    value = np.asarray(g[name][()], dtype=np.float64)
                    if value.shape in ((n_time, k.size), (n_time,)):
                        values[name] = value[keep]
                    else:
                        raise ValueError(
                            f"additional variable {name!r} has shape {value.shape} in "
                            f"session {group}: only [times x heights] ({n_time} x "
                            f"{k.size}) or [times] ({n_time}) can be loaded.")
                sessions.append({"t": t[keep], "t_datenum": t_datenum[keep], "k": k,
                                 "values": values,
                                 "source_profile": first_row + keep + 1,
                                 "source_file": file_index})
    return _on_one_grid(sessions, [temp_name, unc_name] + extra), file_list


# supporting functions


def _year(day_number):
    return int(str(datenum_to_datetime(day_number))[:4])


def _lattice_index(z, group):
    """Lattice index k (height k * LEVEL_STEP) of each height; ValueError if
    any height is off the lattice or repeated."""
    k = np.round(z / LEVEL_STEP)
    off = np.abs(z - k * LEVEL_STEP)
    if not np.all(np.isfinite(z)) or np.any(off > LATTICE_TOLERANCE):
        raise ValueError(f"ALOMAR session {group}: heights are not all on the "
                         f"{LEVEL_STEP} km grid (largest offset {np.nanmax(off):.3g} "
                         "km), so they cannot be put on one grid.")
    k = k.astype(np.int64)
    if np.unique(k).size != k.size:
        raise ValueError(f"ALOMAR session {group}: repeated heights.")
    return k


def _on_one_grid(sessions, names):
    """Concatenate the sessions' rows on the union of their height grids."""
    level_index = (np.unique(np.concatenate([s["k"] for s in sessions]))
                   if sessions else np.zeros(0, dtype=np.int64))
    n = sum(s["t"].size for s in sessions)
    n_level = level_index.size
    rows = {
        "time_posix": np.concatenate([s["t"] for s in sessions]) if sessions
        else np.zeros(0),
        "time_datenum": np.concatenate([s["t_datenum"] for s in sessions])
        if sessions else np.zeros(0),
        "level_index": level_index,
        "present": np.zeros((n, n_level), dtype=bool),
        "source_profile": np.zeros(n, dtype=np.int64),
        "source_file": np.zeros(n, dtype=np.int64),
        "session": np.zeros(n, dtype=np.int64),
    }
    for name in names:
        n_dims = {s["values"][name].ndim for s in sessions}
        if len(n_dims) > 1:
            raise ValueError(f"variable {name!r} is per time in some sessions "
                             "and per point in others.")
        per_time = n_dims == {1}
        rows[name] = np.full((n,) if per_time else (n, n_level), np.nan)
    row = 0
    for number, s in enumerate(sessions, start=1):
        rows_here = slice(row, row + s["t"].size)
        columns = np.searchsorted(level_index, s["k"])
        rows["present"][rows_here, columns] = True
        rows["source_profile"][rows_here] = s["source_profile"]
        rows["source_file"][rows_here] = s["source_file"]
        rows["session"][rows_here] = number
        for name in names:
            value = s["values"][name]
            if value.ndim == 1:
                rows[name][rows_here] = value
            else:
                rows[name][rows_here, columns] = value
        row += s["t"].size
    return rows
