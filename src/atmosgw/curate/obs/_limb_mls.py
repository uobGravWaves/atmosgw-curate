"""
Limb-sounder reader for MLS temperature (NASA Level 2 HDF-EOS5 files,
``Temperature`` swath). Port of ``module_load_MLS.m`` (uobGravWaves/
MultiInstrumentFunctions, commit 78eef78), with the parts of ``get_MLS.m``
(corwin365/MatlabFunctions, commit c943399) that it uses.

Files: one per day, the file below the MLS folder whose path contains
``_<yyyy>d<ddd>``. If several match (e.g. v05-04 and v06-04 of the same day,
as for many days of 2026 in the group's MLS folder), the newest version is
used (``_limb_tools.newest_mls_version``), and the day is listed in the
output's ``multiple_matches`` attribute. The MATLAB loader does the same
since 2026-11-25; before, both took the first match in name order, there the
older version (Q8 in ``docs/quirks.md``).

The fields ``get_MLS`` reads are available: data fields
``L2gpValue`` (temperature), ``L2gpPrecision``, ``Quality``, ``Status``,
``Convergence``; geolocation fields ``ChunkNumber``, ``Latitude``,
``LineOfSightAngle``, ``LocalSolarTime``, ``Longitude``, ``OrbitGeodeticAngle``,
``Pressure``, ``SolarZenithAngle``, ``Time`` (seconds since 1993-01-01).

Each field keeps its stored number type (as ``get_MLS``), except that a field
containing its ``_FillValue`` is converted to float64 with those values NaN.
Negative ``L2gpPrecision`` values are meaningful (a priori influence) and are
kept. Heights are estimated from the pressure levels with
``atmosphere.pressure_to_height`` (U.S. Standard Atmosphere 1976), in double
precision (MLS pressures are float32; until 2026-11-25 the heights were computed
in single precision, as the old MATLAB ``p2h`` did).
"""

import h5py
import numpy as np

from atmosgw.curate.atmosphere import pressure_to_height
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs._limb_tools import (
    choose_file,
    conform_vars,
    day_label,
    loop_days,
    newest_mls_version,
    note,
    note_sources,
    year_and_day,
)
from atmosgw.curate.readers import find_files

_SWATH = "Temperature"
DATA_FIELDS = ("L2gpValue", "L2gpPrecision", "Quality", "Status", "Convergence")
GEOLOCATION_FIELDS = ("ChunkNumber", "Latitude", "LineOfSightAngle", "LocalSolarTime",
                      "Longitude", "OrbitGeodeticAngle", "Pressure", "SolarZenithAngle",
                      "Time")


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    note_sources(inst_info, {"Temp": "L2gpValue", "Lat": "Latitude",
                             "Lon": "Longitude", "Pres": "Pressure", "Time": "Time"})
    note(inst_info, "Alt", "estimated from Pressure with "
                           "atmosphere.pressure_to_height (U.S. Standard "
                           "Atmosphere 1976 [R1])")
    note(inst_info, "Time", "file Time: TAI seconds since 1993-01-01, converted "
                            "without leap seconds")

    def find_file(day):
        year, doy = year_and_day(day)
        found = find_files(inst_info["path"], f"_{year:04d}d{doy:03d}")
        # several versions of a day: the newest (Q8; MATLAB takes the first)
        return choose_file(found, inst_info, day_label(day), newest_mls_version)

    def read_file(path, file_count, day):
        working = _get_mls(path)
        store = {}
        for name in var_names:
            if name == "Temp":
                store["Temp"] = working["L2gpValue"]
            elif name == "Lat":
                store["Lat"] = working["Latitude"]
            elif name == "Lon":
                store["Lon"] = working["Longitude"]
            elif name == "Alt":
                store["Alt"] = pressure_to_height(working["Pressure"])
            elif name == "Pres":
                store["Pres"] = working["Pressure"]
            elif name == "Time":
                store["Time"] = datenum(1993, 1, 1, 0, 0, working["Time"])
            elif name == "SourceProf":
                n = store["Lat"].size
                store["SourceProf"] = np.arange(1.0, n + 1)[:, None]  # 1-based
            elif name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            elif name in working:
                store[name] = working[name]
            else:
                raise KeyError(f"Variable {name} not found in {path}")

        # per-profile and per-level variables onto [profiles x levels]
        n_prof, n_lev = store["Temp"].shape
        for name in ("Lat", "Lon", "Time", "SourceProf", "SourceFile"):
            store[name] = np.tile(store[name], (1, n_lev))
        store["Pres"] = np.tile(store["Pres"].T, (n_prof, 1))
        store["Alt"] = np.tile(store["Alt"].T, (n_prof, 1))
        for name in settings["additional_vars"]:
            field = store[name]
            if field.size == n_lev:
                store[name] = np.tile(field.reshape(1, -1, order="F"), (n_prof, 1))
            elif field.size == n_prof:
                store[name] = np.tile(field.reshape(-1, 1, order="F"), (1, n_lev))
        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, "MLS", inst_info)


def _get_mls(path):
    """The fields get_MLS reads, in their stored types, 1-D fields as columns;
    a field containing its _FillValue becomes float64 with NaN there."""
    working = {}
    with h5py.File(path, "r") as f:
        for group, names in (("Data Fields", DATA_FIELDS),
                             ("Geolocation Fields", GEOLOCATION_FIELDS)):
            for name in names:
                key = f"HDFEOS/SWATHS/{_SWATH}/{group}/{name}"
                if key not in f:
                    working[name] = np.empty((0, 0))
                    continue
                dataset = f[key]
                value = dataset[()]
                fill = dataset.attrs.get("_FillValue")
                if fill is not None:
                    bad = value == np.asarray(fill, dtype=value.dtype).ravel()[0]
                    if bad.any():
                        value = value.astype(np.float64)
                        value[bad] = np.nan
                working[name] = value[:, None] if value.ndim == 1 else value
    return working
