"""
Limb-sounder reader for HIRDLS (NASA Level 2 HDF-EOS5 files). Port of
``module_load_HIRDLS.m`` and its ``get_HIRDLS`` (uobGravWaves/
MultiInstrumentFunctions, commit 78eef78).

Files: one per day, the first file below the HIRDLS folder whose path matches
``*_<yyyy>d<ddd>*``. Variables are read from
``HDFEOS/SWATHS/HIRDLS/Geolocation Fields/`` or ``.../Data Fields/`` as float64,
with the HIRDLS fill value -999 set to NaN. Times are seconds since
1993-01-01; altitudes are converted from m to km.

Geolocation correction [R2]: HIRDLS gives each profile's position at 30 km,
but the measurement point moves about 0.6 km along track for each km of
vertical scan. Each point is moved by 0.6 * (z - 30 km) along the direction
from its profile to the next (at the level nearest 30 km on average), forwards
for upward scans and backwards for downward scans (``ScanUpFlag``; the
direction rule is the group's own reasoning, see ``docs/theory/relations.md``).
"""

import numpy as np

from atmosgw.curate.geodesy import azimuth, km2deg, reckon
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    note,
    note_sources,
    note_version,
    record,
    year_and_day,
)
from atmosgw.curate.obs._limb_mls import pge_version
from atmosgw.curate.obs._limb_tools import conform_vars, loop_days
from atmosgw.curate.readers import find_files, read_hdf5

_SWATH = "HDFEOS/SWATHS/HIRDLS"

# variables of an HDF-EOS5 HIRDLS file, by group (as listed in get_HIRDLS)
GEOLOCATION_FIELDS = (
    "Altitude", "Latitude", "LocalSolarTime", "Longitude", "OrbitAscendingFlag",
    "OrbitNumber", "Pressure", "ProfileID", "ScanAzimuthAtNominalAltitude",
    "ScanElevationAtNominalAltitude", "ScanTable", "ScanUpFlag", "ScienceScanMode",
    "SecondsInDay", "SolarZenithAngle", "SpacecraftAltitude", "SpacecraftLatitude",
    "SpacecraftLongitude", "TangentHeightAtNominalAltitude", "Time",
    "ViewDirectionAtNominalAltitude",
)
DATA_FIELDS = tuple(
    [f"{band}Micron{part}" for band in ("10.8", "12.1", "17.4", "7.1", "8.3")
     for part in ("CloudAerosolFlag", "Extinction", "ExtinctionNormChiSq",
                  "ExtinctionPrecision", "ExtinctionQuality")]
    + [f"{gas}{part}" for gas in ("CFC11", "CFC12", "CH4", "ClONO2", "H2O", "HNO3",
                                  "N2O", "N2O5", "NO2", "O3")
       for part in ("", "NormChiSq", "Precision", "Quality")]
    + ["CloudTopPressure", "GPH", "GPHPrecision", "RawGPH", "RawGPHPrecision",
       "Temperature", "TemperatureNormChiSq", "TemperaturePrecision",
       "TemperatureQuality"]
)


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    note_sources(inst_info, {"Temp": "Temperature", "Lat": "Latitude",
                             "Lon": "Longitude", "Alt": "Altitude",
                             "Pres": "Pressure", "Time": "Time"})
    note(inst_info, "Lat", "moved along track by 0.6 km per km from 30 km [R2]")
    note(inst_info, "Lon", "moved along track by 0.6 km per km from 30 km [R2]")
    note(inst_info, "Alt", "file Altitude (m) / 1000")
    note(inst_info, "Time", "file Time: TAI seconds since 1993-01-01, converted "
                            "without leap seconds")

    def find_file(day):
        year, doy = year_and_day(day)
        found = find_files(inst_info["path"], f"*_{year:04d}d{doy:03d}*")
        return choose_file(found, inst_info, day_label(day))

    def read_file(path, file_count, day):
        def get(name):
            return _get_hirdls(path, name)

        # kept apart from the variables, so a requested ScanUpFlag keeps its values
        scan_sign = get("ScanUpFlag")
        note_version(inst_info, path, pge_version(path))

        store = {}
        for name in var_names:
            if name == "Temp":
                store["Temp"] = get("Temperature")
            elif name == "Lat":
                store["Lat"] = get("Latitude")
            elif name == "Lon":
                store["Lon"] = get("Longitude")
            elif name == "Alt":
                store["Alt"] = get("Altitude") / 1000
            elif name == "Pres":
                store["Pres"] = get("Pressure")
            elif name == "Time":
                store["Time"] = datenum(1993, 1, 1, 0, 0, get("Time"))
            elif name == "SourceProf":
                n = store["Lat"].size
                store["SourceProf"] = np.arange(1, n + 1, dtype=np.float32)[:, None]
            elif name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            else:
                store[name] = get(name)

        # per-profile and per-level variables onto [profiles x levels]
        n_prof, n_lev = store["Temp"].shape
        for name in ("Lat", "Lon", "Time", "SourceProf", "SourceFile"):
            if name not in store:  # (source indices only if asked for)
                continue
            store[name] = np.tile(store[name], (1, n_lev))
        store["Pres"] = np.tile(store["Pres"].T, (n_prof, 1))
        if "ScanUpFlag" in store:
            store["ScanUpFlag"] = np.tile(store["ScanUpFlag"], (1, n_lev))

        # geolocation correction [R2]
        scan_sign = np.where(scan_sign == 0, -1.0, scan_sign)  # +1 up, -1 down
        scan_sign = np.tile(scan_sign, (1, n_lev))
        m = 0.6  # km along track per km of vertical scan
        dz = store["Alt"] - 30
        dx = dz * scan_sign * m
        with np.errstate(invalid="ignore"):
            mean_alt = _nanmean_columns(store["Alt"])
        zidx = int(np.nanargmin(np.abs(mean_alt - 30)))
        lat, lon = store["Lat"], store["Lon"]
        if n_prof < 2:  # (no direction to move along: positions as in the file)
            record(inst_info)["notes"]["Lat"] = (  # (replaces "moved along track")
                "moved along track by 0.6 km per km from 30 km [R2], except in "
                "files with one profile (no direction to move along)")
            return conform_vars(store, settings)
        theta = azimuth(lat[:-1, zidx], lon[:-1, zidx], lat[1:, zidx], lon[1:, zidx])
        theta = np.tile(np.atleast_1d(theta)[:, None], (1, n_lev))
        theta = np.vstack([theta, theta[-1:]])
        store["Lat"], store["Lon"] = reckon(lat, lon, km2deg(dx), theta)

        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, "HIRDLS", inst_info)


def _get_hirdls(path, name):
    """One variable as float64, fill values (-999) as NaN; 1-D as a column.
    A name that is not a HIRDLS variable raises KeyError (get_HIRDLS prints a
    message and returns NaN, which then failed later with a misleading size
    error, or passed through as NaN with ignore_non_modal)."""
    if name in GEOLOCATION_FIELDS:
        group = "Geolocation Fields"
    elif name in DATA_FIELDS:
        group = "Data Fields"
    else:
        raise KeyError(f"Variable {name} not found in {path} (not a HIRDLS "
                       "variable: see GEOLOCATION_FIELDS and DATA_FIELDS in "
                       "atmosgw.curate.obs._limb_hirdls)")
    value = read_hdf5(path, f"{_SWATH}/{group}/{name}")
    value[value == -999] = np.nan
    return value[:, None] if value.ndim == 1 else value


def _nanmean_columns(x):
    """Mean of each column (level) over profiles, ignoring NaN."""
    finite = ~np.isnan(x)
    count = finite.sum(axis=0)
    total = np.where(finite, x, 0).sum(axis=0)
    return np.where(count > 0, total / np.maximum(count, 1), np.nan)
