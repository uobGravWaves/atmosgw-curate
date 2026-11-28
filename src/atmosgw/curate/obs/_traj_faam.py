"""
Trajectory reader for FAAM (Facility for Airborne Atmospheric Measurements,
BAe-146-301 research aircraft) core data files: one netCDF file per flight.

Files: ``core_faam_<yyyymmdd>_v<nnn>_r<rev>_<flight>.nc`` (e.g.
``core_faam_20250709_v005_r0_c420.nc``), anywhere below the FAAM folder.
Variables (checked on the 27 files of 2025-07-09 to 2025-08-01 on
eepc-0184, v005 r0): ``Time`` (s since the flight day, one per second);
measurements as (Time, sps<N>): N samples per second (32 for the ones read
here), with ``<name>_FLAG`` (bit masks; 0 is good; e.g. ``TAT_DI_R_FLAG``
1 mach_out_of_range, 2 discrepancy_threshold_exceeded, 4 in_cloud) and
-9999 as fill. Read: ``LAT_GIN``, ``LON_GIN``, ``ALT_GIN`` (GPS-aided
inertial navigation: geometric altitude, m), ``PS_RVSM`` (hPa),
``TAT_DI_R`` (K, the deiced Rosemount "true air temperature", which in
FAAM's naming is the ambient air temperature), ``U_C``, ``V_C``, ``W_C``
(m s-1, from the turbulence probe and GIN). Attributes ``flight_number``,
``platform``, ``time_coverage_start``, ``processing_software_version``,
``revision_number``.

Choices: ``TAT_DI_R`` rather than the non-deiced ``TAT_ND_R`` (available
as an extra), as it is usable in icing too (UNVERIFIED which FAAM
recommends for general use); ``Alt`` = ``ALT_GIN`` (GPS, decided by C.
Wright, 2026-11-27), pressure altitude ``PALT_RVS`` as an extra. A sample
whose flag is non-zero is ``instrument_flag`` (``<key>Bad``; the rule is in
each variable's ``comment``).

Sample rate (``sample_rate``, decided by C. Wright, 2026-11-27): "1hz"
(default): for each second, the mean of that second's unflagged samples
(longitude: the mean on the circle, so a second across the dateline stays
there; a position sample outside -90..90 or -180..360 gives the second
that value, so it is bad_geolocation, not averaged into a plausible one;
NaN, and the flag kept, if none is unflagged), at the centre of the
second's samples (Time + 15.5/32 s for 32 Hz); "full": every sample, at
Time + k/N s. Extras at other rates are averaged the same way ("1hz") or
repeated to 32 Hz ("full"; rates that do not divide 32, e.g. 20 or 64 Hz,
are refused there). Flagged positions (``LAT_GIN_FLAG`` ...) have no value.
Of several copies of a flight in the folder, the highest core version
(``_v<nnn>_``), then the highest revision, is read.

Flights are selected by their start (``time_coverage_start``, else the
first time).

Origin: new code (no MATLAB equivalent).

Change log:
2026/11/27 [Claude] initial version
2026/11/27 [Claude] review fixes: highest revision per flight; flagged positions NaN;
  extras of rates other than divisors of 32 refused at sample_rate='full'
2026/11/27 [Claude] flag rule in the variables' comments
2026/11/27 [Claude] review: highest core version before revision; scale/offset applied
2026/11/28 [Claude] review: 1 Hz longitude averaged on the circle; out-of-range
  position samples not averaged
"""

import math
import os
import re
import warnings

import netCDF4
import numpy as np

from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs._common import find_in, note, note_sources, note_version, record

_NAME = re.compile(r"core_faam_(\d{8})_v(\d+)_r(\d+)_(\w+?)(?:_1hz)?\.nc$")
_UNITS = re.compile(r"seconds since (\d{4})-(\d\d)-(\d\d)[ T](\d\d):(\d\d):(\d\d)")
_ISO = re.compile(r"(\d{4})-(\d\d)-(\d\d)T(\d\d):(\d\d):(\d\d)")
FULL_RATE = 32
POSITION_RANGE = {"Lat": (-90.0, 90.0), "Lon": (-180.0, 360.0)}  # (file values)
FILE_NAMES = {"Lat": "LAT_GIN", "Lon": "LON_GIN", "Alt": "ALT_GIN", "Pres": "PS_RVSM",
              "Temp": "TAT_DI_R", "U": "U_C", "V": "V_C", "W": "W_C"}


def load(settings, inst_info):
    """Flights starting in settings["time_range"] (datenums). Returns
    (flights, file_list); see trajectory.py for the flight dicts."""
    note_sources(inst_info, dict(FILE_NAMES, Time="Time (+ sample within second)"))
    note(inst_info, "Alt", "ALT_GIN (GPS-aided inertial navigation, geometric "
                           "altitude) / 1000; pressure altitude PALT_RVS as an extra")
    note(inst_info, "Temp", "TAT_DI_R (deiced; the non-deiced TAT_ND_R as an extra)")
    rate = settings.get("sample_rate", "1hz")
    for key in ("Pres", "Temp", "U", "V", "W"):
        note(inst_info, f"{key}Bad",
             f"{FILE_NAMES[key]}_FLAG non-zero" + (" for every sample of the second"
                                                   if rate == "1hz" else ""))
    note(inst_info, "Time", "1 Hz: mean of each second's unflagged 32 Hz samples, at "
                            "their centre" if rate == "1hz" else
                            "every 32 Hz sample, at Time + k/32 s")
    start, end = settings["time_range"]
    files = []
    for path in find_in(inst_info, "core_faam_*.nc"):
        found = _NAME.search(os.path.basename(path))
        if not found or not os.path.isfile(path) or path.endswith("_1hz.nc"):
            continue
        day = found.group(1)
        when = datenum(int(day[:4]), int(day[4:6]), int(day[6:8]))
        if math.floor(start) - 1 <= when <= end:
            files.append((found.group(4), int(found.group(2)),
                          int(found.group(3)), path))
    # one copy of each flight: the highest core version, then the highest
    # revision (recorded if there are several)
    newest = {}
    for flight, version, revision, path in sorted(files):
        if flight in newest:  # (files named by their path below the folder)
            earlier, later = (os.path.relpath(p, inst_info["path"])
                              for p in (newest[flight][2], path))
            record(inst_info)["multiple_matches"].append(
                f"flight {flight}: {earlier} and {later} -> used {later}")
        newest[flight] = (version, revision, path)
    files = [path for _, _, path in newest.values()]
    if not find_in(inst_info, "core_faam_*.nc"):
        record(inst_info)["files_not_found"].append("core_faam files")
    flights, file_list = [], []
    for path in sorted(files):
        with netCDF4.Dataset(path) as handle:
            flight = read_flight(handle, settings)
        if flight is None or not start <= flight["start_time"] <= end:
            continue
        file_list.append(os.path.basename(path))
        note_version(inst_info, path, flight.pop("version"))
        flight["source_file"] = len(file_list)
        flights.append(flight)
    return flights, file_list


def _unpacked(var, value):
    """value with the variable's scale_factor and add_offset applied, if it
    has any (automatic scaling is off, so that fill values are seen as
    stored; none of the 2025 core variables is packed)."""
    scale = float(getattr(var, "scale_factor", 1.0))
    offset = float(getattr(var, "add_offset", 0.0))
    return value if (scale == 1.0 and offset == 0.0) else value * scale + offset


def read_flight(handle, settings):
    """One flight's per-sample arrays and per-flight values (None if the
    file has no times)."""
    handle.set_auto_maskandscale(False)
    rate = settings.get("sample_rate", "1hz")
    seconds = np.asarray(handle["Time"][:], dtype=float)
    found = _UNITS.search(str(getattr(handle["Time"], "units", "")))
    if not found or seconds.size == 0:
        return None
    base = datenum(*(int(g) for g in found.groups()))
    seconds[seconds == -9999] = np.nan
    seconds = _unpacked(handle["Time"], seconds)

    def per_sample(name, circular=False, valid=None):
        """(values, flagged) of a variable at the output rate; circular: a
        longitude (degrees), averaged on the circle where a second spans
        more than 180 degrees; valid: (low, high) of a position: at 1 Hz, a
        second with an unflagged sample outside it takes that sample's value
        (so it is bad_geolocation, as at full rate) instead of a mean."""
        var = handle[name]
        value = np.asarray(var[:], dtype=float)
        value[value == float(getattr(var, "_FillValue", -9999.0))] = np.nan
        value = _unpacked(var, value)
        flag_name = f"{name}_FLAG"
        if flag_name in handle.variables:
            flagged = np.asarray(handle[flag_name][:]) != 0
        else:
            flagged = np.zeros(value.shape, dtype=bool)
        if value.ndim == 1:
            value, flagged = value[:, None], flagged.reshape(-1, 1)
        if rate == "1hz":
            good = np.where(flagged | np.isnan(value), np.nan, value)
            with warnings.catch_warnings():  # (seconds with no good sample: NaN)
                warnings.simplefilter("ignore", RuntimeWarning)
                mean = np.nanmean(good, axis=1)
                if circular:  # (only seconds spanning > 180 degrees: the dateline)
                    across = (np.nanmax(good, axis=1) - np.nanmin(good, axis=1)) > 180
                    rad = np.deg2rad(good[across])
                    mean[across] = np.rad2deg(np.arctan2(
                        np.nanmean(np.sin(rad), axis=1),
                        np.nanmean(np.cos(rad), axis=1)))
            if valid is not None:
                with np.errstate(invalid="ignore"):
                    out = (good < valid[0]) | (good > valid[1])
                rows = np.flatnonzero(out.any(axis=1))
                mean[rows] = good[rows, np.argmax(out[rows], axis=1)]
            no_good = np.all(flagged | np.isnan(value), axis=1)
            return mean, no_good & np.any(flagged, axis=1)
        if FULL_RATE % value.shape[1]:  # (e.g. 20 or 64 Hz: no exact placing)
            raise ValueError(
                f"{name} is at {value.shape[1]} Hz, which cannot be placed on the "
                f"{FULL_RATE} Hz samples of sample_rate='full'; use "
                "sample_rate='1hz' for it.")
        repeat = FULL_RATE // value.shape[1]
        return (np.repeat(value, repeat, axis=1).ravel(),
                np.repeat(flagged, repeat, axis=1).ravel())

    if rate == "1hz":
        time = base + (seconds + (FULL_RATE - 1) / (2 * FULL_RATE)) / 86400.0
    else:
        time = (base + (seconds[:, None] + np.arange(FULL_RATE)[None, :] / FULL_RATE)
                / 86400.0).ravel()
    out = {"Time": time}
    for key, name in FILE_NAMES.items():
        value, flagged = per_sample(name, circular=key == "Lon",
                                    valid=POSITION_RANGE.get(key))
        if key == "Alt":
            value = value / 1000.0
        if key in ("Lat", "Lon", "Alt"):
            value = np.where(flagged, np.nan, value)  # (flagged position: none)
        out[key] = value
        if key in ("Pres", "Temp", "U", "V", "W"):
            out[f"{key}Bad"] = flagged
    for extra in settings["additional_vars"]:
        if extra not in handle.variables:
            raise KeyError(f"Variable {extra} not found in {handle.filepath()}")
        out[extra] = per_sample(extra)[0]
    attrs = {a: str(handle.getncattr(a)) for a in handle.ncattrs()}
    began = _ISO.search(attrs.get("time_coverage_start", ""))
    start = (datenum(*(int(g) for g in began.groups())) if began
             else np.nanmin(time))
    out.update(
        flight_id=attrs.get("flight_number", ""),
        platform_id="FAAM BAe-146-301" if "146" in attrs.get("platform", "")
        else attrs.get("platform", ""),
        platform_name=attrs.get("platform", ""), platform_type="research aircraft",
        start_time=start,
        version=(f"core v{_NAME.search(os.path.basename(handle.filepath())).group(2)} "
                 f"revision {attrs.get('revision_number', '')}, processing "
                 f"{attrs.get('processing_software_version', '')}"),
    )
    return out

