"""
Trajectory reader for IAGOS (In-service Aircraft for a Global Observing
System, including the MOZAIC era) Level 2 time series: one netCDF file per
flight, 1994-2024.

Files: ``IAGOS_timeseries_<yyyymmddHHMMSS><nn>[_L2_<version>].nc[4]``
(departure time and a two-digit number), in ``TimeSeries/iagos_<yyyy>/``.
Unzipped files anywhere below the folder are preferred; zip archives (on
eepc-0184 ``IAGOS/iagos.zip``, the version used in Wright et al. (2024,
ACP; UNVERIFIED: reference not checked here)) are the backup for flights
not unzipped, read in memory. Copies are matched by flight (departure time
and number), so any unzipped copy (e.g. a newer version) replaces the
archive's.

Two generations (checked on four flights, 1998, 2010, 2019, 2023):

- old (``.nc``, NETCDF3, IAGOS-MOZAIC): validity variables ``<var>_val``,
  all -9999 in the files checked (no screening possible); no GPS altitude
  (all -9999); package temperature ``air_temp_PM``; aircraft in ``title``
  (e.g. "Instrument-75, MOZAIC-05, A340, OE-LAG, Austrian, ...");
- new (``.nc4``, L2 3.1.x): ``<var>_validity_flag`` with flag values 0
  good, 2 limited, 3 erroneous, 4 not_validated, 7 missing_values;
  package temperature ``air_temp_P1`` (not on every flight); aircraft in
  ``platform`` (e.g. "Aircraft, IAGOS-01, MSN-304, D-AIGT, Lufthansa").

Both: ``UTC_time`` (s since the departure day, from ``units``), ``lat``,
``lon``, ``baro_alt_AC``, ``gps_alt_AC`` (m), ``air_press_AC`` (Pa),
``air_temp_AC`` (K, the aircraft's sensor), ``zon_wind_AC``,
``mer_wind_AC`` (m s-1), ozone, CO, water vapour and others, every 4 s;
-9999 (and 0 for new files' time) as fill; ``departure_airport``,
``arrival_airport``, ``flight_name`` attributes.

Output (decided by C. Wright, 2026-11-27): ``Temp`` = the IAGOS package
temperature (P1 or PM) where the flight has values not flagged 3 or 7 (3
of 150 L2 3.1 flights reviewed had none), else the aircraft's
(``air_temp_AC``), per flight in ``temperature_source``; ``Alt`` = GPS
altitude (geometric) likewise, else barometric (pressure altitude), per
flight in ``altitude_source``; ``Pres`` =
``air_press_AC`` / 100 (hPa); ``U``, ``V`` = ``zon_wind_AC``,
``mer_wind_AC``. Screening (``<var>Bad``): validity flag 3 (erroneous) or 7
(missing) is ``instrument_flag``; 0, 2 and 4 are kept (old files: no
screening). Every other file variable is available through
``additional_vars`` by its name (NaN where a flight lacks it). Pressures
<= 0 are missing too (a 2023-11-15 flight of D-AIKE has 79 values of
-999900 Pa); its 82 take-off samples of ``air_temp_AC`` exactly 173.16 K
(flag 0) look like a disguised fill value but are kept. ``platform_id`` is
the registration, the fourth field of ``platform`` (new) or ``title``
(old).

Flights are selected by their departure time (``departure_UTC_time``, else
the first time in the file); only files whose name time is within a day of
the request are opened.

Origin: new code (no MATLAB equivalent).

Change log:
2026/11/27 [Claude] initial version
2026/11/27 [Claude] review fixes: sources from usable values; copies only among the
  flights asked for; -9999.9 fill
2026/11/27 [Claude] unzipped files preferred, the zip as backup, silently;
  copies matched by flight, not file name; newest version among copies
"""

import math
import os
import re
import zipfile

import netCDF4
import numpy as np

from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs._common import note, note_sources, note_version, record
from atmosgw.curate.readers import list_tree

_NAME = re.compile(r"IAGOS_timeseries_(\d{14})(\d{2})(?:_L2_([\d.]+))?\.nc4?$")
_UNITS = re.compile(r"seconds since (\d{4})-(\d\d)-(\d\d)[ T](\d\d):(\d\d):(\d\d)")
_ISO = re.compile(r"(\d{4})-(\d\d)-(\d\d)T(\d\d):(\d\d):(\d\d)")
REMOVED_FLAGS = (3, 7)  # erroneous, missing_values (new files)
FILL = -9999.0
PACKAGE_TEMPS = ("air_temp_P1", "air_temp_PM")


def load(settings, inst_info):
    """Flights departing in settings["time_range"] (datenums). Returns
    (flights, file_list); see trajectory.py for the flight dicts."""
    note_sources(inst_info, {"Lat": "lat", "Lon": "lon", "Time": "UTC_time",
                             "Alt": "gps_alt_AC, else baro_alt_AC",
                             "Pres": "air_press_AC", "Temp": "air_temp_P1/PM, else "
                             "air_temp_AC", "U": "zon_wind_AC", "V": "mer_wind_AC"})
    note(inst_info, "Pres", "air_press_AC (Pa) / 100")
    note(inst_info, "Temp", "IAGOS package temperature (air_temp_P1 or air_temp_PM) "
                            "where the flight has values not flagged 3 or 7, else "
                            "the aircraft's (air_temp_AC); per flight in "
                            "temperature_source")
    note(inst_info, "Alt", "GPS altitude (geometric, gps_alt_AC) where the flight "
                           "has values not flagged 3 or 7, else barometric "
                           "altitude (baro_alt_AC, pressure altitude); per "
                           "flight in altitude_source")
    for key in ("Pres", "Temp", "U", "V"):
        note(inst_info, f"{key}Bad", "validity flag 3 (erroneous) or 7 (missing)")
    start, end = settings["time_range"]
    sources = _sources(inst_info["path"])
    # the files whose name time is within a day of the request; one copy of
    # each flight (same departure time and number, whatever the version):
    # unzipped files are preferred and the zip is the backup for flights not
    # unzipped (silently: that is how the archive is meant to be used); among
    # unzipped copies, or among zip copies (29 flights of January 2023 are in
    # two year folders), the newest version, else the first in search order,
    # recorded (only for the flights looked at)
    copies = {}
    for order, (source, member, stamp, flight, version) in enumerate(sources):
        when = datenum(int(stamp[:4]), int(stamp[4:6]), int(stamp[6:8]),
                       int(stamp[8:10]), int(stamp[10:12]), int(stamp[12:14]))
        if math.floor(start) - 1 <= when <= end + 1:
            copies.setdefault(flight, []).append(
                (source.endswith(".zip"), tuple(-v for v in version), order,
                 when, source, member))
    candidates = []
    for flight, entries in copies.items():
        entries.sort()
        candidates.append(entries[0][3:])
        same_kind = [e for e in entries if e[0] == entries[0][0]]
        if len(same_kind) > 1:
            record(inst_info)["multiple_matches"].append(
                f"flight {flight}: {len(same_kind)} copies -> used "
                + "/".join(entries[0][5].replace(os.sep, "/").split("/")[-2:]))
    candidates.sort()
    if not sources:
        record(inst_info)["files_not_found"].append(
            "IAGOS_timeseries files (or a zip of them)")
    flights, file_list = [], []
    archives = {}
    try:
        for k, (_, source, member) in enumerate(candidates, start=1):
            if settings["verbose"]:
                print(f"----> Loading IAGOS flight {k} of {len(candidates)}")
            name = os.path.basename(member)
            if source.endswith(".zip"):
                if source not in archives:
                    archives[source] = zipfile.ZipFile(source)
                handle = netCDF4.Dataset(name, memory=archives[source].read(member))
            else:
                handle = netCDF4.Dataset(member)
            with handle:
                flight = read_flight(handle, settings)
            if flight is None or not start <= flight["start_time"] <= end:
                continue
            file_list.append(name)
            note_version(inst_info, name, flight.pop("version"))
            flight["source_file"] = len(file_list)
            flights.append(flight)
    finally:
        for archive in archives.values():
            archive.close()
    return flights, file_list


def read_flight(handle, settings):
    """One flight's per-sample arrays and per-flight values from an open
    netCDF file (None if it has no times)."""
    handle.set_auto_maskandscale(False)
    names = set(handle.variables)

    def get(name):
        if name not in names:
            return None
        value = np.asarray(handle[name][:], dtype=float)
        fill = getattr(handle[name], "_FillValue", FILL)
        # (-9999.9 also appears, e.g. in baro_alt_AC: anything <= -9999 is fill)
        value[(value == float(fill)) | (value <= FILL)] = np.nan
        return value

    def bad(name):
        """True where the variable's validity flag is 3 or 7 (new files)."""
        flag = handle[f"{name}_validity_flag"][:] if f"{name}_validity_flag" in names \
            else None
        if flag is None:
            return np.zeros(n, dtype=bool)
        return np.isin(np.asarray(flag), REMOVED_FLAGS)

    seconds = np.asarray(handle["UTC_time"][:], dtype=float)
    found = _UNITS.search(str(getattr(handle["UTC_time"], "units", "")))
    if not found or seconds.size == 0:
        return None
    base = datenum(*(int(g) for g in found.groups()))
    keep = seconds > 0 if str(getattr(handle["UTC_time"], "_FillValue", "")) == "0.0" \
        else np.isfinite(seconds)
    n = seconds.size
    time = np.where(keep, base + seconds / 86400.0, np.nan)

    def usable(name):
        """True if the file has values of name not flagged 3 or 7."""
        value = get(name)
        return value is not None and bool(np.any(np.isfinite(value) & ~bad(name)))

    # (a package temperature or GPS altitude whose every value is flagged
    # erroneous or missing is not used: about 2% of L2 3.1 flights)
    temp_source = next((t for t in PACKAGE_TEMPS if usable(t)), "air_temp_AC")
    alt_source = "gps_alt_AC" if usable("gps_alt_AC") else "baro_alt_AC"
    out = {
        "Lat": get("lat"), "Lon": get("lon"), "Time": time,
        "Alt": get(alt_source) / 1000.0,
        # (pressure <= 0 is also missing: one 2023 flight has -999900 Pa, the
        # -9999 fill in hPa scaled to Pa)
        "Pres": np.where(get("air_press_AC") > 0, get("air_press_AC"), np.nan) / 100.0,
        "Temp": get(temp_source),
        "U": get("zon_wind_AC"), "V": get("mer_wind_AC"),
        "PresBad": bad("air_press_AC"), "TempBad": bad(temp_source),
        "UBad": bad("zon_wind_AC"), "VBad": bad("mer_wind_AC"),
    }
    for key, value in list(out.items()):
        if value is None:  # (a variable the file does not have)
            out[key] = np.full(n, np.nan)
    for extra in settings["additional_vars"]:
        value = get(extra)
        out[extra] = np.full(n, np.nan) if value is None else value
    attrs = {a: str(handle.getncattr(a)) for a in handle.ncattrs()}
    departure = _ISO.search(attrs.get("departure_UTC_time", ""))
    start = (datenum(*(int(g) for g in departure.groups())) if departure
             else np.nanmin(time) if np.any(np.isfinite(time)) else np.nan)
    platform = attrs.get("platform") or attrs.get("title", "")
    parts = [p.strip() for p in platform.split(",")]
    registration = parts[3] if len(parts) > 3 else ""  # (e.g. D-AIGT, N384HA)
    out.update(
        flight_id=attrs.get("flight_name", ""),
        platform_id=registration, platform_name=platform,
        platform_type=("aircraft (IAGOS-MOZAIC)" if "MOZAIC" in attrs.get("source", "")
                       else "aircraft (IAGOS)"),
        departure_airport=attrs.get("departure_airport", ""),
        arrival_airport=attrs.get("arrival_airport", ""),
        temperature_source=temp_source, altitude_source=alt_source,
        start_time=start,
        version=(f"L2 {attrs['file_version']}" if "file_version" in attrs
                 else f"MOZAIC L2 revision {attrs.get('file_revision', '')}".strip()),
    )
    return out


def _sources(folder):
    """[(source, member or path, yyyymmddHHMMSS, flight, version)] of every
    IAGOS file below folder: members of zip archives, and files. flight is
    the departure time and number (the file's ``flight_name``); version the
    L2 version as a tuple of integers (() if the name has none)."""
    out = []

    def add(source, member):
        found = _NAME.search(member)
        if found:
            version = tuple(int(v) for v in found.group(3).split(".") if v.isdigit()) \
                if found.group(3) else ()
            out.append((source, member, found.group(1),
                        found.group(1) + found.group(2), version))

    for path in list_tree(folder):
        if path.endswith(".zip") and os.path.isfile(path):
            with zipfile.ZipFile(path) as archive:
                for member in archive.namelist():
                    add(path, member)
        elif _NAME.search(path) and os.path.isfile(path):
            add(path, path)
    return out
