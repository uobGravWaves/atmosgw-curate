"""
Trajectory reader for the DWD (Deutscher Wetterdienst) high-resolution
radiosondes ("Sekundenwerte aero", CDC open data, CC BY as noted beside
the data: credit DWD).

Files (as unzipped from
``opendata.dwd.de/climate_environment/CDC/observations_germany/radiosondes/
high_resolution/historical/<yyyy>/sekundenwerte_aero_<station>_..._hist.zip``),
anywhere below the folder:

- ``produkt_sec_aero_<yyyymmdd>_<yyyymmdd>_<station>.txt``: one per station
  and period (usually a year); ``;``-separated, one row per sample:
  ``STATIONS_ID; BEZUGSDATUM_SYNOP`` (reference time, yyyymmddhh);
  ``MESSZEITPUNKT`` (s); ``QN_1``; ``AE_GB_POS``, ``AE_GL_POS`` (degrees);
  ``AE_GPM`` (gpm); ``AE_P`` (hPa); ``AE_TT`` (degrees C); ``AE_RF`` (%);
  ``AE_TD`` (degrees C); ``AE_FF`` (m/s); ``AE_DD`` (degrees); ``eor``.
  -999 is missing. Rows are in order of reference time (every file checked
  is; a file found out of order is read in full instead).
- ``Metadaten_Sekunde_Aero_<station>.txt`` (beside it): per ascent
  ``Bezugsdatum``, ``AE_Start``, ``AE_Ende`` (yyyymmddHHMM), the burst
  height and ``Geraete_Name`` (sonde type).
- ``Metadaten_Geographie_<station>.txt``: station position, height and name
  over periods (``von_datum``, ``bis_datum``).

No field description is published beside the data (checked 2026-11-27: the
CDC folders hold only a station list), so the meanings below are inferred
from the files (2026-11-27, stations 01757 in 2004 and 2020 and 01303 in
1993):

- a flight's samples have one ``BEZUGSDATUM_SYNOP``; ``MESSZEITPUNKT`` is
  seconds since launch (its last value is within 50 s of AE_Ende - AE_Start
  for 95% of flights, the metadata times being to the minute); sample time
  = ``AE_Start`` + ``MESSZEITPUNKT``. Where AE_Start is not a valid time
  (43 of 383,195 metadata rows) or missing, the reference time is used
  (recorded). Times are UTC (not stated in the files; against ERA5 on four
  days of 2020, the wind differences are smallest, or within 0.01 m/s of
  the smallest, with no time shift, and grow for whole-hour shifts either
  way; an offset of less than an hour would be implausible: C. Wright,
  2026-11-27). Launches are 13-76 min
  before the reference time.
- ``AE_GPM`` is geopotential height: height increments equal the
  hypsometric thickness to 0.06% at 10-40 km (2020-01-07).
- ``AE_DD`` is the direction the wind blows from, clockwise from north:
  with u = -FF sin(DD), v = -FF cos(DD), the balloons' drift over 300
  flights correlates 0.98 (east) and 0.99 (north) with the mean wind.
- positions are the sonde's (they drift); early years give none (-999),
  e.g. 2004. There the station position (``Metadaten_Geographie``, for the
  flight's date) is used for every sample, and this is recorded.
- a position of exactly 0.00 is a second fill value (e.g. 7,245 samples of
  station 01757 in 2024 have latitude 0.00, most with longitude -999; the
  stations are all in Germany): set to NaN. Samples without a position
  inside a flight that has some are left NaN.

Output names: ``Lat``, ``Lon``, ``Time`` (datenum), ``Alt`` = AE_GPM / 1000
(km, geopotential), ``Pres`` = AE_P, ``Temp`` = AE_TT + 273.15, ``U``, ``V``
(m s-1, from AE_FF, AE_DD); extras by their file names (``AE_RF``,
``AE_TD``, ``AE_FF``, ``AE_DD``, ``QN_1``, ``MESSZEITPUNKT``). Samples are in
file order (the path order).

Origin: new code (no MATLAB equivalent).

Change log:
2026/11/27 [Claude] initial version
2026/11/27 [Claude] UTC and geopotential heights confirmed (ERA5, hypsometric checks)
"""

import math
import os
import re

import numpy as np
import pandas as pd

from atmosgw.curate.matlab_time import datenum, datenum_to_datetime
from atmosgw.curate.obs._common import find_in, note, note_sources, record

KELVIN = 273.15
MISSING = -999.0
COLUMNS = ("STATIONS_ID", "BEZUGSDATUM_SYNOP", "MESSZEITPUNKT", "QN_1", "AE_GB_POS",
           "AE_GL_POS", "AE_GPM", "AE_P", "AE_TT", "AE_RF", "AE_TD", "AE_FF", "AE_DD")
EXTRAS = ("AE_RF", "AE_TD", "AE_FF", "AE_DD", "QN_1", "MESSZEITPUNKT")
_PRODUCT = re.compile(r"produkt_sec_aero_(\d{8})_(\d{8})_(\d{5})\.txt$")
_ENCODING = "latin-1"  # (the metadata files use ISO 8859-1)


def load(settings, inst_info):
    """
    The flights launched in settings["time_range"] (datenums) at the
    stations asked for. Returns (flights, file_list): a list of dicts, one
    per flight, with per-sample arrays (Lat, Lon, Time, Alt, Pres, Temp, U,
    V and the extras asked for) and per-flight values (flight_id,
    platform_id, platform_name, platform_type, start_time, reference_time,
    source_file); file_list: the product files read (names).
    """
    note_sources(inst_info, {"Lat": "AE_GB_POS", "Lon": "AE_GL_POS",
                             "Time": "AE_Start (metadata) + MESSZEITPUNKT",
                             "Alt": "AE_GPM", "Pres": "AE_P", "Temp": "AE_TT",
                             "U": "AE_FF, AE_DD", "V": "AE_FF, AE_DD"})
    note(inst_info, "Alt", "geopotential height (AE_GPM, gpm / 1000)")
    note(inst_info, "Temp", "AE_TT (degrees C) + 273.15")
    note(inst_info, "U", "-AE_FF sin(AE_DD): AE_DD is the direction the wind "
                         "blows from (inferred from the balloons' drift)")
    note(inst_info, "V", "-AE_FF cos(AE_DD)")
    start, end = settings["time_range"]
    # reference times read: a day either side (launch is before the reference)
    ref_lo = _ref_number(math.floor(start) - 1)
    ref_hi = _ref_number(math.floor(end) + 2) - 1
    first_day, last_day = (int(_ref_number(math.floor(t)) // 100)
                           for t in (start - 1, end + 1))
    wanted = _station_filter(settings.get("stations"))
    by_name = {}
    for path in find_in(inst_info, "produkt_sec_aero_*.txt"):
        found = _PRODUCT.search(os.path.basename(path))
        if not found or not os.path.isfile(path):
            continue
        begin, finish = int(found.group(1)), int(found.group(2))
        if finish < first_day or begin > last_day:
            continue
        by_name.setdefault(os.path.basename(path), []).append(
            (found.group(3), begin, path))
    products = []
    for file_name, copies in sorted(by_name.items()):
        # the same file in several folders (e.g. 1995/ holds older, shorter
        # copies of 1993 files): the copy in the folder of its own year, else
        # the largest; recorded
        def rank(copy):
            in_year = os.path.basename(os.path.dirname(copy[2])) == str(copy[1])[:4]
            return (not in_year, -os.path.getsize(copy[2]))
        chosen = sorted(copies, key=rank)[0]
        if len(copies) > 1:
            record(inst_info)["multiple_matches"].append(
                f"{file_name}: " + ", ".join(os.path.relpath(c[2], inst_info["path"])
                                             for c in copies)
                + f" -> used {os.path.relpath(chosen[2], inst_info['path'])}")
        products.append(chosen)
    products.sort()
    seen_flights = set()  # (station, reference time): each ascent once
    flights, file_list = [], []
    counts = {"filled_positions": 0, "no_launch_time": 0, "flights": 0}
    names_seen = {}
    for station, _, path in products:
        folder = os.path.dirname(path)
        geography = _read_geography(os.path.join(
            folder, f"Metadaten_Geographie_{station}.txt"))
        name = geography[0][3] if geography else ""
        names_seen[station] = name
        if wanted is not None and station not in wanted and name.lower() not in wanted:
            continue
        if settings["verbose"]:
            print(f"----> Loading DWD radiosondes: {station} {name} "
                  f"({os.path.basename(path)})")
        rows = read_rows(path, ref_lo, ref_hi)
        if rows.shape[0] == 0:
            continue
        launches = _read_launches(os.path.join(
            folder, f"Metadaten_Sekunde_Aero_{station}.txt"))
        refs = rows[:, 1].astype(np.int64)
        breaks = np.flatnonzero(np.diff(refs) != 0) + 1
        file_used = False
        for chunk in np.split(np.arange(refs.size), breaks):
            ref = int(refs[chunk[0]])
            if (station, ref) in seen_flights:  # (in two overlapping files)
                counts["repeated"] = counts.get("repeated", 0) + 1
                continue
            seen_flights.add((station, ref))
            launch, sonde = launches.get(ref, (None, ""))
            ref_time = _ref_datenum(ref)
            if launch is None:
                counts["no_launch_time"] += 1
                launch = ref_time
            if not start <= launch <= end:
                continue
            if not file_used:
                file_list.append(os.path.basename(path))
                file_used = True
            flight = _flight(rows[chunk], launch, settings)
            if np.all(np.isnan(flight["Lat"])) or np.all(np.isnan(flight["Lon"])):
                lat, lon = _station_position(geography, ref // 100)
                missing = np.isnan(flight["Lat"]) | np.isnan(flight["Lon"])
                flight["Lat"] = np.where(missing, lat, flight["Lat"])
                flight["Lon"] = np.where(missing, lon, flight["Lon"])
                counts["filled_positions"] += 1
            flight.update(
                flight_id=f"{station}-{str(datenum_to_datetime(launch))[:16]}",
                platform_id=station, platform_name=name, platform_type=sonde,
                start_time=launch, reference_time=ref_time,
                source_file=len(file_list))
            flights.append(flight)
            counts["flights"] += 1
    if wanted is not None:
        unknown = sorted(w for w in wanted if w not in names_seen
                         and w not in {n.lower() for n in names_seen.values()})
        if unknown and products:
            record(inst_info)["files_not_found"].extend(f"station {w}" for w in unknown)
    if not products:
        record(inst_info)["files_not_found"].append(
            f"produkt_sec_aero files for {str(datenum_to_datetime(start))[:10]} to "
            f"{str(datenum_to_datetime(end))[:10]}")
    position = (f"the sonde's position; the station's where the file gives none "
                f"({counts['filled_positions']} of {counts['flights']} flights)")
    note(inst_info, "Lat", position)
    note(inst_info, "Lon", position)
    note(inst_info, "Time", f"launch time (AE_Start) + MESSZEITPUNKT (s), UTC "
                            f"(checked against ERA5); reference time used as the "
                            f"launch time for {counts['no_launch_time']} flights "
                            f"without a valid AE_Start")
    if counts.get("repeated"):
        record(inst_info)["multiple_matches"].append(
            f"{counts['repeated']} ascent(s) in more than one file: read once")
    return flights, file_list


def read_rows(path, ref_lo, ref_hi):
    """
    The sample rows of one product file whose reference time (yyyymmddhh)
    is in [ref_lo, ref_hi], as a float array [rows x 13] (COLUMNS), -999 as
    NaN. The file is sorted by reference time, so the first row is found by
    bisection on byte offsets; if the rows read are out of order, the whole
    file is read instead.
    """
    with open(path, "rb") as handle:
        header = handle.readline()
        body = len(header)
        size = os.fstat(handle.fileno()).st_size
        lo, hi = body, size
        while hi - lo > 65536:  # (then scan forward from lo)
            mid = (lo + hi) // 2
            handle.seek(mid)
            handle.readline()  # (the rest of a partial line)
            line = handle.readline()
            if not line or _ref_of(line) >= ref_lo:
                hi = mid
            else:
                lo = mid
        handle.seek(lo)
        if lo > body:
            handle.readline()
        lines = []
        for line in handle:
            ref = _ref_of(line)
            if ref < ref_lo:
                continue
            if ref > ref_hi:
                break
            lines.append(line)
    rows = _parse(lines)
    if rows.shape[0] and np.any(np.diff(rows[:, 1]) < 0):  # (not sorted)
        with open(path, "rb") as handle:
            handle.readline()
            lines = [line for line in handle if ref_lo <= _ref_of(line) <= ref_hi]
        rows = _parse(lines)
        order = np.argsort(rows[:, 1], kind="stable")
        rows = rows[order]
    return rows


def _flight(rows, launch, settings):
    """Per-sample arrays of one flight (rows of read_rows)."""
    col = {name: rows[:, k] for k, name in enumerate(COLUMNS)}
    speed, direction = col["AE_FF"], np.deg2rad(col["AE_DD"])
    # (exactly 0 is a second fill value for positions: German stations)
    lat = np.where(col["AE_GB_POS"] == 0, np.nan, col["AE_GB_POS"])
    lon = np.where(col["AE_GL_POS"] == 0, np.nan, col["AE_GL_POS"])
    out = {
        "Lat": lat, "Lon": lon,
        "Time": launch + col["MESSZEITPUNKT"] / 86400.0,
        "Alt": col["AE_GPM"] / 1000.0, "Pres": col["AE_P"],
        "Temp": col["AE_TT"] + KELVIN,
        "U": -speed * np.sin(direction), "V": -speed * np.cos(direction),
    }
    for name in settings["additional_vars"]:
        if name not in col:
            raise KeyError(f"Variable {name} not found in the DWD radiosonde files "
                           f"(available: {', '.join(EXTRAS)})")
        out[name] = col[name]
    return out


def _parse(lines):
    """Rows of ``;``-separated values (first 13 fields) as float, -999 as NaN."""
    if not lines:
        return np.empty((0, len(COLUMNS)))
    text = b"".join(lines).decode(_ENCODING)
    values = np.array([line.split(";")[:len(COLUMNS)] for line in text.splitlines()],
                      dtype=float)
    values[values == MISSING] = np.nan
    return values


def _ref_of(line):
    """The reference time (yyyymmddhh) of one product line, as an integer."""
    return int(line.split(b";", 2)[1])


def _ref_number(day):
    """yyyymmddhh (hour 00) of an integer datenum, as an integer."""
    return int(str(datenum_to_datetime(day))[:10].replace("-", "")) * 100


def _ref_datenum(ref):
    """Datenum of a reference time yyyymmddhh."""
    text = str(ref)
    return datenum(int(text[:4]), int(text[4:6]), int(text[6:8]), int(text[8:10]))


def _read_launches(path):
    """{reference time: (launch datenum or None, sonde type)} from a
    Metadaten_Sekunde_Aero file ({} if there is none)."""
    if not os.path.isfile(path):
        return {}
    table = pd.read_csv(path, sep=";", encoding=_ENCODING, dtype=str,
                        index_col=False, skipinitialspace=True)
    table.columns = [c.strip() for c in table.columns]
    out = {}
    for ref, start, sonde in zip(table["Bezugsdatum"], table["AE_Start"],
                                 table.get("Geraete_Name", [""] * len(table))):
        ref, start = str(ref).strip(), str(start).strip()
        if not ref.isdigit():
            continue
        launch = None
        if len(start) == 12 and start.isdigit():
            launch = datenum(int(start[:4]), int(start[4:6]), int(start[6:8]),
                             int(start[8:10]), int(start[10:12]))
        out[int(ref)] = (launch, "" if pd.isna(sonde) else str(sonde).strip())
    return out


def _read_geography(path):
    """[(from yyyymmdd, to yyyymmdd or 99999999, (lat, lon), name)] from a
    Metadaten_Geographie file ([] if there is none)."""
    if not os.path.isfile(path):
        return []
    table = pd.read_csv(path, sep=";", encoding=_ENCODING, dtype=str,
                        index_col=False, skipinitialspace=True)
    table.columns = [c.strip() for c in table.columns]
    out = []
    for _, row in table.iterrows():
        begin = str(row.get("von_datum", "")).strip()
        finish = str(row.get("bis_datum", "")).strip()
        out.append((int(begin) if begin.isdigit() else 0,
                    int(finish) if finish.isdigit() else 99999999,
                    (float(row["Geogr.Breite"]), float(row["Geogr.Laenge"])),
                    str(row.get("Stationsname", "")).strip()))
    return out


def _station_position(geography, day):
    """(lat, lon) of the station on day (yyyymmdd); NaN if not known."""
    for begin, finish, position, _ in geography:
        if begin <= day <= finish:
            return position
    return (np.nan, np.nan)


def _station_filter(stations):
    """Set of station ids (5 digits) and lower-case names, or None for all."""
    if stations is None:
        return None
    if isinstance(stations, (str, int)):
        stations = [stations]
    out = set()
    for station in stations:
        text = str(station).strip()
        out.add(f"{int(text):05d}" if text.isdigit() else text.lower())
    return out
