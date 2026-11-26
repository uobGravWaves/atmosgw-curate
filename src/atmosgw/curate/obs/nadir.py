"""
Load nadir-sounder swaths (AIRS, IASI; later CrIS, IASI-NG) in a common format.

Load swaths from cross-track scanning nadir sounders (retrieved temperature,
or brightness temperature for 2-D products) and return them in the common
**swath** format (``docs/design/data_formats.md``), on a regular grid in
distance by default. Port of the MATLAB ``prep_airs_3d`` (corwin365/
MatlabFunctions, commit f1a78c0), without its detrending and smoothing,
which are analysis steps (``docs/design/nadir_loader.md``). Implemented:
AIRS 3-D, AIRS 2-D (brightness temperature; not compared with MATLAB) and
IASI 3-D.

Inputs and outputs
------------------
``load_nadir_sounders(time_range, instrument, **options)`` -> ``xarray.Dataset``

``time_range``: as ``load_limb_sounders`` (one time or two; datetime64, ISO
strings or datetimes), with the same ``time_handling`` rules. Whole granules
whose mean time lies in the range are loaded. May be None if ``granules`` is
given.

``instrument``: ``"AIRS"`` or ``"IASI"`` (any case).

Options (type, default; MATLAB name):

- ``product`` (str, "3d"): "3d" retrieval; "2d" brightness temperature (AIRS:
  L. Hoffmann's granule files, ``obs._nadir_airs2d``).
- ``time_handling`` (str, "fuzzy"): as ``load_limb_sounders``.
- ``granules`` (list of (date, number), None): load exactly these granules
  (``DateNum``, ``GranuleId``); dates as for ``time_range``.
- ``lat_range``, ``lon_range`` ((min, max), whole globe): keep along-track
  rows with any point in the box.
- ``original_grid`` (bool, False): no horizontal regridding (``LoadOnly``,
  which also implies ``original_z``).
- ``n_cross`` (int, 90; ``NXT``), ``n_along`` (int per granule, 135;
  ``NAT``, which in MATLAB is the total), ``dx_cross``, ``dx_along`` (km,
  None; ``dXT``, ``dAT``; override the numbers).
- ``original_z`` (bool, False; not ``VerticallyInterpolate``): keep the file
  levels.
- ``height_scale`` (array, km, None): heights to interpolate to; default a
  3 km grid from the lowest level (``VerticalSpacing`` 3, as MATLAB).
- ``height_range`` ((min, max), km, (0, inf)): clips the heights.
- ``interp_type`` (str, "linear"; ``Interpolant``): "linear" or "nearest",
  horizontal only (vertical is always linear, as MATLAB).
- ``extrapolate`` (bool, False; ``Extrapolant`` "linear" if True): horizontal
  only.
- ``keep_outliers`` (bool, False): keep values whose state marks them bad
  (the states are recorded either way).
- ``channels`` (list of str, None): 2-D products: which channels or channel
  groups, in that order (AIRS: "4mu", "8mu", "15mu_low", "15mu_high"; default
  all).
- ``screen`` (None): the statistical outlier screen; None uses the
  instrument's default (IASI ``"rolling"``, AIRS ``"none"``). A name, a
  ``(name, settings)`` pair, a function ``screen(temp, time, lat, lon,
  **settings) -> bool``, or ``"none"``; see ``obs._screen``. Recorded in the
  ``screen`` attribute. Data read beyond the request for the screen are
  listed in ``files_read``.
- ``day_night`` (bool, False; ``DayNightFlag``): add ``day_retrieval``.
- ``keep_native_geolocation`` (bool, False; ``KeepOldTime`` for time): also
  return the positions and times before regridding.
- ``additional_vars`` (list of str, ()): extra file variables, [cross x
  along] or [cross x along x level], regridded like temperature.
- ``file_source`` (bool, False): also return ``source_file``, ``source_row``
  and ``original_files``.
- ``file_path`` (str or Path, None): read the files from this folder.
- ``date_warning`` (bool, True), ``verbose`` (bool, False): as
  ``load_limb_sounders``.

**Output**: a Dataset in the swath format, version 1 (see
``data.formats``): dimensions ``(along_track, cross_track, level)``;
``latitude``, ``longitude``, ``time``, ``altitude`` (level), ``temperature``
(float32) and ``temperature_flag`` (data states); ``granule``,
``granule_day``, ``segment`` along track; ``distance_along`` (km from the
start of each segment) and ``distance_cross`` on a regular grid; optional
variables as requested. Variables that do not exist are listed in the
``availability`` attribute (e.g. the satellite position: ``not_in_files``).

Segments: the swath is split where successive rows are more than 180 s apart
(missing granules, or separate time ranges; not 1.5 scans, because the AIRS
daily files store times as float32); each segment is
regridded on its own along-track axis, all on one cross-track axis (that of
the first segment, as MATLAB assumes every scan has the same cross-track
geometry). Rows of different segments are not interpolated across. A scan
with any position missing or out of range, or with no valid time (fill
values in the files, see the readers), cannot be placed: such scans form
segments of their own, kept on the native grid output (unless a
latitude/longitude box is given: they are never inside it) and left out of
the regridded one (as a missing granule is); the cross-track axis is that
of the first segment that can be
placed. If no two successive scans can be placed, the native grid is
returned, with a warning. A stray pixel time alone is ``bad_geolocation``
but does not stop its scan being regridded (its NaN time spreads to the
regridded times around it).

MATLAB -> Python names (private layer ``_load_nadir_arrays``, compared with
MATLAB): ``l1_lat`` Lat, ``l1_lon`` Lon, ``l1_time`` Time (datenums),
``ret_temp`` Temp [cross x along x level], ``ret_z`` Alt, ``Spacing``
Spacing ([cross, along, vertical] km), ``DayNightFlag`` DayNight.

Science
-------
- Regridding and vertical grid: ``swath.regrid``, ``swath.regular_levels``,
  ``swath.interpolate_levels`` (see ``swath``), as MATLAB.
- States (``data.states``): fill values in the files are ``missing``;
  positions outside [-90, 90] / [-180, 180] or not finite, or a missing
  time, are ``bad_geolocation``; temperatures below 100 K or above 400 K are
  ``outside_physical_range`` (IASI also at exactly 100 K, its retrieval's
  lower bound; proposal: MATLAB applies this only to IASI; for AIRS it
  changes nothing where no such values exist, as in the test granules);
  values the instrument's statistical screen marks are
  ``statistical_outlier`` (IASI default ``"rolling"``: more than 5
  interquartile ranges beyond the quartiles of their height level within a
  10-minute window, day and night separately; AIRS none; see
  ``obs._screen``). Unless ``keep_outliers``, bad values are NaN before
  regridding, so a regridded point is NaN if any input it uses is. Regridded
  points take the states of the inputs they use (the OR of the flags of
  every input point with a non-zero weight); NaN points with no other state
  (outside the input) are ``no_observation``.
- Day/night: ``solar.day_retrieval`` (solar zenith angle < 96 degrees, where
  L. Hoffmann's retrieval switches).

Origin
------
Port of ``prep_airs_3d.m`` (corwin365/MatlabFunctions, commit f1a78c0; C.
Wright, 2019-2026), the group's own code.

Validation
----------
``tests/test_nadir_matlab.py``: the private layer ``_load_nadir_arrays``
(MATLAB-like names and arrays, datenum times) against the MATLAB reference
outputs of ``test_prep_airs_3d.m`` (reference file, MATLAB commit and
tolerances in the test's docstring): AIRS 2007-01-01 cases; IASI 2007-11-06
and 07 cases with ``screen="pooled"``, identical wherever Python has a value
(MATLAB fills in every NaN from the nearest point along track, screened
points and NaNs already in the files alike; Python leaves them NaN with their
state). A by-hand sweep at 0.4.1 compared the same layer on 878 real-data
calls (AIRS 3-D 2002-2025, IASI 3-D 2007-2018): no unexplained differences.
Not compared with MATLAB: the AIRS granule-file reader, the default IASI
screen, AIRS 2-D, several options, and the fill handling added after the
sweep. Full record: ``docs/validation.md`` in the repository.
``tests/test_nadir_known.py``, ``tests/test_nadir_iasi.py``,
``tests/test_nadir_airs2d.py``, ``tests/test_nadir_2d.py`` and
``tests/test_nadir_fill.py``: known answers and format checks on synthetic
files.

Status
------
Preliminary. AI-ported (Claude) from the MATLAB; to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.obs import load_nadir_sounders  # doctest: +SKIP
>>> ds = load_nadir_sounders("2007-01-01T05:00", "AIRS",
...                          time_handling="strict")  # doctest: +SKIP

Change log
----------
2026/11/25 [Claude] initial version (AIRS 3-D)
2026/11/25 [Claude] IASI 3-D: orbit files, pseudo-granules, outlier screen (no refill)
2026/11/25 [Claude] swappable screens; 2-D products; AIRS 2-D reader; channels option
2026/11/25 [Claude] shared helpers moved to obs/_common.py (outputs unchanged)
2026/11/25 [Claude] docstring: IASI screen, segment rule, duplicate option fixed
2026/11/25 [Claude] IASI: exactly 100 K (retrieval lower bound) outside_physical_range
2026/11/25 [Claude] extrapolated state; record-based attributes; heights and
  preprocessing attributes (preprocessing as the readers record it)
2026/11/26 [Claude] missing time is bad_geolocation; unplaceable scans not regridded
2026/11/26 [Claude] docstring: AIRS 2-D listed as implemented; Validation scope

"""

import re
import warnings

import numpy as np
import xarray as xr

from atmosgw.curate import swath
from atmosgw.curate.data import states
from atmosgw.curate.data.formats import FORMATS
from atmosgw.curate.matlab_time import datenum, datenum_to_datetime
from atmosgw.curate.obs import _nadir_airs, _nadir_airs2d, _nadir_iasi, _screen
from atmosgw.curate.obs._common import (
    COORD_ATTRS,
    check_additional_vars,
    check_flags,
    check_time_handling,
    day_label,
    expand_time_range,
    file_versions,
    find_instrument,
    load_attrs,
    outside_physical_range,
    record,
    record_attrs,
    to_datenums,
    warn_if_no_data,
    warn_mixed_versions,
)
from atmosgw.curate.paths import dataset_dir
from atmosgw.curate.solar import day_retrieval

_FAR = datenum(9999, 999, 999)

# instruments: dates covered (UNVERIFIED start dates), reader and geometry
INSTRUMENTS = {
    "AIRS": {"time_range": (datenum(2002, 8, 30), _FAR),
             "reader": {"3d": _nadir_airs, "2d": _nadir_airs2d},
             "dataset_kind": {"3d": "raw", "2d": "2d"},
             "dataset": "AIRS", "n_cross": 90, "rows_per_granule": 135,
             "row_seconds": 360 / 135, "max_row_gap_seconds": 180.0,
             "products": ("3d", "2d"),
             "characterisation": {"3d": "AIRS-3D-placeholder/0",
                                  "2d": "AIRS-2D-placeholder/0"},
             "screen": ("none", {}),
             "not_in_files": ("satellite_latitude", "satellite_longitude",
                              "satellite_altitude", "view_zenith_angle",
                              "view_azimuth_angle")},
}

# IASI: n_cross and rows_per_granule are the regridding defaults (MATLAB NXT 90,
# NAT 135 per granule, as for AIRS), not the native 60 x 90 (proposal: as MATLAB)
INSTRUMENTS["IASI"] = {
    "time_range": (datenum(2007, 5, 29), _FAR), "reader": _nadir_iasi,
    "dataset": "IASI", "n_cross": 90, "rows_per_granule": 135,
    "row_seconds": 4.0, "max_row_gap_seconds": 180.0, "products": ("3d",),
    "characterisation": {"3d": "IASI-3D-placeholder/0"},
    "screen": ("rolling", {}),
    # exactly 100 K is the retrieval's lower bound, not a measurement: the
    # files are float32 and never below 100.0 K, and 100.0 occurs in whole
    # bad profiles (2010-07-01: 30,805 points in 6,467 profiles); removed as
    # outside_physical_range (C. Wright, 2026-11-25; UNVERIFIED with L. Hoffmann)
    "floor_is_bad": True,
    "not_in_files": INSTRUMENTS["AIRS"]["not_in_files"],
}


def load_nadir_sounders(
    time_range,
    instrument,
    *,
    product="3d",
    time_handling="fuzzy",
    granules=None,
    lat_range=(-90, 90),
    lon_range=(-180, 180),
    original_grid=False,
    n_cross=None,
    n_along=None,
    dx_cross=None,
    dx_along=None,
    original_z=False,
    height_scale=None,
    height_range=(0, np.inf),
    interp_type="linear",
    extrapolate=False,
    keep_outliers=False,
    screen=None,
    channels=None,
    day_night=False,
    keep_native_geolocation=False,
    additional_vars=(),
    file_source=False,
    file_path=None,
    date_warning=True,
    verbose=False,
):
    """Load nadir-sounder swaths as an xarray.Dataset. See the module docstring."""
    options = {k: v for k, v in locals().items() if k not in ("time_range",
                                                              "instrument")}
    call = {"time_range": time_range, "instrument": instrument, **options}
    data, settings = _load_nadir_arrays(time_range, instrument, **options)
    return _to_dataset(data, settings, call)


def _load_nadir_arrays(time_range, instrument, **options):
    """
    The MATLAB prep_airs_3d, step for step (without detrending): returns
    (data, settings). data has MATLAB-like names and [cross x along (x level)]
    arrays, times as datenums; compared directly with MATLAB by the tests.
    """
    settings, inst_info, info = _prepare(time_range, instrument, options)
    say = print if settings["verbose"] else (lambda *a, **k: None)

    say("--> Loading data")
    var_names = (["Lat", "Lon", "Time", "Temp"]
                 + (["Alt"] if settings["product"] == "3d" else [])
                 + settings["additional_vars"])
    screen, screen_settings, screen_text, pad = _screen.resolve(settings["screen"],
                                                               info["screen"])
    settings["screen_text"] = screen_text
    # with a screen that needs neighbouring data, read beyond the request
    # (then trim), so a scan's flags do not depend on how it was requested
    load_settings, padded = settings, screen is not None and pad > 0
    if padded:
        load_settings = dict(settings)
        if settings["granules"] is not None:
            load_settings["granules"] = _with_neighbours(settings["granules"])
        else:
            t0, t1 = settings["time_range"]
            load_settings["time_range"] = np.array([t0 - pad / 1440, t1 + pad / 1440])
    data, file_list = inst_info["reader"].load(load_settings, inst_info, var_names)
    settings["record"] = record(inst_info)
    settings["files_read"] = list(file_list)
    warn_mixed_versions(settings, stacklevel=3)
    if padded and settings["granules"] is not None:  # padding is not a request
        extra = {f"{day_label(d)} granule {g}"
                 for d, g in set(load_settings["granules"]) - set(settings["granules"])}
        settings["record"]["files_not_found"][:] = [
            f for f in settings["record"]["files_not_found"]
            if _granule_label(f) not in extra]

    # states, then bad values to NaN (unless keep_outliers)
    temp = np.asarray(data["Temp"], dtype=float)
    lat, lon = np.asarray(data["Lat"], float), np.asarray(data["Lon"], float)
    flags = states.new_flags(temp.shape)
    states.set_state(flags, np.isnan(temp), "missing")
    bad_geo = ~(np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90)
                & (np.abs(lon) <= 180) & np.isfinite(np.asarray(data["Time"], float)))
    states.set_state(flags, np.broadcast_to(bad_geo[..., None], temp.shape),
                     "bad_geolocation")
    # instrument statistical screen (swappable: see obs._screen), on the padded
    # data, before the range check as MATLAB
    if screen is not None:
        states.set_state(flags, screen(temp, np.asarray(data["Time"], float), lat, lon,
                                       **screen_settings), "statistical_outlier")
    with np.errstate(invalid="ignore"):
        implausible = outside_physical_range(
            temp, floor_is_bad=info.get("floor_is_bad", False))
        states.set_state(flags, implausible, "outside_physical_range")
    if not settings["keep_outliers"]:
        temp[states.no_value_mask(flags)] = np.nan
    data["Temp"] = temp
    if padded:
        keep = _requested_rows(data, settings)
        for name in _along_track_names(data):
            data[name] = np.take(data[name], keep, axis=1 if np.ndim(data[name]) > 1
                                 else 0)
        flags = flags[:, keep]
        lat, lon = data["Lat"], data["Lon"]
    if (settings["scalar_time"] and settings["time_handling"] == "strict"
            and data["Lat"].shape[1] == 0):
        warnings.warn("A single time_range value with time_handling='strict' selects "
                      "only granules whose mean time is exactly that instant, and none "
                      "were found. Use 'fuzzy' (default) or two times.", stacklevel=3)
    if settings["time_range"] is None:  # granules given: record their times
        t = data["Time"]
        settings["time_range"] = ((np.nanmin(t), np.nanmax(t)) if t.size
                                  else (np.nan, np.nan))

    # segments: break where successive rows are more than max_row_gap_seconds
    # apart (180 s for AIRS: the daily files store times as float32, in steps of
    # 16 s in 2007 to 64 s in 2025, so a 1.5-scan rule would split every
    # granule; a missing granule leaves a 6 min gap)
    n_rows = lat.shape[1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        row_time = np.nanmean(np.asarray(data["Time"], float), axis=0)
    gap = np.diff(row_time) * 86400 > info["max_row_gap_seconds"]
    # and around scans that cannot be placed (any position or time missing)
    placeable = _placeable(data)
    gap |= placeable[1:] != placeable[:-1]
    data["Segment"] = 1.0 + np.concatenate([[0.0], np.cumsum(gap)]) if n_rows else (
        np.empty(0))
    if settings["keep_native_geolocation"]:
        data["NativeLat"], data["NativeLon"] = lat.copy(), lon.copy()
        data["NativeTime"] = np.asarray(data["Time"], float).copy()
        data["ScanPosition"] = np.arange(1.0, lat.shape[0] + 1)

    # horizontal regridding, segment by segment
    spacing = [np.nan, np.nan, np.nan]
    regrid = not settings["original_grid"] and _regriddable(data["Segment"],
                                                             placeable)
    if not settings["original_grid"] and not regrid and np.any(~placeable):
        warnings.warn("No two successive scans have positions and times (fill "
                      "values in the files): the native grid is returned.",
                      stacklevel=3)
    if regrid:
        say("--> Regridding")
        data, flags, spacing[:2] = _regrid(data, flags, settings, info)
    settings["grid"] = "regular" if regrid else "native"

    # vertical grid (3-D products; a 2-D product's last dimension is channels)
    z = (np.asarray(data["Alt"], dtype=float).ravel() if settings["product"] == "3d"
         else None)
    if z is None:
        pass
    elif not settings["original_z"]:
        new_z = (swath.regular_levels(z, 3.0) if settings["height_scale"] is None
                 else settings["height_scale"])
        low, high = settings["height_range"]
        new_z = new_z[(new_z >= low) & (new_z <= high)]
        temp, flags = _interpolate_with_states(z, data["Temp"], flags, new_z)
        data["Temp"], data["Alt"] = temp, new_z
        steps = np.diff(new_z)
        if steps.size and np.allclose(steps, steps[0], rtol=0, atol=1e-9):
            spacing[2] = float(steps[0])
    else:
        low, high = settings["height_range"]
        keep = (z >= low) & (z <= high)
        data["Temp"], data["Alt"], flags = data["Temp"][..., keep], z[keep], (
            flags[..., keep])
    data["Spacing"] = np.array(spacing)

    # rows entering the latitude/longitude box
    (la0, la1), (lo0, lo1) = settings["lat_range"], settings["lon_range"]
    if (la0, la1, lo0, lo1) != (-90, 90, -180, 180) and data["Lat"].shape[1]:
        inside = ((data["Lat"] >= min(la0, la1)) & (data["Lat"] <= max(la0, la1))
                  & (data["Lon"] >= min(lo0, lo1)) & (data["Lon"] <= max(lo0, lo1)))
        keep = np.flatnonzero(np.any(inside, axis=0))
        for name in _along_track_names(data):
            data[name] = np.take(data[name], keep, axis=1 if data[name].ndim > 1
                                 else 0)
        flags = flags[:, keep]

    if settings["day_night"]:
        data["DayNight"] = day_retrieval(data["Time"], data["Lat"], data["Lon"],
                                         datenum=True)
    data["Flags"] = {"Temp": flags}
    return data, settings


# supporting functions


def _prepare(time_range, instrument, options):
    """Check the arguments; return (settings, inst_info, info)."""
    check_flags(options, ("original_grid", "original_z", "extrapolate",
                          "keep_outliers", "day_night", "keep_native_geolocation",
                          "file_source", "date_warning", "verbose"))
    check_time_handling(options)
    name = find_instrument(instrument, INSTRUMENTS)
    info = INSTRUMENTS[name]
    product = str(options["product"]).lower()
    if product not in info["products"]:
        raise NotImplementedError(f"{name} product {product!r} is not available "
                                  f"yet (available: {', '.join(info['products'])}).")
    if options["interp_type"].lower() not in ("linear", "nearest"):
        raise ValueError("interp_type must be 'linear' or 'nearest'.")
    _screen.resolve(options["screen"], info["screen"])  # (checks the choice)
    channels = options["channels"]
    if channels is not None and (isinstance(channels, str)
                                 or not all(isinstance(c, str) for c in channels)):
        raise TypeError("channels must be a list of channel names.")
    check_additional_vars(options)
    additional = options["additional_vars"]

    granules = options["granules"]
    if granules is not None:
        granules = [(int(np.floor(to_datenums(d)[0])), int(g)) for d, g in granules]
        expanded, scalar = None, False
    else:
        if time_range is None:
            raise ValueError("Give time_range, or granules.")
        raw = to_datenums(time_range)
        if options["date_warning"]:
            warn_if_no_data(raw, name, info["time_range"], stacklevel=3)
        expanded, scalar = expand_time_range(raw, options["time_handling"])

    height_scale = options["height_scale"]
    if height_scale is not None:
        height_scale = np.asarray(height_scale, dtype=float).ravel()
    settings = dict(options, instrument=name, product=product, granules=granules,
                    time_range=expanded, scalar_time=scalar,
                    height_scale=height_scale,
                    height_range=tuple(float(h) for h in options["height_range"]),
                    lat_range=tuple(float(x) for x in options["lat_range"]),
                    lon_range=tuple(float(x) for x in options["lon_range"]),
                    additional_vars=list(additional),
                    interp_type=options["interp_type"].lower())
    if settings["original_grid"]:
        settings["original_z"] = True  # as MATLAB LoadOnly
    reader = info["reader"]
    if isinstance(reader, dict):  # a reader per product
        reader = reader[product]
    kind = info.get("dataset_kind", {}).get(product, "raw")
    inst_info = {"inst": name, "reader": reader}
    inst_info["path"] = (str(options["file_path"]) if options["file_path"] is not None
                         else str(dataset_dir(info["dataset"], kind,
                                              must_exist=False)))
    return settings, inst_info, info


def _regrid(data, flags, settings, info):
    """Regrid every segment; returns (data, flags, (cross km, along km))."""
    lat, lon = data["Lat"], data["Lon"]
    segment = data["Segment"]
    placeable = _placeable(data)
    first = np.flatnonzero(segment == segment[placeable][0])
    d0 = swath.swath_distances(lat[:, first], lon[:, first])
    cross_axis = swath.regular_axis(d0["cross"],
                                    settings["n_cross"] or info["n_cross"],
                                    settings["dx_cross"])
    extra = [n for n in settings["additional_vars"] + list(data.get("DerivedVars", {}))
             if np.ndim(data[n]) >= 2 and np.shape(data[n])[:2] == lat.shape]
    parts, spacing = [], (np.nan, np.nan)
    for seg in np.unique(segment):
        rows = np.flatnonzero(segment == seg)
        if rows.size < 2 or not placeable[rows[0]]:
            continue  # a single row, or scans without positions: not regridded
        n_gran = len({(d, g) for d, g in zip(data["Day"][rows], data["Granule"][rows],
                                             strict=True)})
        n_along = (None if settings["dx_along"] is not None
                   else (settings["n_along"] or info["rows_per_granule"]) * n_gran)
        fields = {"Temp": data["Temp"][:, rows], "Time": data["Time"][:, rows],
                  **{n: data[n][:, rows] for n in extra}}
        out = swath.regrid(lat[:, rows], lon[:, rows], fields, n_along=n_along,
                           dx_along=settings["dx_along"],
                           method=settings["interp_type"],
                           extrapolate="linear" if settings["extrapolate"] else None,
                           cross_axis=cross_axis, flags=flags[:, rows])
        # states: those of every input point used (swath.regrid), and
        # no_observation for NaN points outside the input
        new_flags = out["flags"]
        unexplained = np.isnan(out["fields"]["Temp"]) & ~states.no_value_mask(new_flags)
        states.set_state(new_flags, unexplained, "no_observation")
        d_seg = swath.swath_distances(lat[:, rows], lon[:, rows])
        if settings["extrapolate"]:
            # outside this segment's own distances (e.g. cross-track, where the
            # common cross axis is wider than the segment): extrapolated, where
            # there is a value
            cross_out = ((cross_axis < np.nanmin(d_seg["cross"]))
                         | (cross_axis > np.nanmax(d_seg["cross"])))
            along_out = ((out["along"] < np.nanmin(d_seg["along"]))
                         | (out["along"] > np.nanmax(d_seg["along"])))
            outside = cross_out[:, None] | along_out[None, :]
            outside = np.broadcast_to(
                outside.reshape(outside.shape + (1,) * (new_flags.ndim - 2)),
                new_flags.shape)
            states.set_state(new_flags,
                             outside & np.isfinite(out["fields"]["Temp"]),
                             "extrapolated")
        # per-row information: from the nearest input row along track
        d_in = d_seg["along"]
        nearest = rows[np.clip(np.searchsorted(d_in, out["along"]), 0, rows.size - 1)]
        part = {"Lat": out["lat"], "Lon": out["lon"],
                "Time": out["fields"]["Time"], "Temp": out["fields"]["Temp"],
                **{n: out["fields"][n] for n in extra},
                "Flags": new_flags, "DistAlong": out["along"],
                "Segment": np.full(out["along"].size, float(seg))}
        for name in ("Granule", "Day", "SourceFile", "SourceRow"):
            part[name] = np.asarray(data[name])[nearest]
        parts.append(part)
        if np.isnan(spacing[1]):
            spacing = out["spacing"]
    new = dict(data)
    for name in parts[0]:
        if name in ("Flags",):
            continue
        new[name] = np.concatenate([p[name] for p in parts],
                                   axis=1 if parts[0][name].ndim > 1 else 0)
    new["DistCross"] = cross_axis
    for name in ("NativeLat", "NativeLon", "NativeTime", "ScanPosition"):
        if name in data:
            new[name] = data[name]
    return new, np.concatenate([p["Flags"] for p in parts], axis=1), spacing


def _placeable(data):
    """True for each along-track row whose positions are all valid (finite and
    in range, as for bad_geolocation) and which has a valid time (a stray
    pixel time alone does not stop regridding)."""
    lat, lon = np.asarray(data["Lat"], float), np.asarray(data["Lon"], float)
    with np.errstate(invalid="ignore"):
        ok = (np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90)
              & (np.abs(lon) <= 180))
    has_time = np.any(np.isfinite(np.asarray(data["Time"], float)), axis=0)
    return np.all(ok, axis=0) & has_time


def _granule_label(entry):
    """'<day> granule g' at the start of a files_not_found entry, else entry."""
    match = re.match(r"^(\S+ granule \d+)(?=$| |:)", entry)
    return match.group(1) if match else entry


def _regriddable(segment, placeable):
    """True if some segment of placeable rows has at least two rows."""
    ids, counts = np.unique(np.asarray(segment)[placeable], return_counts=True)
    return bool(np.any(counts >= 2))


def _interpolate_with_states(z, values, flags, new_z):
    """Linear interpolation to new_z; each new level takes the states of the
    levels used (swath.interpolate_level_flags)."""
    out = swath.interpolate_levels(z, values, new_z)
    new_flags = swath.interpolate_level_flags(z, flags, new_z)
    unexplained = np.isnan(out) & ~states.no_value_mask(new_flags)
    states.set_state(new_flags, unexplained, "no_observation")
    return out, new_flags


def _with_neighbours(granules):
    """The granules and the one before and after each (for screen padding)."""
    out = set(granules)
    for day, g in granules:
        for step in (-1, 1):
            n = g + step
            if n < 1:
                out.add((day - 1, 240))
            elif n > 240:
                out.add((day + 1, 1))
            else:
                out.add((day, n))
    return sorted(out)


def _requested_rows(data, settings):
    """Indices of the along-track rows of the granules actually requested."""
    key = np.asarray(data["Day"]) * 1000 + np.asarray(data["Granule"])
    if settings["granules"] is not None:
        wanted = [d * 1000 + g for d, g in settings["granules"]]
        return np.flatnonzero(np.isin(key, wanted))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        row_time = np.nanmean(np.asarray(data["Time"], float), axis=0)
    ids, inverse = np.unique(key, return_inverse=True)
    sums = np.bincount(inverse, weights=np.nan_to_num(row_time))
    counts = np.bincount(inverse, weights=np.isfinite(row_time).astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sums / counts
    t0, t1 = settings["time_range"]
    return np.flatnonzero(((mean >= t0) & (mean <= t1))[inverse])


def _along_track_names(data):
    """Names of the arrays with an along-track dimension."""
    n = data["Lat"].shape[1]
    names = []
    for name, value in data.items():
        skip = ("Flags", "DistCross", "Alt", "Spacing", "Channel", "ChannelName",
                "DerivedVars")
        if name in skip or name.startswith(
                "Native") or name == "ScanPosition":
            continue
        shape = np.shape(value)
        if (len(shape) >= 2 and shape[1] == n) or (len(shape) == 1 and shape[0] == n):
            names.append(name)
    return names


_ATTRS = {
    **COORD_ATTRS,
    "temperature": {"units": "K", "standard_name": "air_temperature",
                    "origin": "retrieved", "ancillary_variables": "temperature_flag"},
    "brightness_temperature": {"units": "K", "long_name": "brightness temperature",
                               "origin": "measured",
                               "ancillary_variables": "brightness_temperature_flag"},
    "channel_wavenumber": {"units": "cm-1", "long_name": "channel centre wavenumber"},
    "distance_along": {"units": "km",
                       "long_name": "along-track distance from the segment start"},
    "distance_cross": {"units": "km",
                       "long_name": "cross-track distance from the middle of the "
                                    "first scan of the first segment"},
    "granule": {"long_name": "granule number within its day (1-based)"},
    "granule_day": {"long_name": "day of the granule"},
    "segment": {"long_name": "unbroken stretch of the swath (1-based)"},
    "source_file": {"long_name": "index into original_files (1-based)"},
    "source_row": {"long_name": "nearest along-track row within its granule "
                                "(1-based)"},
    "day_retrieval": {"long_name": "1: retrieval's day set-up (solar zenith "
                                   "angle < 96 degrees), 0: night"},
    "native_latitude": {"units": "degrees_north"},
    "native_longitude": {"units": "degrees_east"},
    "scan_position": {"long_name": "cross-track position in the scan (1-based)"},
}


def _to_dataset(data, settings, call=None):
    """Arrange the arrays as a Dataset in the swath format."""
    ac = ("along_track", "cross_track")
    flags = data["Flags"]["Temp"].transpose(1, 0, 2)
    regular = settings["grid"] == "regular"
    three_d = settings["product"] == "3d"
    # the main variable: temperature on levels (3-D), or brightness temperature
    # (2-D), with a channel dimension unless there is one channel and no list
    main = "temperature" if three_d else "brightness_temperature"
    if three_d:
        main_dims = ac + ("level",)
    elif data["Temp"].shape[2] == 1 and "Channel" not in data:
        main_dims = ac
        flags = flags[..., 0]
    else:
        main_dims = ac + ("channel",)

    def t(x):  # [cross x along (x ...)] -> [along x cross (x ...)]
        x = np.asarray(x)
        return np.swapaxes(x, 0, 1) if x.ndim >= 2 else x

    def times(x):
        return datenum_to_datetime(np.asarray(x, dtype=float), "us")

    v = {
        "latitude": (ac, t(data["Lat"]).astype(np.float64)),
        "longitude": (ac, t(data["Lon"]).astype(np.float64)),
        "time": (ac, times(t(data["Time"]))),
        main: (main_dims, t(data["Temp"]).astype(np.float32).reshape(
            flags.shape)),
        f"{main}_flag": (main_dims, flags, states.flag_attrs()),
        "granule": (("along_track",), np.asarray(data["Granule"]).astype(np.int32)),
        "granule_day": (("along_track",), times(data["Day"])),
        "segment": (("along_track",), np.asarray(data["Segment"]).astype(np.int32)),
    }
    if not three_d and "ChannelName" in data and main_dims[-1] == "channel":
        v["channel_name"] = (("channel",), np.asarray(data["ChannelName"], dtype=str),
                             {"long_name": "channel or channel group"})
    # variables the files provide already derived (e.g. a producer's
    # perturbation): the main variable's states, plus "derived"
    derived_outputs = []
    for key, spec in data.get("DerivedVars", {}).items():
        out = spec["output"]
        value = t(data[key]).astype(np.float32).reshape(flags.shape)
        if not settings["keep_outliers"]:
            value[states.no_value_mask(flags)] = np.nan
        dflag = flags.copy()
        # (qualified states only where this variable has a value)
        for name in ("extrapolated", "filled"):
            bit = states.FLAG_DTYPE(states.state(name).bit)
            dflag[np.isnan(value)] &= ~bit
        states.set_state(dflag, np.isfinite(value), "derived")
        states.set_state(dflag, np.isnan(value) & ~states.no_value_mask(dflag),
                         "no_observation")
        v[out] = (main_dims, value, {"units": spec["units"],
                                     "long_name": spec["long_name"],
                                     "origin": "derived",
                                     "ancillary_variables": f"{out}_flag",
                                     "comment": spec["comment"]})
        v[f"{out}_flag"] = (main_dims, dflag, states.flag_attrs())
        derived_outputs.append(out)
    if three_d:
        v["altitude"] = (("level",), np.asarray(data["Alt"], dtype=np.float64))
    elif "Channel" in data:
        v["channel_wavenumber"] = (("channel",), np.asarray(data["Channel"],
                                                           dtype=np.float64))
    if regular:
        v["distance_along"] = (("along_track",), np.asarray(data["DistAlong"],
                                                            dtype=np.float64))
        v["distance_cross"] = (("cross_track",), np.asarray(data["DistCross"],
                                                            dtype=np.float64))
    if settings["file_source"]:
        v["source_file"] = (("along_track",),
                            np.asarray(data["SourceFile"]).astype(np.int32))
        v["source_row"] = (("along_track",),
                           np.asarray(data["SourceRow"]).astype(np.int32))
    if settings["day_night"]:
        v["day_retrieval"] = (ac, t(data["DayNight"]).astype(np.float32))
    if settings["keep_native_geolocation"]:
        nat = ("native_along", "native_cross")
        v["native_latitude"] = (nat, t(data["NativeLat"]).astype(np.float64))
        v["native_longitude"] = (nat, t(data["NativeLon"]).astype(np.float64))
        v["native_time"] = (nat, times(t(data["NativeTime"])))
        v["scan_position"] = (("native_cross",),
                              np.asarray(data["ScanPosition"]).astype(np.int32))
    for name in settings["additional_vars"]:
        value = np.asarray(data[name])
        if value.ndim >= 2 and value.shape[:2] == np.shape(data["Lat"]):
            v[name] = (ac + tuple(f"{name}_dim{d}" for d in range(2, value.ndim)),
                       t(value))
    ds = xr.Dataset({k: (val[0], val[1], val[2] if len(val) > 2 else
                         dict(_ATTRS.get(k, {}))) for k, val in v.items()})
    # where each variable comes from, from the readers' record (as limb)
    for name, out in (("Lat", "latitude"), ("Lon", "longitude"), ("Time", "time"),
                      ("Alt", "altitude"), ("Temp", main)):
        if out in ds:
            ds[out].attrs.update(record_attrs(
                name, settings["record"],
                interpolated=three_d and not settings["original_z"]))
    if settings["file_source"]:
        ds = ds.assign_coords(
            original_files=("file", np.array(settings["files_read"], dtype=str)),
            original_file_version=("file", np.array(file_versions(
                settings["record"], settings["files_read"]), dtype=str),
                {"long_name": "data version stated by each file"}))

    info = INSTRUMENTS[settings["instrument"]]
    tracked = [n for n, s in FORMATS["swath"]["variables"].items() if s.get("tracked")]
    availability = {n: ("not_in_files" if n in info["not_in_files"]
                        else "not_requested") for n in tracked if n not in ds}
    availability[f"{main}_uncertainty"] = "not_provided"
    for out in derived_outputs:
        availability[f"{out}_uncertainty"] = "not_provided"
    flags_by_variable = {main: flags, **{out: v[f"{out}_flag"][1]
                                         for out in derived_outputs}}
    common = load_attrs(settings, call, "atmosgw.curate.obs.load_nadir_sounders")
    # what was done to the files before loading, as the readers record it
    preprocessing = "; ".join(sorted(settings["record"].get("preprocessing", ())))
    ds.attrs = {
        "format": "swath", "format_version": FORMATS["swath"]["version"],
        "instrument": settings["instrument"], "product": settings["product"],
        "characterisation": info["characterisation"][settings["product"]],
        "data": "preprocessed" if preprocessing else "raw",
        "preprocessing": preprocessing,
        "heights": (("original" if settings["original_z"] else "interpolated")
                    if three_d else "none"),
        "time_range": common.pop("time_range"),
        "availability": states.format_availability(availability),
        **common,
        "state_counts": states.counts_text(flags_by_variable),
        "spatial_dims": "along_track cross_track" + (" level" if three_d else ""),
        "grid_spacing_km": [float(data["Spacing"][1]), float(data["Spacing"][0]),
                            float(data["Spacing"][2])],
        "grid": settings["grid"],
        "screen": settings["screen_text"],
    }
    return ds


__all__ = ["load_nadir_sounders", "INSTRUMENTS"]
