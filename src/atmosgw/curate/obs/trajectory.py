"""
Load in-situ measurements along a path (radiosondes; later aircraft and
balloons) in a common format.

Each flight (an ascent, a flight) is one row of a [flight x sample] grid,
its samples in path order with their own position, time and height, so
that drift is kept. This is the "trajectory" format
(``docs/design/data_formats.md``), version 1: it uses the variable names,
state flags and attributes of the other formats, without a level grid.
First instrument: the DWD high-resolution radiosondes ("DWDsondes").

Inputs and outputs
------------------
``load_trajectories(time_range, instrument, **options)`` -> ``xarray.Dataset``

``time_range``: one time (that day) or two (start and end), as
``numpy.datetime64``, ISO strings or ``datetime`` objects, as for
``load_limb_sounders`` (numbers are refused). A flight is selected by its
start (launch) time.

``instrument``: ``"DWDsondes"`` (any case, e.g. ``"dwdsondes"``): the DWD
high-resolution radiosondes. Each sonde dataset gets its own name (e.g. a
later ``"IGRAsondes"``), as the files and their meanings differ.

Options (type, default):

- ``stations`` (station ids or names, or a list of them, None): only these
  stations (e.g. ``1757``, ``"01757"`` or ``"Greifswald"``); all if None.
  Stations asked for but not in the files are listed in ``files_not_found``.
- ``time_handling`` (str, "fuzzy"): "strict", "generous" or "fuzzy", as for
  ``load_limb_sounders``.
- ``keep_outliers`` (bool, False): keep temperatures outside the physical
  range (their state is still set).
- ``additional_vars`` (list of str, ()): extra per-sample variables, by their
  names in the files (DWD: ``AE_RF`` relative humidity %, ``AE_TD`` dew point
  degrees C, ``AE_FF``, ``AE_DD``, ``QN_1``, ``MESSZEITPUNKT``), as read
  (-999 as NaN).
- ``file_source`` (bool, False): also return ``source_file`` and
  ``original_files``.
- ``file_path`` (str or Path, None): read the files from this folder.
- ``verbose`` (bool, False), ``date_warning`` (bool, True).

**Output**: an ``xarray.Dataset`` with dimensions ``flight`` and ``sample``
(and ``file`` with ``file_source``), flights in order of station then start
time, samples padded with NaN (NaT) to the longest flight:

- ``latitude``, ``longitude`` (flight, sample): degrees;
- ``time`` (flight, sample): ``datetime64[us]``;
- ``altitude`` (flight, sample): km (DWD: geopotential height);
- ``pressure`` (hPa), ``temperature`` (K), ``eastward_wind``,
  ``northward_wind`` (m s-1) (flight, sample), each with ``<name>_flag``;
- ``flight_id``, ``platform_id`` (station), ``platform_name``,
  ``platform_type`` (sonde type), ``start_time`` (launch),
  ``reference_time`` (the nominal synoptic time) (flight);
- ``source_file`` (flight): 1-based index into ``original_files``.

Data variables are float64. States (shared rules, ``data_formats.md`` 6.3):
``no_observation`` (padding after a flight's last sample, and before its
first value), ``missing`` (NaN between a variable's first and last value in
the flight), ``outside_physical_range`` (temperature < 100 K, or > 400 K
below 100 km; removed unless ``keep_outliers``; only the temperature is
removed, as the sensors are separate). ``origin`` "measured".
``availability``: ``upward_wind`` not_in_files, the uncertainties
not_provided. Attributes: as the other loaders (``format`` "trajectory",
``format_version``, ``instrument``, ``product`` "sounding",
``characterisation``, ``data`` "raw", ``preprocessing``, ``time_range``,
``availability``, ``source``, ``atmosgw_curate_version``, ``options``,
``files_read``, ``source_versions``, ``files_not_found``,
``files_without_data``, ``multiple_matches``, ``state_counts``).

Science
-------
DWD: see ``obs/_traj_dwd_sondes.py`` for the meanings of the file fields,
which are inferred from the files (no description is published with them):
times from the launch time plus seconds since launch (UTC, checked against
ERA5); winds from speed and the direction the wind blows from (checked
against the balloons' drift); positions of the sonde, or of the station
where the files give none (early years). IAGOS and FAAM: see
``obs/_traj_iagos.py`` and ``obs/_traj_faam.py``. Comparisons with ERA5
are in ``docs/validation.md``.

Origin
------
New code (no MATLAB equivalent).

Validation
----------
``tests/test_new_loaders.py`` (synthetic files: parsing, bisection against
a full read, times, winds, states, options, repeated files). Real files
checked by hand on 1993-06-15, 2004-01-01, 2020-01-07 and 2024-12-31.

Status
------
Proposal (format version 1 is a draft for agreement). AI-originated
(Claude), from the design agreed with C. Wright on 2026-11-27; to be
reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.obs import load_trajectories  # doctest: +SKIP
>>> ds = load_trajectories("2020-01-07", "DWDsondes",  # doctest: +SKIP
...                        stations="Greifswald")

Change log
----------
2026/11/27 [Claude] initial version (DWD radiosondes)
2026/11/27 [Claude] IAGOS and FAAM; per-variable instrument flags; upward_wind

"""

import os
import warnings

import numpy as np
import xarray as xr

from atmosgw.curate.data import states
from atmosgw.curate.data.formats import FORMATS
from atmosgw.curate.matlab_time import datenum, datenum_to_datetime
from atmosgw.curate.obs import _traj_dwd_sondes, _traj_faam, _traj_iagos
from atmosgw.curate.obs._common import (
    COORD_ATTRS,
    check_additional_vars,
    check_flags,
    check_time_handling,
    expand_time_range,
    file_versions,
    find_instrument,
    load_attrs,
    outside_physical_range,
    record,
    record_attrs,
    to_datenums,
    warn_if_no_data,
    warn_if_nothing_in_range,
    warn_missing_files,
    warn_mixed_versions,
    warn_multiple_matches,
)
from atmosgw.curate.obs.limb import _state_flags
from atmosgw.curate.paths import dataset_dir

_FAR = datenum(9999, 999, 999)  # "still running"

_SONDE_TEXT = {
    "flight_id": "flight identifier: <station>-<launch time>",
    "platform_id": "station identifier",
    "platform_name": "station name",
    "platform_type": "sonde type, as in the files' metadata",
}
_AIRCRAFT_TEXT = {
    "flight_id": "flight identifier, as in the files",
    "platform_id": "aircraft registration or identifier",
    "platform_name": "aircraft, as described in the files",
    "platform_type": "kind of platform",
}

# instruments: dates with data, reader, product, the main variables the
# reader gives (reader names, see _OUTPUTS), the format's variables the files
# do not have, per-flight text and times, and the options only they take
INSTRUMENTS = {
    "DWDsondes": {"time_range": (datenum(1990, 1, 1), _FAR),
                  "reader": _traj_dwd_sondes, "product": "sounding",
                  "variables": ("Pres", "Temp", "U", "V"),
                  "not_in_files": ("upward_wind",),
                  "per_flight": _SONDE_TEXT,
                  "flight_times": {"reference_time": "nominal (synoptic) reference "
                                                     "time of the ascent"},
                  "options": ("stations",)},
    "IAGOS": {"time_range": (datenum(1994, 8, 1), _FAR),
              "reader": _traj_iagos, "product": "flight",
              "variables": ("Pres", "Temp", "U", "V"),
              "not_in_files": ("upward_wind",),
              "per_flight": {**_AIRCRAFT_TEXT,
                             "departure_airport": "departure airport",
                             "arrival_airport": "arrival airport",
                             "temperature_source": "file variable used for "
                                                   "temperature",
                             "altitude_source": "file variable used for altitude"},
              "flight_times": {}, "options": ()},
    "FAAM": {"time_range": (datenum(2025, 7, 1), _FAR),
             "reader": _traj_faam, "product": "flight",
             "variables": ("Pres", "Temp", "U", "V", "W"),
             "not_in_files": (),
             "per_flight": {**_AIRCRAFT_TEXT},
             "flight_times": {}, "options": ("sample_rate",)},
}

# reader names -> output names and attributes
_OUTPUTS = {
    "Pres": ("pressure", {"units": "hPa", "standard_name": "air_pressure"}),
    "Temp": ("temperature", {"units": "K", "standard_name": "air_temperature"}),
    "U": ("eastward_wind", {"units": "m s-1", "standard_name": "eastward_wind"}),
    "V": ("northward_wind", {"units": "m s-1", "standard_name": "northward_wind"}),
    "W": ("upward_wind", {"units": "m s-1", "standard_name": "upward_air_velocity"}),
}
# options only some instruments take, with their defaults
_SPECIFIC = {"stations": None, "sample_rate": "1hz"}


def load_trajectories(
    time_range,
    instrument,
    *,
    stations=None,
    sample_rate="1hz",
    time_handling="fuzzy",
    keep_outliers=False,
    additional_vars=(),
    file_source=False,
    file_path=None,
    verbose=False,
    date_warning=True,
):
    # (the module docstring, so that help(load_trajectories) shows it all)
    options = {k: v for k, v in locals().items() if k not in ("time_range",
                                                              "instrument")}
    call = {"time_range": time_range, "instrument": instrument, **options}
    settings, inst_info, reader = _prepare(time_range, instrument, options)
    flights, file_list = reader.load(settings, inst_info)
    settings["record"] = record(inst_info)
    settings["files_read"] = list(file_list)
    warn_missing_files(settings, stacklevel=2)
    warn_multiple_matches(settings, stacklevel=2)
    # (IAGOS: flights of one day are often L2 3.1.1 and 3.1.2: not warned)
    warn_mixed_versions(settings, stacklevel=2, exempt=("IAGOS",))
    warn_if_nothing_in_range(settings, len(flights), stacklevel=2, what="flights")
    if (settings["scalar_time"] and settings["time_handling"] == "strict"
            and not flights):
        warnings.warn(
            "A single time_range value with time_handling='strict' selects only "
            "flights starting exactly at that instant, and none were found. Use "
            "time_handling 'generous' or 'fuzzy' (default) to get the whole day.",
            stacklevel=2,
        )
    return _to_dataset(flights, file_list, settings, call)


# supporting functions


def _prepare(time_range, instrument, options):
    """Check the arguments; return (settings, inst_info, reader)."""
    check_flags(options, ("keep_outliers", "file_source", "verbose", "date_warning"))
    check_time_handling(options)
    name = find_instrument(instrument, INSTRUMENTS,
                           " (satellite profiles: load_limb_sounders)")
    info = INSTRUMENTS[name]
    reserved = ({"Lat", "Lon", "Time", "Alt", *_OUTPUTS}
                | {out for out, _ in _OUTPUTS.values()} | set(info["per_flight"])
                | set(info["flight_times"])
                | {"latitude", "longitude", "time", "altitude", "start_time",
                   "source_file", "original_files"})
    check_additional_vars(options, reserved)
    def given(option):  # (a value other than the default)
        value = options[option]
        if option == "sample_rate":
            return str(value).lower() != _SPECIFIC[option]
        return value is not None

    for option in _SPECIFIC:  # (options of other instruments)
        if option not in info["options"] and given(option):
            users = [n for n, i in INSTRUMENTS.items() if option in i["options"]]
            raise TypeError(f"{option} applies to {', '.join(users)} only, not {name}.")
    if str(options["sample_rate"]).lower() not in ("1hz", "full"):
        raise ValueError("sample_rate must be '1hz' or 'full'.")
    raw = to_datenums(time_range)
    if options["date_warning"]:
        warn_if_no_data(raw, name, info["time_range"], stacklevel=3)
    expanded, scalar = expand_time_range(raw, options["time_handling"])
    settings = dict(options, instrument=name, time_range=expanded,
                    scalar_time=scalar, sample_rate=str(options["sample_rate"]).lower(),
                    additional_vars=list(options["additional_vars"]))
    inst_info = {"path": "", "inst": name}
    if options["file_path"] is not None:
        inst_info["path"] = str(options["file_path"])
    else:
        inst_info["path"] = str(dataset_dir(name, "raw", must_exist=False))
    if not os.path.isdir(inst_info["path"]):
        given = ("file_path" if options["file_path"] is not None
                 else "the data root and datasets.toml")
        raise FileNotFoundError(f"The {name} data folder {inst_info['path']} does not "
                                f"exist (from {given}; see atmosgw.curate.paths).")
    settings["data_path"] = inst_info["path"]
    return settings, inst_info, info["reader"]


def _to_dataset(flights, file_list, settings, call):
    """The reader's flights as a Dataset in the trajectory format."""
    name = settings["instrument"]
    info = INSTRUMENTS[name]
    rec = settings["record"]
    n_flight = len(flights)
    n_sample = max((np.size(f["Time"]) for f in flights), default=0)

    def grid(key, fill=np.nan, dtype=float):
        out = np.full((n_flight, n_sample), fill, dtype=dtype)
        for k, flight in enumerate(flights):
            if key in flight:
                value = np.asarray(flight[key], dtype=dtype)
                out[k, :value.size] = value
        return out

    present = np.zeros((n_flight, n_sample), dtype=bool)
    for k, flight in enumerate(flights):
        present[k, :np.size(flight["Time"])] = True
    alt = grid("Alt")
    lat, lon = grid("Lat"), grid("Lon")
    with np.errstate(invalid="ignore"):
        bad_geo = (np.abs(lat) > 90) | (lon < -180) | (lon > 180)

    ts = ("flight", "sample")
    variables, flags = {}, {}
    for key in info["variables"]:
        out, attrs = _OUTPUTS[key]
        value = grid(key)
        masks = {"no_observation": ~present, "bad_geolocation": bad_geo}
        removed = bad_geo.copy()
        if key == "Temp":
            masks["outside_physical_range"] = outside_physical_range(value, alt)
            removed |= masks["outside_physical_range"]
        # the provider's own screening, per sample (e.g. IAGOS validity flags)
        provider = grid(f"{key}Bad", fill=False, dtype=bool)
        if provider.any():
            masks["instrument_flag"] = provider
            removed |= provider
        if not settings["keep_outliers"]:
            value = np.where(removed, np.nan, value)
        flag = _state_flags(value, value.shape, masks)
        flags[out] = flag
        variables[out] = (ts, value, {
            **attrs, **record_attrs(key, rec), "origin": "measured",
            "ancillary_variables": f"{out}_flag"})
        variables[f"{out}_flag"] = (ts, flag, states.flag_attrs())
    for extra in settings["additional_vars"]:
        variables[extra] = (ts, grid(extra), {"source_variable": extra})

    for key, long_name in info["per_flight"].items():
        variables[key] = (("flight",), np.array([str(f.get(key, "")) for f in flights],
                                                dtype=str), {"long_name": long_name})
    for key, long_name in info["flight_times"].items():
        variables[key] = (("flight",), datenum_to_datetime(
            np.array([f[key] for f in flights], dtype=float), "ms"),
            {"long_name": long_name})
    coords = {
        "latitude": (ts, lat, {**COORD_ATTRS["latitude"], **record_attrs("Lat", rec)}),
        "longitude": (ts, lon, {**COORD_ATTRS["longitude"],
                                **record_attrs("Lon", rec)}),
        # (to the ms; FAAM at its full 32 Hz to the us)
        "time": (ts, datenum_to_datetime(grid("Time"), "us" if name == "FAAM" and
                                         settings["sample_rate"] == "full" else "ms"), {
            **COORD_ATTRS["time"], **record_attrs("Time", rec)}),
        "altitude": (ts, alt, {**COORD_ATTRS["altitude"], **record_attrs("Alt", rec)}),
        "start_time": (("flight",), datenum_to_datetime(
            np.array([f["start_time"] for f in flights], dtype=float), "ms"),
            {"long_name": "start (launch, take-off) time of the flight"}),
    }
    if settings["file_source"]:
        variables["source_file"] = (("flight",), np.array(
            [f["source_file"] for f in flights], dtype=np.int32),
            {"long_name": "index into original_files (1-based)"})
        coords["original_files"] = (("file",), np.array(file_list, dtype=str))
        coords["original_file_version"] = (("file",), np.array(
            file_versions(rec, file_list), dtype=str),
            {"long_name": "data version stated by each file"})

    spec = FORMATS["trajectory"]["variables"]
    availability = {n: ("not_in_files" if n in info["not_in_files"]
                        else "not_requested")
                    for n, s in spec.items() if s.get("tracked")
                    and n not in variables and n not in coords}
    availability.update({f"{out}_uncertainty": "not_provided" for out in flags})
    common = load_attrs(settings, call, "atmosgw.curate.obs.load_trajectories")
    product = info["product"]
    attrs = {
        "format": "trajectory",
        "format_version": FORMATS["trajectory"]["version"],
        "instrument": name, "product": product,
        "characterisation": f"{name}-{product.upper()}-placeholder/0",
        "data": "raw", "preprocessing": rec.get("preprocessing", ""),
        "time_range": common.pop("time_range"),
        "availability": states.format_availability(availability),
        **common,
        "state_counts": states.counts_text(flags),
    }
    if name == "FAAM":
        attrs["sample_rate"] = settings["sample_rate"]
    return xr.Dataset(variables, coords=coords, attrs=attrs)


load_trajectories.__doc__ = __doc__

__all__ = ["load_trajectories", "INSTRUMENTS"]
