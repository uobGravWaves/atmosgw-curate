"""
Load limb-sounder and other profile data in a common format.

Purpose
-------
Load temperature profiles from several instruments (ACE-FTS, GNSS radio
occultation, HIRDLS, MIPAS, MLS, SABER, SOFIE, the ALOMAR lidar), from
planetary-wave-filtered versions of them, from the output of the Wright and
Hindley (2018) sampler, or from any files in the same MATLAB format ("Misc"),
and return them on one [profile x level] grid, by default interpolated to a
common height scale, with outliers removed.

Port of the MATLAB ``get_limbsounders`` (uobGravWaves/
MultiInstrumentFunctions, commit 78eef78), the group's own code. AIRS and
IASI, which the MATLAB function can also load, are loaded by
``load_nadir_sounders`` instead.

Status of the port: complete for every instrument above. It matches all 33
cases of the MATLAB reference outputs (see Validation).

Inputs and outputs
------------------
``load_limb_sounders(time_range, instrument, **options)`` -> ``xarray.Dataset``

``time_range``: one time (that day) or two (start and end), as
``numpy.datetime64``, ISO strings or ``datetime`` objects. (Numbers are
refused: convert MATLAB datenums with
``atmosgw.curate.matlab_time.datenum_to_datetime``.) Strings must be real
dates: MATLAB-style carry-over such as "2022-01-34" is refused, with the date
it would mean in MATLAB (use ``matlab_time.make_datetime(2022, 1, 34)`` for
that).

``instrument``: ``"ACE"``, ``"GNSS"``, ``"HIRDLS"``, ``"MIPAS"``, ``"MLS"``,
``"SABER"``, ``"SOFIE"``, ``"ALOMAR"``, ``"Sampler"`` or ``"Misc"`` (any case).

Options (type, default):

- ``original_z`` (bool, False): keep the original levels instead of
  interpolating to ``height_scale``.
- ``keep_outliers`` (bool, False): do not remove outliers (see Science).
- ``height_scale`` (array, km, 18:0.5:60): heights to interpolate to.
- ``height_range`` ((min, max), km, (0, inf)): clips ``height_scale``, or the
  data if ``original_z``.
- ``lat_range`` ((min, max), (-90, 90)): keep profiles entering this latitude
  range at any height.
- ``lon_range`` ((min, max), (-180, 180)): the same for longitude.
- ``time_handling`` (str, "fuzzy"): "strict", "generous" or "fuzzy"; see below.
- ``verbose`` (bool, False): print progress.
- ``additional_vars`` (list of str, ()): extra variables to load, by their
  names in the files.
- ``date_warning`` (bool, True): warn if the instrument has no data for the
  dates asked for.
- ``hindley_pws`` (bool, False): load the planetary-wave-filtered data instead
  of the raw data.
- ``ignore_non_modal`` (bool, False): pass through ``additional_vars`` that
  are not per point, profile or level, unprocessed.
- ``file_source`` (bool, False): also return ``source_profile``,
  ``source_file`` and ``original_files``.
- ``dx`` (bool, False): also return the distance from the previous profile at
  each level.
- ``misc_path`` (str or Path, None): folder of the files (Misc only).
- ``misc_identifier`` (str, ""): text in the file names before the date (Misc
  only).
- ``file_path`` (str or Path, None): read the instrument's files from this
  folder instead.
- ``interp_type`` (str, "linear"): interpolation method (see
  ``matlab_compat.interp1``).
- ``extrapolate`` (bool, False): extrapolate beyond each profile's heights
  when interpolating.
- ``sample_path`` (str or Path, None): folder of the sampler output (Sampler
  only).

``time_handling``:

- ``"strict"``: exactly the times given; one time selects only profiles at
  exactly that instant (with a warning if there are none);
- ``"generous"``: whole days, from the start of the first day touched to the
  end of the last;
- ``"fuzzy"``: one time means its whole day; an end time at exactly midnight
  includes that whole day (e.g. 2010-01-01 to 2010-01-02 gives both days).

Profiles are selected by their mean time, so profiles are never split.

**Output**: an ``xarray.Dataset`` with dimensions ``profile`` and ``level``
(and ``mode`` for ``temperature_pw``):

- ``latitude``: degrees north;
- ``longitude``: degrees east, in [-180, 180];
- ``time``: ``datetime64[us]``, to the nearest millisecond;
- ``altitude``: km (the height scale, unless ``original_z``);
- ``pressure``: hPa (estimated from height for some instruments);
- ``temperature``: K (not for PW data);
- ``temperature_residual``, ``temperature_pw`` (with ``mode``): K (PW data);
- ``source_profile``: 1-based profile number in its file (``file_source``);
- ``source_file``: 1-based index into ``original_files`` (``file_source``);
- ``original_files`` (dimension ``file``): file names (``file_source``);
- ``dx``: km from the previous profile at the same level (``dx``).

All but ``original_files`` have dimensions (profile, level).

``additional_vars`` keep their names in the files. A 3-D one has a third
dimension ``<name>_dim2``. A MATLAB table (e.g. GNSS ``MetaData``) becomes
one variable per column, ``<name>_<column>``, along ``profile``. Variables
passed through by ``ignore_non_modal`` have their own dimensions
(``<name>_dim0`` ...), or none if they hold a single value. Missing values
are NaN (NaT for times). With no profiles, ``profile`` has length 0.

The output is in the common "profile" format, version 1
(``docs/design/data_formats.md``; check with ``data.validate(ds,
"profile")``); the values are as before the format was applied. Each data
variable (``temperature``, ``pressure``, ``temperature_residual``,
``temperature_pw``) has a state flag ``<name>_flag`` of the same dimensions
(``data.states``) and an ``origin`` attribute. States set: ``bad_geolocation``
(latitude or longitude out of range), ``outside_physical_range`` (temperature
< 100 K, or > 400 K below 100 km), ``outside_requested_range`` (outside
``height_range``, with ``original_z``; interpolated data are only on the
height scale); as in
MATLAB, a point removed for one variable is removed (and flagged) in all of
them. With ``keep_outliers`` the values are kept and the states still set.
Every other NaN is ``no_observation``: above or below a profile, padding
between profiles of different lengths, or NaN in the file (the loader cannot
tell these apart). With ``extrapolate``, points outside a profile's own
heights are ``extrapolated`` (where they have a value). ``origin``:
"retrieved", except "derived" for planetary-wave data, for Sampler data (from
a model) and for pressure estimated from height, and "mixed" when only some
of the pressures in the output are estimated (files without pressure);
values computed this way also have the ``derived`` state. With no profiles,
pressure counts as "retrieved".
Data variables keep their type (float64 from most readers; the profile format
allows float32 or float64).

Attributes of the Dataset, recording the load: ``format`` ("profile"),
``format_version`` (1), ``product`` ("temperature"), ``characterisation``
(``<instrument>-TEMPERATURE-placeholder/0``, until the characterisation
packages exist), ``availability`` (optional variables of the format not in
the Dataset, and why), ``state_counts`` (points in each state, per
variable), ``instrument``, ``data``
("raw" or "planetary-wave-filtered"), ``heights`` ("interpolated" or
"original"), ``time_range`` (as expanded by ``time_handling``), ``source``,
``atmosgw_curate_version``, ``options`` (the call's arguments as a JSON
string; infinities as "inf"), ``files_read`` (file names, comma-separated),
``files_not_found`` (dates, SABER months or ALOMAR years that had no file)
and ``multiple_matches`` (each day for which more than one file matched,
with the one used; a warning is also given). All are text, so they survive
``to_netcdf``.

Attributes of the variables: ``units`` and ``standard_name`` or
``long_name``; ``source_variable``, the variable's name in the files; and
``comment`` where the values are more than read from the file (e.g.
"estimated from height with atmosphere.height_to_pressure ..." for pressure
from files without one, the HIRDLS geolocation correction, SOFIE's time
units). With interpolation, ``altitude`` is the height scale, and its
comment names the file heights the data were interpolated from.

Files: where more than one file matches a day, the first in name order is
used, as in MATLAB, except for MLS, where the newest data version is used,
as now also in MATLAB (Q8 in ``docs/quirks.md``).

MATLAB -> Python names: ``Lat`` latitude, ``Lon`` longitude, ``Time`` time,
``Temp`` temperature, ``Pres`` pressure, ``Alt`` altitude, ``Temp_Residual``
temperature_residual, ``Temp_PW`` temperature_pw, ``SourceProf``
source_profile, ``SourceFile`` source_file, ``OriginalFiles`` original_files,
``dx`` dx. Options: ``OriginalZ`` original_z, ``KeepOutliers`` keep_outliers,
``HeightScale`` height_scale, ``HeightRange`` height_range, ``LatRange``
lat_range, ``LonRange`` lon_range, ``TimeHandling`` 1/2/3 time_handling
"strict"/"generous"/"fuzzy", ``Verbose`` verbose (default now False),
``AdditionalVars`` additional_vars, ``DateWarning`` date_warning,
``GetHindleyPWs`` hindley_pws, ``IgnoreNonModal`` ignore_non_modal,
``FileSource`` file_source, ``dx`` dx, ``MiscInfo{1}``/``{2}`` misc_path/
misc_identifier, ``FilePath`` file_path, ``InterpType`` interp_type,
``ExtrapData`` extrapolate, ``SamplePath`` sample_path.

Data folders come from ``atmosgw.curate.paths.dataset_dir`` (see
``datasets.toml``).

Science
-------
- Interpolation to the height scale: each profile separately, in height,
  linear by default (``matlab_compat.interp1``); only levels with finite,
  distinct heights are used, and profiles with fewer than three are NaN.
  ``additional_vars`` whose values are all whole numbers (flags, counts) use
  nearest-neighbour instead. Longitudes are unwrapped along each profile
  before interpolating and wrapped back afterwards, so profiles crossing the
  dateline interpolate correctly. Pressure is interpolated in log(pressure)
  [R3] (exact for an isothermal layer in geopotential height; in geometric
  height within about 5e-6 relative for 1 km levels), except with
  ``nearest``, ``next`` or ``previous``; pressures <= 0 become NaN. (Linear
  interpolation, as before 2026-11-25, overestimated it between levels: Q6
  in ``docs/quirks.md``.)
- Outliers (unless ``keep_outliers``) are set to NaN at that point in every
  variable: latitude outside [-90, 90] or longitude outside [-180, 180];
  height outside ``height_scale`` (when interpolating); temperature below
  100 K, or above 400 K below 100 km.
- Pressure, where a file has none, is estimated from height with
  ``atmosphere.height_to_pressure`` (U.S. Standard Atmosphere 1976 [R1]).

Origin
------
Port of ``get_limbsounders.m`` and its modules (uobGravWaves/
MultiInstrumentFunctions, commit 78eef78; C. Wright, with changes by Claude
on 2026-11-22), the group's own code.

Validation
----------
``tests/test_limb_matlab.py``: the private layer against the MATLAB reference
outputs (``get_limbsounders_reference.mat``, made 2026-11-22 with MATLAB
R2024b Update 6 from MultiInstrumentFunctions@78eef78 and
MatlabFunctions@c943399), case by case: all 33 cases (14 synthetic, 19 real
data, for every instrument). 11 cases are bit-identical in every field; the
rest differ only in HIRDLS positions (within 7.1e-14 degrees, from the
geolocation correction's trig functions), in distances derived from them,
and in pressure computed from height (within 4.9e-16 relative).
``tests/test_limb_known.py``: the known-answer checks of the MATLAB test, on
synthetic files and real data. Tolerances are stated in the tests.

AI involvement
--------------
AI-ported (Claude) from the MATLAB. Checked by C. Wright on 2026-11-25 by
comparing its outputs with the MATLAB function's (SABER, HIRDLS); the code
has not been reviewed line by line.

Usage
-----
>>> from atmosgw.curate.obs import load_limb_sounders  # doctest: +SKIP
>>> ds = load_limb_sounders("2007-01-01", "SABER")  # doctest: +SKIP
>>> ds.temperature.isel(level=0)  # doctest: +SKIP

Status
------
Matches the MATLAB reference outputs (see Validation); checked by C. Wright
(see AI involvement). First released in 0.1.0.

Change log
----------
2026/11/24 [Claude] initial version: core, ACE/GNSS/Misc, PW and Sampler readers
2026/11/24 [Claude] added HIRDLS, SABER, MLS, MIPAS, SOFIE and ALOMAR readers
2026/11/25 [Claude] verbose: a progress line for each variable interpolated
2026/11/25 [Claude] source indices not interpolated unless file_source (faster)
2026/11/25 [Claude] status: checked by C. Wright; released in 0.1.0
2026/11/25 [Claude] shared helpers moved to obs/_common.py (outputs unchanged)
2026/11/25 [Claude] output in the profile format: state flags, origins, attributes
2026/11/25 [Claude] states: estimated pressure marked per point; extrapolated
2026/11/25 [Claude] pressure interpolated in log(pressure), not linearly (Q6)
2026/11/25 [Claude] height <-> pressure: U.S. Standard Atmosphere 1976 (Q1-Q3)
2026/11/25 [Claude] MLS: newest version if several match a day (Q8); warn on duplicates
2026/11/25 [Claude] attributes: version, options, files read/not found, variable sources
2026/11/25 [Claude] clearer error for dates that do not exist (e.g. "2022-01-34")
2026/11/25 [Claude] less memory (about -35% peak): fewer copies; outputs unchanged

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import warnings

import numpy as np
import pandas as pd
import xarray as xr

from atmosgw.curate.data import states
from atmosgw.curate.data.formats import FORMATS
from atmosgw.curate.geodesy import haversine
from atmosgw.curate.matlab_compat import unwrap, wrap_to_pi
from atmosgw.curate.matlab_time import datenum, datenum_to_datetime
from atmosgw.curate.obs import (
    _limb_ace_gnss,
    _limb_alomar,
    _limb_hirdls,
    _limb_mipas,
    _limb_mls,
    _limb_pwdata,
    _limb_saber,
    _limb_sampler,
    _limb_sofie,
)
from atmosgw.curate.obs._common import (
    COORD_ATTRS,
    check_additional_vars,
    check_flags,
    check_time_handling,
    expand_time_range,
    find_instrument,
    load_attrs,
    outside_physical_range,
    record,
    to_datenums,
    warn_if_no_data,
)
from atmosgw.curate.obs._limb_tools import ESTIMATED_PRES
from atmosgw.curate.paths import dataset_dir
from atmosgw.curate.profiles import (
    clip_to_height_range,
    interpolate_profiles,
    mask_points,
    profiles_in_box,
    select_profiles,
)

_FAR = datenum(9999, 999, 999)  # "still running" in the MATLAB table

# instruments, with the dates and heights they have data for, and their reader
# (None: not ported yet). Heights are for information only, as in MATLAB.
INSTRUMENTS = {
    "ACE": {"time_range": (datenum(2004, 1, 35), _FAR), "height_range": (18, 125),
            "reader": _limb_ace_gnss},
    "ALOMAR": {"time_range": (datenum(2003, 1, 1), _FAR), "height_range": (0, 80),
               "reader": _limb_alomar},
    "GNSS": {"time_range": (datenum(2000, 1, 1), _FAR), "height_range": (0, 40),
             "reader": _limb_ace_gnss},
    "HIRDLS": {"time_range": (datenum(2005, 1, 29), datenum(2008, 1, 77)),
               "height_range": (0, 80), "reader": _limb_hirdls},
    "MIPAS": {"time_range": (datenum(2002, 7, 2), datenum(2012, 4, 8)),
              "height_range": (0, 80), "reader": _limb_mipas},
    "MLS": {"time_range": (datenum(2004, 1, 275), _FAR), "height_range": (0, 100),
            "reader": _limb_mls},
    "SABER": {"time_range": (datenum(2002, 1, 1), _FAR), "height_range": (0, 120),
              "reader": _limb_saber},
    "SOFIE": {"time_range": (datenum(2007, 1, 135), _FAR), "height_range": (10, 110),
              "reader": _limb_sofie},
    "Sampler": {"time_range": (-_FAR, _FAR), "height_range": (-9999, 9999),
                "reader": _limb_sampler},
    "Misc": {"time_range": (-_FAR, _FAR), "height_range": (-9999, 9999),
             "reader": _limb_ace_gnss},
}

# MATLAB variable names (used inside) -> output names
OUTPUT_NAMES = {
    "Lat": "latitude", "Lon": "longitude", "Time": "time", "Temp": "temperature",
    "Pres": "pressure", "Alt": "altitude", "Temp_Residual": "temperature_residual",
    "Temp_PW": "temperature_pw", "SourceProf": "source_profile",
    "SourceFile": "source_file", "OriginalFiles": "original_files", "dx": "dx",
}

# the data variables of the profile format (each gets a state flag)
_DATA_VARIABLES = ("temperature", "pressure", "temperature_residual",
                   "temperature_pw")
# product name, for the attributes (proposal: one product, temperature
# profiles, for every limb instrument)
_PRODUCT = "temperature"

_ATTRS = {
    **{k: COORD_ATTRS[k] for k in ("latitude", "longitude", "time")},
    "temperature": {"units": "K", "standard_name": "air_temperature"},
    "temperature_residual": {"units": "K",
                             "long_name": "temperature minus planetary waves"},
    "temperature_pw": {"units": "K", "long_name": "planetary-wave temperature"},
    "pressure": {"units": "hPa", "standard_name": "air_pressure"},
    "altitude": COORD_ATTRS["altitude"],
    "source_profile": {"long_name": "profile number in the source file (1-based)"},
    "source_file": {"long_name": "index into original_files (1-based)"},
    "dx": {"units": "km",
           "long_name": "great-circle distance from the previous profile, same level"},
}


def load_limb_sounders(
    time_range,
    instrument,
    *,
    original_z=False,
    keep_outliers=False,
    height_scale=None,
    height_range=(0, np.inf),
    lat_range=(-90, 90),
    lon_range=(-180, 180),
    time_handling="fuzzy",
    verbose=False,
    additional_vars=(),
    date_warning=True,
    hindley_pws=False,
    ignore_non_modal=False,
    file_source=False,
    dx=False,
    misc_path=None,
    misc_identifier="",
    file_path=None,
    interp_type="linear",
    extrapolate=False,
    sample_path=None,
):
    """Load profile data as an xarray.Dataset. See the module docstring."""
    options = {k: v for k, v in locals().items() if k not in ("time_range",
                                                              "instrument")}
    call = {"time_range": time_range, "instrument": instrument, **options}
    data, settings = _load_limb_arrays(time_range, instrument, **options)
    return _to_dataset(data, settings, call)


def _load_limb_arrays(time_range, instrument, **options):
    """
    The MATLAB get_limbsounders, step for step: returns (data, settings).

    data is a dict of NumPy arrays with the MATLAB names (Lat, Temp, ...) and
    times as datenums, as MATLAB's output struct. Used by load_limb_sounders
    and compared directly with the MATLAB reference outputs by the tests.
    """
    settings, inst_info, reader = _prepare(time_range, instrument, options)
    verbose = settings["verbose"]
    height_scale = settings["height_scale"]

    # loading: instrument-specific, see the reader modules
    if verbose:
        print("--> Loading data")
    if settings["hindley_pws"]:
        var_names = ["Lat", "Lon", "Time", "Temp_Residual", "Temp_PW", "Pres", "Alt",
                     "SourceProf", "SourceFile"]
    else:
        var_names = ["Lat", "Lon", "Time", "Temp", "Pres", "Alt", "SourceProf",
                     "SourceFile"]
    # source indices are made only if asked for (MATLAB makes them, interpolates
    # them and drops them at the end; the other results are the same)
    if not settings["file_source"]:
        var_names = [name for name in var_names
                     if name not in ("SourceProf", "SourceFile")]
    var_names = var_names + list(settings["additional_vars"])
    # which pressures the reader estimated from height, point by point (for the
    # output's data states; not in the MATLAB, so taken out of data at the end)
    if "Pres" in var_names and reader in (_limb_ace_gnss, _limb_pwdata):
        var_names = var_names + ["PresEstimated"]
    data, file_list = reader.load(settings, inst_info, var_names)
    settings["record"] = record(inst_info)
    settings["files_read"] = list(file_list)
    multiple = settings["record"]["multiple_matches"]
    if multiple:
        warnings.warn(
            f"{len(multiple)} day(s) had more than one matching "
            f"{settings['instrument']} file; one was used for each (e.g. "
            f"{multiple[0]}). All are listed in "
            "the output attribute 'multiple_matches'.",
            stacklevel=3,
        )

    # variables not on the [profiles x levels] grid of Lat are left untouched
    # (only if ignore_non_modal; otherwise the readers have already refused them)
    skip = []
    if settings["ignore_non_modal"]:
        lat_shape = np.shape(data["Lat"])
        skip = [name for name, value in data.items()
                if not isinstance(value, pd.DataFrame)
                and np.shape(value)[:2] != lat_shape]

    # time range, by the mean time of each profile so profiles are not split
    mean_time = _nanmean_rows(data["Time"])
    keep = (mean_time >= settings["time_range"][0]) & (
        mean_time <= settings["time_range"][1])
    if not keep.all():  # (all kept: nothing to do, and no copy of everything)
        data = select_profiles(data, np.flatnonzero(keep), skip)

    if (settings["scalar_time"] and settings["time_handling"] == "strict"
            and np.shape(data["Time"])[0] == 0):
        warnings.warn(
            "A single time_range value with time_handling='strict' selects only "
            "profiles at exactly that instant, and none were found. Use "
            "time_handling 'generous' or 'fuzzy' (default) to get the whole day.",
            stacklevel=3,
        )

    # (the ALOMAR reader always makes source indices: dropped here if not asked for)
    if not settings["file_source"]:
        data.pop("SourceProf", None)
        data.pop("SourceFile", None)

    if not settings["original_z"]:
        # unwrap longitudes along each profile (finite values only), so that
        # interpolation is safe across the dateline
        lon = np.deg2rad(np.asarray(data["Lon"], dtype=float))
        for i in range(lon.shape[0]):
            finite = np.isfinite(lon[i])
            lon[i, finite] = unwrap(lon[i, finite])

        if verbose:
            print("----> Interpolating to output z")
        extrap = True if settings["extrapolate"] else None
        # heights outside each profile's own heights: extrapolated (data state)
        extrapolated = None
        if settings["extrapolate"] and np.size(data["Alt"]):
            alt_in = np.asarray(data["Alt"], dtype=float)
            with warnings.catch_warnings():  # (all-NaN profiles: no extrapolation)
                warnings.simplefilter("ignore", RuntimeWarning)
                low, high = np.nanmin(alt_in, axis=1), np.nanmax(alt_in, axis=1)
            extrapolated = ((height_scale[None, :] < low[:, None])
                            | (height_scale[None, :] > high[:, None]))
        # (only for progress messages: the variables interpolated below)
        n_interp = sum(name not in ("Alt", "PresEstimated") and name not in skip
                       and not isinstance(data[name], pd.DataFrame)
                       for name in var_names)
        done = 0
        new = {}
        for name in var_names:
            value = lon if name == "Lon" else data[name]
            if name in skip:
                new[name] = value
            elif name == "Alt":
                new["Alt"] = np.tile(height_scale, (np.shape(new["Lat"])[0], 1))
            elif isinstance(value, pd.DataFrame):
                new[name] = value
            else:
                if name != "PresEstimated":
                    done += 1
                if verbose and name != "PresEstimated":
                    print(f"------> {name} ({done} of {n_interp})", flush=True)
                method = settings["interp_type"]
                if name in settings["additional_vars"] and _is_whole_numbers(value):
                    method = "nearest"  # flags and counts keep valid values
                # pressure in log(pressure) [R3], except for methods that only
                # pick existing values; pressures <= 0 have no log: NaN
                log_pres = name == "Pres" and method.lower() not in (
                    "nearest", "next", "previous")
                if log_pres:
                    value = np.asarray(value)  # (keeps float32 as MATLAB)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        value = np.log(np.where(value > 0, value, np.nan))
                new[name] = interpolate_profiles(
                    data["Alt"], value, height_scale, method, extrap
                )
                if log_pres:
                    new[name] = np.exp(new[name])
            # free each input once done (memory); heights are needed to the end
            if name != "Alt":
                data.pop(name, None)
            if name == "Lon":
                lon = None
            value = None
        data = new
        data["Lon"] = np.rad2deg(wrap_to_pi(data["Lon"]))
    else:
        # clip to the height range: NaN outside it, then drop empty levels
        # (the points set to NaN are recorded for the data states, on the
        # levels kept, as clip_to_height_range chooses them)
        alt = np.asarray(data["Alt"])
        with np.errstate(invalid="ignore"):
            outside = ((alt < min(settings["height_range"]))
                       | (alt > max(settings["height_range"])))
        levels = np.flatnonzero(np.any(np.isfinite(alt) & ~outside, axis=0))
        outside_requested = outside[:, levels]
        data = clip_to_height_range(data, data["Alt"], settings["height_range"], skip,
                                    copy=False)  # (in place: memory)

    # longitudes into [-180, 180]
    lon = np.array(data["Lon"], dtype=float)
    lon[lon > 180] = lon[lon > 180] - 360
    lon[lon < -180] = lon[lon < -180] + 360
    data["Lon"] = lon

    if verbose:
        print("--> Postprocessing")

    # why points are bad (data states, recorded with keep_outliers too)
    masks = _find_outliers(data, settings)
    bad = (masks["bad_geolocation"] | masks["outside_requested_range"]
           | masks["outside_physical_range"])
    if settings["original_z"]:  # (already NaN: clipped above)
        masks["outside_requested_range"] = outside_requested
    elif extrapolated is not None:
        masks["extrapolated"] = extrapolated
    if not settings["keep_outliers"]:
        data = mask_points(data, bad, skip, copy=False)  # (in place: memory)
        if verbose:
            print("----> Outliers removed")

    lat_range, lon_range = settings["lat_range"], settings["lon_range"]
    if (min(lat_range) != -90 or max(lat_range) != 90
            or min(lon_range) != -180 or max(lon_range) != 180):
        keep = profiles_in_box(data["Lat"], data["Lon"], lat_range, lon_range)
        data = select_profiles(data, np.flatnonzero(keep), skip)
        masks = select_profiles(masks, np.flatnonzero(keep))
        if verbose:
            print("----> Geographically filtered")

    if settings["file_source"]:
        data["OriginalFiles"] = list(file_list)

    if settings["dx"]:
        lat, lon = np.asarray(data["Lat"]), np.asarray(data["Lon"])
        distance = np.full(lat.shape, np.nan)
        if lat.shape[0] > 1:
            distance[1:] = haversine(lat[:-1], lon[:-1], lat[1:], lon[1:])
        data["dx"] = distance
        if verbose:
            print("----> Pair distances computed")

    settings["state_masks"] = masks  # (for the output's data states)
    estimated = data.pop("PresEstimated", None)
    settings["pres_estimated"] = (None if estimated is None else
                                  np.asarray(estimated, dtype=float) > 0.5)
    return data, settings


# supporting functions


def _prepare(time_range, instrument, options):
    """Check the arguments; return (settings, inst_info, reader)."""
    check_flags(options, ("original_z", "keep_outliers", "verbose", "date_warning",
                          "hindley_pws", "ignore_non_modal", "file_source", "dx",
                          "extrapolate"))
    check_time_handling(options)
    check_additional_vars(options)
    additional = options["additional_vars"]

    # instrument: exact name, any case
    name = find_instrument(instrument, INSTRUMENTS,
                           " (for AIRS and IASI use load_nadir_sounders)")
    info = INSTRUMENTS[name]

    # heights
    if options["height_scale"] is None:
        height_scale = np.arange(18, 60.5, 0.5)
    else:
        height_scale = np.asarray(options["height_scale"], dtype=float).ravel()
    for label, values in (("height_scale", height_scale),
                          ("height_range", options["height_range"])):
        if not np.all(np.asarray(values, dtype=float) >= 0):
            raise ValueError(f"{label} must be >= 0.")
    height_range = tuple(float(h) for h in options["height_range"])
    height_scale = height_scale[(height_scale >= min(height_range))
                                & (height_scale <= max(height_range))]
    lat_range = tuple(float(x) for x in options["lat_range"])
    lon_range = tuple(float(x) for x in options["lon_range"])
    if not all(-90 <= x <= 90 for x in lat_range):
        raise ValueError("lat_range must be within [-90, 90].")
    if not all(-180 <= x <= 180 for x in lon_range):
        raise ValueError("lon_range must be within [-180, 180].")

    # times: datenums inside, so the MATLAB time rules apply exactly
    raw = to_datenums(time_range)
    if np.any(raw < datenum(1979, 1, 1)):
        raise ValueError("time_range must be in the satellite era (1979 or later).")
    if options["date_warning"]:
        warn_if_no_data(raw, name, info["time_range"], stacklevel=3)
    expanded, scalar = expand_time_range(raw, options["time_handling"])

    settings = dict(options, instrument=name, instrument_name=name,
                    time_range=expanded, scalar_time=scalar,
                    height_scale=height_scale, height_range=height_range,
                    lat_range=lat_range, lon_range=lon_range,
                    additional_vars=list(additional))

    # where the files are
    inst_info = {"path": "", "identifier": "", "inst": name}
    reader = info["reader"]
    # (with file_path given, the data root is not needed, so is not looked up)
    if name not in ("Misc", "Sampler") and options["file_path"] is None:
        kind = "pws" if options["hindley_pws"] else "raw"
        try:
            inst_info["path"] = str(dataset_dir(name, kind, must_exist=False))
        except KeyError as error:
            raise ValueError(f"No {kind} data are defined for {name}.") from error
    if options["hindley_pws"]:
        settings["instrument"] = "HindleyPWs"
        reader = _limb_pwdata
    if name == "Misc":
        if options["misc_path"] is None:
            raise ValueError("instrument 'Misc' needs misc_path (and usually "
                             "misc_identifier).")
        inst_info["path"] = str(options["misc_path"])
        inst_info["identifier"] = options["misc_identifier"]
    if options["file_path"] is not None:
        inst_info["path"] = str(options["file_path"])
    if reader is None:
        raise NotImplementedError(f"The {name} reader is not ported yet.")
    return settings, inst_info, reader


def _nanmean_rows(x):
    """Mean of each row ignoring NaN (NaN for rows with no finite values)."""
    x = np.asarray(x, dtype=float)
    if x.ndim < 2:
        return np.full(x.shape[:1], np.nan)
    finite = ~np.isnan(x)
    count = finite.sum(axis=1)
    total = np.where(finite, x, 0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(count > 0, total / np.maximum(count, 1), np.nan)


def _is_whole_numbers(value):
    """True if the finite values are all whole numbers (and there are some)."""
    value = np.asarray(value)
    if value.dtype.kind not in "biuf":
        return False
    finite = value[np.isfinite(value)]
    return finite.size > 0 and bool(np.all(finite == np.round(finite)))


def _find_outliers(data, settings):
    """[profiles x levels] masks of points outside physical or requested
    limits, by data state: bad_geolocation, outside_requested_range (heights
    outside the height scale; with original_z all False here) and
    outside_physical_range (temperature; all False without Temp)."""
    lat, lon, alt = (np.asarray(data[k]) for k in ("Lat", "Lon", "Alt"))
    none = np.zeros(lat.shape, dtype=bool)
    with np.errstate(invalid="ignore"):
        masks = {"bad_geolocation": (lat < -90) | (lat > 90) | (lon < -180)
                 | (lon > 180), "outside_requested_range": none,
                 "outside_physical_range": none}
        if not settings["original_z"]:
            scale = settings["height_scale"]
            masks["outside_requested_range"] = (alt < scale.min()) | (alt > scale.max())
        if "Temp" in data:
            masks["outside_physical_range"] = outside_physical_range(data["Temp"], alt)
    return masks


def _variable_attrs(name, out, settings):
    """Standard attributes of one output variable, plus where it came from."""
    attrs = dict(_ATTRS.get(out, {}))
    rec = settings["record"]
    if name in rec["source_variables"]:
        attrs["source_variable"] = rec["source_variables"][name]
    if name in rec["notes"]:
        attrs["comment"] = rec["notes"][name]
    if name == "Alt" and not settings["original_z"]:
        attrs.pop("source_variable", None)
        origin = rec["notes"].get(name) or rec["source_variables"].get(name)
        attrs["comment"] = "the height scale; the data were interpolated to it" + (
            f" from the file's heights ({origin})" if origin else "")
    return attrs


def _state_flags(value, lat_shape, masks):
    """
    The state flag of one data variable (shape of value): the states of
    state_masks ([profiles x levels] masks by state name, repeated over any
    further dimension), then no_observation at every other NaN (outside a
    profile's heights, padding between profiles of different lengths, or NaN
    in the file: the loader cannot tell these apart).
    """
    def over_value(x):  # [profiles x levels] -> shape of value
        return np.broadcast_to(x.reshape(x.shape + (1,) * (value.ndim - 2)),
                               value.shape)

    flag = states.new_flags(lat_shape)
    for name, mask in masks.items():
        if np.shape(mask) != tuple(lat_shape):  # (a bug, not a data problem)
            raise ValueError(f"state mask {name} has shape {np.shape(mask)}, "
                             f"not that of the data {tuple(lat_shape)}")
        if states.state(name).no_value:
            states.set_state(flag, mask, name)
    flag = over_value(flag).copy()
    # qualified states (e.g. extrapolated) only where there is a value
    for name, mask in masks.items():
        if not states.state(name).no_value:
            states.set_state(flag, over_value(mask) & np.isfinite(value), name)
    states.set_state(flag, np.isnan(value) & ~states.no_value_mask(flag),
                     "no_observation")
    return flag


def _origin(out, value, settings):
    """
    (origin attribute, points to mark derived or None) of a data variable.
    Proposal: limb temperatures and pressures count as retrieved, except
    values computed by the loader or a model (Sampler). Pressure estimated
    from height: "derived" if every value present is, "mixed" if some are
    (those points marked), by the profiles in the output.
    """
    finite = np.isfinite(value)
    if out in ("temperature_residual", "temperature_pw"):
        return "derived", finite  # planetary-wave filtering
    if settings["instrument_name"] == "Sampler":
        return "derived", finite  # sampled from a model
    if out == "pressure":
        if settings["record"]["notes"].get("Pres") == ESTIMATED_PRES:
            return "derived", finite  # every file (ALOMAR)
        estimated = settings.get("pres_estimated")
        if estimated is not None and np.shape(estimated) == value.shape:
            estimated = estimated & finite
            if estimated.any():
                return ("derived" if estimated.sum() == finite.sum() else "mixed",
                        estimated)
    return "retrieved", None


def _to_dataset(data, settings, call=None):
    """Arrange the MATLAB-named arrays as an xarray.Dataset with output names.
    call: the arguments of the public call, recorded in the attributes."""
    lat_shape = np.shape(data["Lat"])
    variables, coords = {}, {}
    for name, value in data.items():
        out = OUTPUT_NAMES.get(name, name)
        if name == "OriginalFiles":
            coords[out] = ("file", np.array(value, dtype=str))
            continue
        if isinstance(value, pd.DataFrame):
            for column in value.columns:
                variables[f"{out}_{column}"] = ("profile", value[column].to_numpy())
            continue
        value = np.asarray(value)
        if name == "Time":
            value = datenum_to_datetime(value, "ms")
        if value.shape[:2] == lat_shape:
            dims = ["profile", "level"]
            dims += ["mode" if out == "temperature_pw" else f"{out}_dim{d}"
                     for d in range(2, value.ndim)]
        elif value.size == 1:
            value, dims = value.reshape(()), []
        else:
            dims = [f"{out}_dim{d}" for d in range(value.ndim)]
        variables[out] = (dims, value, _variable_attrs(name, out, settings))

    # profile format: a state flag and an origin for each data variable
    counts = []
    for out in _DATA_VARIABLES:
        if out not in variables:
            continue
        dims, value, var_attrs = variables[out]
        flag = _state_flags(value, lat_shape, settings.get("state_masks", {}))
        origin, derived = _origin(out, value, settings)
        if derived is not None:
            states.set_state(flag, derived, "derived")
        var_attrs.update(origin=origin, ancillary_variables=f"{out}_flag")
        variables[f"{out}_flag"] = (dims, flag, states.flag_attrs())
        found = states.state_counts(flag)
        counts.append(f"{out}: " + (", ".join(f"{k} {n}" for k, n in found.items())
                                    or "none"))
    tracked = [n for n, s in FORMATS["profile"]["variables"].items()
               if s.get("tracked")]
    availability = {n: "not_requested" for n in tracked
                    if n not in variables and n not in coords}
    availability.update({f"{out}_uncertainty": "not_provided"
                         for out in _DATA_VARIABLES if out in variables})
    instrument = settings["instrument_name"]
    common = load_attrs(settings, call, "atmosgw.curate.obs.load_limb_sounders")
    attrs = {
        "format": "profile", "format_version": FORMATS["profile"]["version"],
        "instrument": instrument, "product": _PRODUCT,
        "characterisation": f"{instrument}-{_PRODUCT.upper()}-placeholder/0",
        "data": "planetary-wave-filtered" if settings["hindley_pws"] else "raw",
        "heights": "original" if settings["original_z"] else "interpolated",
        "time_range": common.pop("time_range"),
        "availability": states.format_availability(availability),
        **common,
        "state_counts": "; ".join(counts),
    }
    return xr.Dataset(variables, coords=coords, attrs=attrs)


__all__ = ["load_limb_sounders", "INSTRUMENTS"]
