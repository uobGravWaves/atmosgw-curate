"""
Remove the background from loader output, leaving gravity-wave perturbations.

``detrend(ds)`` takes the Dataset of any loader and applies the group's
current detrending methods: a Savitzky-Golay filter in height (profiles,
time-height series), a polynomial across the scan (swaths), the Hindley23
planetary-wave files already read by the loader, and an optional box or
Gaussian pre-smoothing. Methods can be chained. Every number is an option,
and the values used are recorded in the output.

Inputs and outputs
------------------
``detrend(ds, method=None, var=None, *, steps=None, skip_nan=False,
mask_tropopause=True, tropopause_km=None, mask_half_width_km=None,
tropopause_options=None, overwrite=False, keep_steps=False, **params)``
-> a new ``xarray.Dataset`` (the input is not altered) with two added
variables, ``<var>_perturbation`` and ``<var>_background``, each with a
``_flag`` variable.

``method`` (str, None): one of the methods below; None chooses by the data:

==============================================  ================================
Data (``format``, ``product``, ``data``)        Default
==============================================  ================================
profile, planetary-wave-filtered                ``"hindley23"``
profile (limb, GNSS ...)                        ``"savgol"``, 18 km, order 2
profile converted from radiosondes              ``"savgol"``, 5 km, order 2
time_height                                     ``"savgol"``, 18 km, order 2
                                                (the limb default; not a
                                                lidar-community method)
swath, 3-D or 2-D                               ``"polynomial"``, order 4,
                                                across track
trajectory                                      refused: convert radiosondes
                                                with ``data.to_profiles`` first
==============================================  ================================

``var`` (str, None): the variable to detrend; None is ``temperature``
(``brightness_temperature`` for 2-D swaths). Any variable on the same grid
can be named (e.g. a wind); only temperature has defaults.

Methods and their options (``**params`` for a single method; defaults in
``DEFAULTS``):

- ``"savgol"``: Savitzky-Golay low-pass along ``level`` as the background.
  ``order`` (2), ``length_km`` (18; 5 for radiosondes), ``min_points``
  (None; ``skip_nan`` only: valid points needed in a frame, default just
  over half the frame). The frame is ``length_km`` divided by the mean
  level spacing, rounded, and raised to the next odd number if even (for a
  descending grid, lowered: MATLAB's arithmetic), so 18 km at 0.5 km
  spacing is 37 levels, 18.5 km. If the spacing varies by 1 % or more
  (``REGULAR_GRID_TOLERANCE``) each profile gets its own frame from its own
  mean spacing; the filter always works in level index, not height.
- ``"polynomial"``: a polynomial fitted along ``dim`` ("cross_track") for
  every other index, as the background. ``order`` (4), ``abscissa``
  ("index": pixel number 1..N as MATLAB; or the name of a 1-D coordinate
  along ``dim``, e.g. "distance_cross"), ``min_points`` (None; ``skip_nan``
  only: default half the row, and at least order + 2).
- ``"hindley23"``: for Datasets loaded with ``hindley_pws=True``: the
  perturbation is ``temperature_residual`` and the background the sum of
  ``temperature_pw`` over modes. ``outlier_removal`` (True, MATLAB's
  ``H23_OutRem``): both are NaN (``outside_physical_range``) where that sum
  is below 100 K, or above 400 K below 100 km. Must be the first step.
- ``"smooth"`` (also ``"boxcar"``, ``"gaussian"``): smooth what the
  earlier steps left; not a background. ``size`` (odd points per dimension,
  as a dict by dimension name, e.g. ``{"cross_track": 3, "along_track":
  3}``, or a tuple in the variable's dimension order; required),
  ``kernel`` ("box" or "gaussian"), ``sigma`` (0.65; Gaussian standard
  deviation as a fraction of the half-width, a number or a dict by
  dimension). Swaths are smoothed within each along-track segment.
  See ``analysis.smooth``.

``steps`` (list, None): a chain, in place of ``method``: names or
``(name, options)`` pairs, e.g. ``[("smooth", {"size": {"cross_track": 3,
"along_track": 3}}), "polynomial"]`` or ``["hindley23", ("savgol",
{"length_km": 18})]``. Each step acts on what the previous one left; the
backgrounds are added. The order matters. After a smoothing step,
perturbation + background is the smoothed input, not the input.

``skip_nan`` (bool, False): how missing values are treated.

- False (MATLAB): a polynomial row with any NaN is all NaN; a
  Savitzky-Golay value is NaN if any value in its frame is NaN (so one NaN
  removes half a frame either side); smoothing spreads NaN over the kernel.
- True: each fit uses the values present in the same row or frame, so only
  the missing values stay missing. Where nothing is missing the result is
  the same as with False (identical for the fits; to rounding, 1e-13 K,
  for smoothing). New behaviour, not in MATLAB.

Tropopause mask (``"savgol"`` steps): a vertical filter turns the kink at
the tropopause into spurious waves, so by default values within
``mask_half_width_km`` (None: half the frame) of each profile's tropopause
are set to NaN in the perturbation and background, with a warning
(``TropopauseMaskWarning``) saying how much was removed. The tropopause is
found in the Dataset's temperature (``find_tropopause``; options through
``tropopause_options``), or given as ``tropopause_km`` (a number, or one
per profile). A tropopause is looked for only in profiles that reach down
to 5 km (``reach_km``): data that start higher (the default 18-60 km limb
grid, lidars) are taken to hold no tropopause and are not masked, because
the winter polar stratosphere and noise otherwise pass for one. If such
data do start below their tropopause, give ``tropopause_km``.
``mask_tropopause=False`` switches the mask off. The heights used are
returned as ``tropopause_altitude``.

``overwrite`` (bool, False): replace an existing ``<var>_perturbation`` or
``<var>_background`` (e.g. the producer's perturbation in AIRS 2-D files);
otherwise that is an error.

``keep_steps`` (bool, False): also return ``<var>_perturbation_step<k>``,
what was left after step k (1-based).

**Output variables** (float32; computed in float64), with the dimensions
and units of ``var``:

- ``<var>_perturbation``: what the chain left; ``<var>_background``: the
  sum of the backgrounds removed. ``origin`` "derived"; the ``detrending``
  attribute is the chain as JSON, with every option and the frame used.
- flags: the input's states (without ``missing``, ``no_observation`` and
  ``outside_requested_range`` where the output has a value), plus ``derived`` on every value, and
  ``not_computed`` where the input has a value but none could be computed
  (a NaN in the frame or row, the tropopause mask).
- ``<var>_perturbation_uncertainty`` and ``<var>_background_uncertainty``
  are listed in ``availability`` as ``not_provided``.

``find_tropopause(altitude, temperature, lapse_rate=2.0, depth_km=2.0,
search_km=(5, 22), resolution_km=0.25, reach_km=None)``: lapse-rate
tropopause height (km) of each profile ([profile x level] arrays, any level
order), NaN where none is found. ``reach_km`` (None: the lower search
limit, 5 km): a profile is searched only if it has a value at or below this
height.

Differences from MATLAB
-----------------------
- ``hindley23``: with ``outlier_removal=False`` the background is NaN where
  every mode is NaN (MATLAB: ``nansum``, so 0). MATLAB's outlier removal
  also removes points with latitude or longitude out of range; here the
  loader has already done that.
- ``savgol`` in a chain subtracts the background from the current residual
  on regular and irregular grids alike (MATLAB's regular branch subtracts
  it from the temperature; the same for a single filter).
- A profile whose heights are all NaN gives NaN (MATLAB stops).
- Outputs are float32.

MATLAB -> Python: ``Tp`` -> ``<var>_perturbation``; ``BG`` ->
``<var>_background``; ``SGOrder``, ``SGLength`` -> ``order``,
``length_km``; ``PreSmooth`` -> a ``"smooth"`` step before ``"polynomial"``.

Science
-------
- Savitzky-Golay: least-squares polynomial of the given order over a
  sliding frame; the first and last half-frames take the fit to the first
  and last full frame (as MATLAB ``sgolayfilt``).
- Tropopause: "the lowest level at which the lapse rate decreases to
  2 degC/km or less, provided that the average lapse-rate, between that
  level and all other higher levels within 2.0 km does not exceed
  2 degC/km" (WMO International Meteorological Vocabulary, 1992, as quoted
  by Wikipedia "Tropopause"; the WMO text itself was not read: UNVERIFIED
  at source). Choices made here, not from the definition: the search is
  limited to 5-22 km; profiles finer than 0.25 km are averaged to 0.25 km
  first, against noise; a profile must reach down to 5 km (``reach_km``)
  and show a lapse rate above the threshold below the level found. (Without
  the first rule, HIRDLS, SABER and ALOMAR data starting at 15-22 km gave
  "tropopauses" at 18.5-22 km in 560 of 5553, 79 of 1412 and 717 of 2429
  profiles, mostly where the winter polar stratosphere cools with height.)
- The mask half-width (half the frame) is the reach of the filter: an
  assumption, not a published value.

Origin
------
The group's own MATLAB (port freely): ``func_filter_sgolay.m`` and
``func_filter_hindley23.m`` (gwanalyse_limb, corwin365/MatlabFunctions-Common
52518fc) and ``airs_4dp_detrend`` in ``prep_airs_3d.m``
(MatlabFunctions-Personal 5112703). The Savitzky-Golay filter itself is
written from its documented behaviour (MathWorks ``sgolayfilt`` was not
copied). ``skip_nan``, chaining, the tropopause mask and ``find_tropopause``
are new code.

Validation
----------
``tests/test_detrend_matlab.py``: Savitzky-Golay, Hindley23 and smoothing
against MATLAB R2024b outputs for HIRDLS and SABER, 2007-01-01; the
polynomial against the ``Tp`` and ``BG`` of ``prep_airs_3d_reference.mat``
(tolerances stated there). ``tests/test_detrend.py``: known answers,
``skip_nan``, chains, the tropopause mask, formats and flags.

Status
------
Matching MATLAB for the ported methods (see Validation); ``skip_nan``, the
tropopause mask and chains are new and preliminary. AI-ported (Claude); to
be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.analysis import detrend  # doctest: +SKIP
>>> out = detrend(load_limb_sounders("2007-01-01", "HIRDLS"))  # doctest: +SKIP
>>> out = detrend(swath, steps=[("smooth", {"size": {"cross_track": 3,
...     "along_track": 3}}), "polynomial"])  # doctest: +SKIP

Change log
----------
2026/11/28 [Claude] initial version
2026/11/28 [Claude] find_tropopause: only profiles reaching down to reach_km (5 km)
2026/11/28 [Claude] hindley23 outlier removal (H23_OutRem); kept outliers keep their states
2026/11/28 [Claude] find_tropopause: all-NaN profile gives NaN (was IndexError)

"""

import json
import warnings

import numpy as np
import scipy.linalg

from atmosgw.curate.analysis.smooth import smooth_array
from atmosgw.curate.data import states

DEFAULTS = {
    "savgol": {"order": 2, "length_km": 18.0, "min_points": None},
    "polynomial": {"order": 4, "dim": "cross_track", "abscissa": "index",
                   "min_points": None},
    "hindley23": {"outlier_removal": True},
    "smooth": {"size": None, "kernel": "box", "sigma": 0.65},
}
# Savitzky-Golay frame for radiosondes (as Wright et al. 2016, AMT 9, 877, are
# understood to have used; not checked against the paper: UNVERIFIED)
SONDE_LENGTH_KM = 5.0
REGULAR_GRID_TOLERANCE = 0.01  # spacing spread below which a grid is regular (MATLAB)
TROPOPAUSE_DEFAULTS = {"lapse_rate": 2.0, "depth_km": 2.0,
                       "search_km": (5.0, 22.0), "resolution_km": 0.25,
                       "reach_km": None}
_ALIASES = {"boxcar": ("smooth", {"kernel": "box"}),
            "gaussian": ("smooth", {"kernel": "gaussian"})}
_BACKGROUND_METHODS = ("savgol", "polynomial", "hindley23")
# "no value" states that cannot stay on a value that exists (the states of
# values removed as bad are kept, as with keep_outliers=True)
_ABSENT_BITS = states.FLAG_DTYPE(sum(s.bit for s in states.STATES if s.no_value
                                     and s.group != "removed as bad"))
_CHUNK_VALUES = 2.0e7  # values per block of sliding frames (memory)


class TropopauseMaskWarning(UserWarning):
    """Values near the tropopause were removed after a vertical filter."""


def detrend(ds, method=None, var=None, *, steps=None, skip_nan=False,
            mask_tropopause=True, tropopause_km=None, mask_half_width_km=None,
            tropopause_options=None, overwrite=False, keep_steps=False, **params):
    """Add ``<var>_perturbation`` and ``<var>_background`` to a copy of ds.
    See the module docstring."""
    kind = ds.attrs.get("format")
    if kind == "trajectory":
        raise ValueError(
            "detrend has no method for data along a path. Radiosondes: put the "
            "ascents on a height grid with data.to_profiles(ds) first. There is "
            "no current method for aircraft.")
    if kind not in ("profile", "swath", "time_height"):
        raise ValueError(f"detrend needs a Dataset from a loader; its 'format' "
                         f"attribute is {kind!r}.")
    chain = _chain(ds, method, steps, params)
    first_is_pw = chain[0][0] == "hindley23"
    if any(name == "hindley23" for name, _ in chain[1:]):
        raise ValueError("'hindley23' must be the first step of a chain.")
    if not any(name in _BACKGROUND_METHODS for name, _ in chain):
        raise ValueError("The chain removes no background: add 'savgol', "
                         "'polynomial' or 'hindley23' (smoothing alone: "
                         "analysis.smooth.smooth_array).")

    # the values the chain starts from (float64 copies)
    if first_is_pw:
        (base, current, flag, attrs, dims, background,
         pw_info) = _planetary_wave_start(ds, var, chain[0][1])
    else:
        base = var or _default_variable(ds)
        if base not in ds:
            hint = (" This Dataset holds planetary-wave-filtered data: use "
                    "method='hindley23' (or make it the first step)."
                    if "temperature_residual" in ds else "")
            raise ValueError(f"Variable {base!r} is not in the Dataset.{hint}")
        dims, attrs = ds[base].dims, ds[base].attrs
        current = np.array(ds[base].values, dtype=np.float64)
        flag = _input_flag(ds, base, current)
        background = np.zeros_like(current)
    names = (f"{base}_perturbation", f"{base}_background")
    present = [n for n in names if n in ds]
    if present and not overwrite:
        raise ValueError(f"{' and '.join(present)} already in the Dataset (e.g. the "
                         "producer's own perturbation); pass overwrite=True to "
                         "replace.")

    # apply the chain
    record, kept, frames_km = [], [], []
    for name, options in chain:
        if name == "hindley23":
            info = pw_info
        elif name == "savgol":
            _need_dim(dims, "level", name)
            removed, info = _savgol_step(_heights(ds, dims), current, dims,
                                         options, skip_nan)
            current, background = current - removed, background + removed
            frames_km.append(info["frame_km"])
        elif name == "polynomial":
            _need_dim(dims, options["dim"], name)
            removed, info = _polynomial_step(ds, current, dims, options, skip_nan)
            current, background = current - removed, background + removed
        else:  # smooth
            current, info = _smooth_step(ds, current, dims, options, skip_nan)
        record.append({"method": name, **_plain(options), **info})
        kept.append(current.copy())

    # the tropopause mask, for vertical filters
    extra = {}
    mask_note = None
    if frames_km and mask_tropopause:
        current, background, extra, mask_note = _mask_tropopause(
            ds, dims, current, background, tropopause_km,
            mask_half_width_km if mask_half_width_km is not None
            else max(frames_km) / 2.0,
            {**TROPOPAUSE_DEFAULTS, **(tropopause_options or {})},
            first_is_pw)

    # assemble the output
    out = ds.copy()
    if overwrite:
        out = out.drop_vars([n for n in out.variables
                             if n in names or n in [f"{m}_flag" for m in names]])
    detrending = json.dumps({"variable": base, "steps": record,
                             "skip_nan": bool(skip_nan),
                             "tropopause_mask": mask_note})
    flags = {}
    for name, values, what in ((names[0], current, "perturbation"),
                               (names[1], background, "background")):
        this_flag = flag.copy()
        # (a value that exists is not "missing": e.g. a Hindley23 background
        # where the residual is missing; "removed as bad" states are kept)
        this_flag[np.isfinite(values)] &= ~_ABSENT_BITS
        states.set_state(this_flag, np.isfinite(values), "derived")
        states.set_state(this_flag,
                         np.isnan(values) & ~states.no_value_mask(flag),
                         "not_computed")
        flags[name] = this_flag
        out[name] = (dims, values.astype(np.float32), {
            "units": attrs.get("units", ""),
            "long_name": f"{base.replace('_', ' ')} {what}",
            "origin": "derived",
            "ancillary_variables": f"{name}_flag",
            "detrending": detrending,
            "comment": f"{what} from analysis.detrend ("
                       + ", ".join(step["method"] for step in record) + ")"})
        out[f"{name}_flag"] = (dims, this_flag, states.flag_attrs())
    if keep_steps:
        for k, values in enumerate(kept, start=1):
            out[f"{names[0]}_step{k}"] = (dims, values.astype(np.float32), {
                "units": attrs.get("units", ""),
                "comment": f"what was left after step {k} "
                           f"({record[k - 1]['method']}) of analysis.detrend; "
                           "before any tropopause mask"})
    for name, (var_dims, values, var_attrs) in extra.items():
        out[name] = (var_dims, values, var_attrs)
    listed = states.parse_availability(out.attrs.get("availability", ""))
    listed.update({f"{n}_uncertainty": "not_provided" for n in names})
    out.attrs["availability"] = states.format_availability(listed)
    counts = [part.strip() for part in str(out.attrs.get("state_counts", "")).split(";")
              if part.strip() and part.split(":")[0].strip() not in names]
    out.attrs["state_counts"] = "; ".join(counts + [states.counts_text(flags)])
    return out


def find_tropopause(altitude, temperature, lapse_rate=2.0, depth_km=2.0,
                    search_km=(5.0, 22.0), resolution_km=0.25, reach_km=None):
    """Lapse-rate tropopause height (km) of each profile; NaN where none.

    altitude (km), temperature (K): [profile x level], levels in any order.
    See the module docstring (Science) for the rule and the choices made.
    """
    altitude = np.atleast_2d(np.asarray(altitude, dtype=np.float64))
    temperature = np.atleast_2d(np.asarray(temperature, dtype=np.float64))
    altitude = np.broadcast_to(altitude, temperature.shape)
    found = np.full(temperature.shape[0], np.nan)
    reach = search_km[0] if reach_km is None else reach_km
    for p in range(temperature.shape[0]):
        z, t = _clean_profile(altitude[p], temperature[p], resolution_km)
        if z.size < 3 or z[0] > reach:  # (no troposphere in the data)
            continue
        lapse = -np.diff(t) / np.diff(z)  # K/km, of the layer above each level
        seen_troposphere = False
        for i in range(z.size - 1):
            if lapse[i] > lapse_rate:
                seen_troposphere = True
                continue
            if (not seen_troposphere or z[i] < search_km[0] or z[i] > search_km[1]
                    or z[-1] < z[i] + depth_km):
                continue
            above = (z > z[i]) & (z <= z[i] + depth_km)
            mean_lapse = -(t[above] - t[i]) / (z[above] - z[i])
            if np.all(mean_lapse <= lapse_rate):
                found[p] = z[i]
                break
    return found


# supporting functions: the chain and its starting values


def _chain(ds, method, steps, params):
    """The chain as [(name, options)], aliases resolved and defaults filled."""
    if steps is not None:
        if method is not None or params:
            raise TypeError("Give either steps, or method with its options, "
                            "not both.")
        raw = [(s, {}) if isinstance(s, str) else (s[0], dict(s[1])) for s in steps]
        if not raw:
            raise ValueError("steps is empty.")
    else:
        raw = [(method if method is not None else _default_method(ds), dict(params))]
    chain = []
    for name, options in raw:
        name = str(name).lower()
        if name in _ALIASES:
            name, fixed = _ALIASES[name]
            options = {**options, **fixed}
        if name not in DEFAULTS:
            raise ValueError(f"Unknown method {name!r}; known: "
                             f"{', '.join(list(DEFAULTS) + list(_ALIASES))}")
        unknown = set(options) - set(DEFAULTS[name])
        if unknown:
            raise TypeError(f"{name!r} has no option(s) {sorted(unknown)}; it has: "
                            f"{', '.join(DEFAULTS[name]) or 'none'}")
        filled = {**DEFAULTS[name], **options}
        if (name == "savgol" and "length_km" not in options
                and ds.attrs.get("converted_from") == "trajectory"):
            filled["length_km"] = SONDE_LENGTH_KM
        chain.append((name, filled))
    return chain


def _default_method(ds):
    """The current method for this kind of data (module docstring table)."""
    kind = ds.attrs.get("format")
    if kind == "swath":
        return "polynomial"
    if kind == "profile" and str(ds.attrs.get("data", "")).lower().startswith(
            "planetary"):
        return "hindley23"
    return "savgol"


def _default_variable(ds):
    """temperature, or brightness_temperature for 2-D swaths."""
    if (ds.attrs.get("format") == "swath"
            and str(ds.attrs.get("product", "")).lower() == "2d"):
        return "brightness_temperature"
    if str(ds.attrs.get("product", "")).lower() == "hlos_wind":
        raise ValueError("No default variable for the 'hlos_wind' product: name "
                         "one (e.g. var='hlos_wind').")
    return "temperature"


def _planetary_wave_start(ds, var, options):
    """Hindley23: the stored residual as the perturbation so far, the sum of
    the stored waves as the background (func_filter_hindley23.m)."""
    if var not in (None, "temperature"):
        raise ValueError("'hindley23' applies to temperature only.")
    if "temperature_residual" not in ds or "temperature_pw" not in ds:
        raise ValueError(
            "'hindley23' needs temperature_residual and temperature_pw: load the "
            "planetary-wave files (load_limb_sounders(..., hindley_pws=True)).")
    current = np.array(ds["temperature_residual"].values, dtype=np.float64)
    waves = np.asarray(ds["temperature_pw"].values, dtype=np.float64)
    some = np.any(np.isfinite(waves), axis=-1)
    background = np.where(some, np.nansum(waves, axis=-1), np.nan)
    flag = _input_flag(ds, "temperature_residual", current)
    removed = 0
    if options["outlier_removal"]:
        # MATLAB H23_OutRem (default true): no perturbation or background where
        # the summed waves are below 100 K, or above 400 K below 100 km (the
        # sum is 0 where every mode is missing, so those points go too)
        total = np.nansum(waves, axis=-1)
        altitude = _heights(ds, ds["temperature_residual"].dims)
        bad = (total < 100.0) | ((total > 400.0) & (altitude < 100.0))
        states.set_state(flag, bad & some & np.isfinite(current),
                         "outside_physical_range")
        removed = int(np.count_nonzero(bad & np.isfinite(current)))
        current[bad] = np.nan
        background = np.where(bad, np.nan, background)
    return ("temperature", current, flag, ds["temperature_residual"].attrs,
            ds["temperature_residual"].dims, background,
            {"values_removed_as_outliers": removed})


def _input_flag(ds, name, values):
    """A copy of the variable's flag; for a variable without one, 'missing'
    at its NaN values."""
    if f"{name}_flag" in ds and ds[f"{name}_flag"].dims == ds[name].dims:
        return np.array(ds[f"{name}_flag"].values, dtype=states.FLAG_DTYPE)
    flag = states.new_flags(values.shape)
    return states.set_state(flag, np.isnan(values), "missing")


def _need_dim(dims, dim, method):
    if dim not in dims:
        raise ValueError(f"{method!r} works along {dim!r}, which this variable "
                         f"(dimensions {dims}) does not have.")


def _plain(options):
    """Options as JSON-ready values."""
    return {k: (list(v) if isinstance(v, tuple) else v) for k, v in options.items()}


def _heights(ds, dims):
    """Altitude (km) broadcast to the variable's shape (float64)."""
    altitude = ds["altitude"]
    shape = tuple(ds.sizes[d] for d in dims)
    order = [d for d in dims if d in altitude.dims]
    values = np.asarray(altitude.transpose(*order).values, dtype=np.float64)
    index = tuple(slice(None) if d in altitude.dims else None for d in dims)
    return np.broadcast_to(values[index], shape)


def _level_last(values, dims):
    """values as [everything else x level], and how to undo it."""
    axis = dims.index("level")
    moved = np.moveaxis(values, axis, -1)
    return moved.reshape(-1, moved.shape[-1]), moved.shape, axis


def _level_back(rows, shape, axis):
    return np.moveaxis(rows.reshape(shape), -1, axis)


# supporting functions: Savitzky-Golay


def _savgol_step(altitude, current, dims, options, skip_nan):
    """The Savitzky-Golay background of current along level, and what was
    used. Frame lengths follow func_filter_sgolay.m."""
    order, length = int(options["order"]), float(options["length_km"])
    alt, shape, axis = _level_last(altitude, dims)
    values = _level_last(current, dims)[0]
    n_levels = values.shape[1]
    with warnings.catch_warnings(), np.errstate(invalid="ignore", divide="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)  # (all-NaN rows)
        steps = np.diff(alt, axis=1)
        spread = np.nanstd(steps, ddof=1) / np.nanmean(steps)
        regular = bool(abs(spread) < REGULAR_GRID_TOLERANCE)
        spacing = (np.full(alt.shape[0], np.nanmean(steps)) if regular
                   else np.nanmean(steps, axis=1))
    frames = _frame_levels(length, spacing)  # one per profile (NaN: no heights)
    background = np.full(values.shape, np.nan)
    for frame in np.unique(frames[np.isfinite(frames)]).astype(int):
        if frame > n_levels or frame <= order:
            raise ValueError(
                f"Savitzky-Golay frame of {frame} levels ({length} km) does not "
                f"suit {n_levels} levels and order {order}: it must be longer "
                "than the order and no longer than the profile.")
        rows = np.flatnonzero(frames == frame)
        background[rows] = _savgol_rows(values[rows], order, frame, skip_nan,
                                        options["min_points"])
    used = frames[np.isfinite(frames)]
    info = {"regular_grid": regular,
            "frame_levels": ([int(used.min()), int(used.max())] if used.size
                             else None),
            "frame_km": (float(used.max() * abs(np.nanmean(steps))) if regular
                         and used.size else length)}
    return _level_back(background, shape, axis), info


def _frame_levels(length_km, spacing_km):
    """Frame length in levels: abs(make_odd(round(length / spacing))), with
    MATLAB's rounding (half away from zero) and make_odd (even: +1)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = length_km / np.asarray(spacing_km, dtype=np.float64)
        frame = np.sign(ratio) * np.floor(np.abs(ratio) + 0.5)
        frame = np.where(np.mod(frame, 2) == 0, frame + 1, frame)
    return np.where(np.isfinite(frame), np.abs(frame), np.nan)


def _savgol_rows(values, order, frame, skip_nan, min_points):
    """Savitzky-Golay smoothing of each row ([row x level]) with one frame.

    Interior points take the centre value of the least-squares polynomial
    over their frame; the first and last half-frames take the fit to the
    first and last full frame. A NaN in a frame gives NaN, unless skip_nan.
    """
    half = (frame - 1) // 2
    offsets = np.arange(-half, half + 1, dtype=np.float64) / max(half, 1)
    design = np.vander(offsets, order + 1, increasing=True)
    project = design @ np.linalg.pinv(design)  # fitted values from frame values
    n_rows, n_levels = values.shape
    out = np.full(values.shape, np.nan)
    block = max(1, int(_CHUNK_VALUES // (n_levels * frame)))
    for start in range(0, n_rows, block):
        part = values[start:start + block]
        frames = np.lib.stride_tricks.sliding_window_view(part, frame, axis=1)
        out[start:start + block, half:n_levels - half] = frames @ project[half]
        out[start:start + block, :half] = part[:, :frame] @ project[:half].T
        out[start:start + block, n_levels - half:] = (part[:, -frame:]
                                                      @ project[half + 1:].T)
    if skip_nan:
        needed = (frame + 1) // 2 if min_points is None else int(min_points)
        needed = max(needed, order + 1)
        for row in np.flatnonzero(np.any(np.isnan(values), axis=1)):
            out[row] = _savgol_gappy(values[row], out[row], order, frame, needed)
    return out


def _savgol_gappy(values, plain, order, frame, needed):
    """One row with gaps: every point whose frame holds a NaN is refitted
    from the frame's valid values (at least `needed` of them); points that
    are NaN stay NaN. plain: the ordinary result, kept where frames are full."""
    n_levels = values.size
    half = (frame - 1) // 2
    valid = np.isfinite(values)
    centre = np.arange(n_levels)
    first = np.clip(centre - half, 0, n_levels - frame)  # the frame MATLAB uses
    index = first[:, None] + np.arange(frame)[None, :]
    present = valid[index]
    redo = np.flatnonzero(valid & ~np.all(present, axis=1))
    out = plain.copy()
    out[~valid] = np.nan
    if redo.size == 0:
        return out
    weight = present[redo].astype(np.float64)
    offsets = (index[redo] - centre[redo, None]) / max(half, 1)
    powers = offsets[:, :, None] ** np.arange(order + 1)[None, None, :]
    filled = np.where(present[redo], values[index[redo]], 0.0)
    normal = np.einsum("nf,nfj,nfk->njk", weight, powers, powers)
    right = np.einsum("nf,nfj->nj", weight * filled, powers)
    enough = weight.sum(axis=1) >= needed
    result = np.full(redo.size, np.nan)
    if np.any(enough):
        coefficients = np.linalg.solve(normal[enough], right[enough][:, :, None])
        result[enough] = coefficients[:, 0, 0]  # the fit at the point itself
    out[redo] = result
    return out


# supporting functions: polynomial


def _polynomial_step(ds, current, dims, options, skip_nan):
    """The polynomial background along options['dim'] (airs_4dp_detrend)."""
    dim, order = options["dim"], int(options["order"])
    axis = dims.index(dim)
    moved = np.moveaxis(current, axis, -1)
    rows = moved.reshape(-1, moved.shape[-1])
    n = rows.shape[1]
    if n <= order:
        raise ValueError(f"A polynomial of order {order} needs more than {order} "
                         f"points along {dim!r}; there are {n}.")
    if options["abscissa"] == "index":
        x = np.arange(1, n + 1, dtype=np.float64)  # pixel number, as MATLAB
    else:
        if options["abscissa"] not in ds or ds[options["abscissa"]].dims != (dim,):
            raise ValueError(f"abscissa {options['abscissa']!r} is not a "
                             f"coordinate along {dim!r} in the Dataset.")
        x = np.asarray(ds[options["abscissa"]].values, dtype=np.float64)
        if not np.all(np.isfinite(x)):
            raise ValueError(f"abscissa {options['abscissa']!r} has missing values.")
    design = x[:, None] ** np.arange(order + 1)[None, :]
    background = np.full(rows.shape, np.nan)
    complete = np.all(np.isfinite(rows), axis=1)
    if np.any(complete):
        coefficients = _least_squares(design, rows[complete].T)
        background[complete] = (design @ coefficients).T
    refitted = 0
    if skip_nan:
        needed = (max(order + 2, (n + 1) // 2) if options["min_points"] is None
                  else max(int(options["min_points"]), order + 1))
        for row in np.flatnonzero(~complete):
            here = np.isfinite(rows[row])
            if here.sum() < needed:
                continue
            coefficients = _least_squares(design[here], rows[row, here])
            background[row, here] = design[here] @ coefficients
            refitted += 1
    info = {"points": int(n), "rows_with_gaps_fitted": int(refitted)}
    return np.moveaxis(background.reshape(moved.shape), -1, axis), info


def _least_squares(design, values):
    """Least-squares coefficients by QR with column pivoting, as MATLAB's
    backslash. (The design matrix of pixel numbers to the 4th power has a
    condition number of 1e8 to 1e9; the SVD route differs from MATLAB by up
    to 2e-7 K in the fitted values, this one by under 2e-12 K.)"""
    return scipy.linalg.lstsq(design, values, lapack_driver="gelsy")[0]


# supporting functions: smoothing


def _smooth_step(ds, current, dims, options, skip_nan):
    """Smooth current (within each along-track segment of a swath)."""
    if options["size"] is None:
        raise TypeError("'smooth' needs size: odd points per dimension, e.g. "
                        "{'cross_track': 3, 'along_track': 3}.")
    size = _by_dimension(options["size"], dims, 1, "size")
    sigma = _by_dimension(options["sigma"], dims, 0.65, "sigma")
    if "profile" in dims and size[dims.index("profile")] != 1:
        raise ValueError("Profiles are not neighbours on a grid: smooth along "
                         "'level' only.")
    kernel = options["kernel"]
    if ("along_track" in dims and "segment" in ds
            and size[dims.index("along_track")] > 1):
        axis = dims.index("along_track")
        segment = np.asarray(ds["segment"].values)
        out = np.full(current.shape, np.nan)
        for value in np.unique(segment):
            index = [slice(None)] * current.ndim
            index[axis] = np.flatnonzero(segment == value)
            out[tuple(index)] = smooth_array(current[tuple(index)], size, kernel,
                                             sigma, skip_nan)
    else:
        out = smooth_array(current, size, kernel, sigma, skip_nan)
    return out, {"size_by_dimension": dict(zip(dims, [int(s) for s in size]))}


def _by_dimension(value, dims, default, name):
    """A per-dimension option as a list in the variable's dimension order."""
    if isinstance(value, dict):
        unknown = set(value) - set(dims)
        if unknown:
            raise ValueError(f"{name}: no dimension(s) {sorted(unknown)} in {dims}")
        return [value.get(d, default) for d in dims]
    if np.ndim(value) == 0:
        return [value] * len(dims)
    if len(value) != len(dims):
        raise ValueError(f"{name} needs one value per dimension {dims}")
    return list(value)


# supporting functions: tropopause


def _mask_tropopause(ds, dims, current, background, tropopause_km, half_width,
                     options, first_is_pw):
    """NaN within half_width (km) of each profile's tropopause, with a warning.
    Returns (current, background, extra variables, a note for the record)."""
    altitude = _heights(ds, dims)
    alt, shape, axis = _level_last(altitude, dims)
    if tropopause_km is not None:
        height = np.broadcast_to(np.asarray(tropopause_km, dtype=np.float64).ravel(),
                                 (alt.shape[0],)).copy()
        source = "given"
    else:
        if "temperature" in ds and ds["temperature"].dims == dims:
            temperature = np.asarray(ds["temperature"].values, dtype=np.float64)
        elif first_is_pw:  # the temperature the planetary-wave fit was made on
            temperature = (np.asarray(ds["temperature_residual"].values, float)
                           + np.nansum(np.asarray(ds["temperature_pw"].values,
                                                  float), axis=-1))
        else:
            warnings.warn(
                "detrend: no temperature on this grid to find the tropopause "
                "in, so nothing was masked near it; pass tropopause_km, or "
                "mask_tropopause=False to silence this.",
                TropopauseMaskWarning, stacklevel=3)
            return current, background, {}, "not applied: no temperature"
        rows = _level_last(temperature, dims)[0]
        height = find_tropopause(alt, rows, **options)
        source = "found"
        reach = (options["search_km"][0] if options["reach_km"] is None
                 else options["reach_km"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # (all-NaN rows)
            lowest = np.nanmin(np.where(np.isfinite(rows), alt, np.nan), axis=1)
        not_searched = int(np.count_nonzero(~(lowest <= reach)))
    near = np.abs(alt - height[:, None]) <= half_width
    values = _level_last(current, dims)[0]
    removed = int(np.count_nonzero(near & np.isfinite(values)))
    affected = int(np.count_nonzero(np.any(near & np.isfinite(values), axis=1)))
    note = {"half_width_km": float(half_width), "tropopause": source,
            "profiles_masked": affected, "values_masked": removed,
            **({"profiles_not_reaching_reach_km": not_searched}
               if source == "found" else {}),
            **({k: (list(v) if isinstance(v, tuple) else v)
                for k, v in options.items()} if source == "found" else {})}
    lead = tuple(d for d in dims if d != "level")
    extra = {"tropopause_altitude": (
        lead, height.reshape(shape[:-1]),
        {"units": "km", "long_name": "tropopause height used by analysis.detrend",
         "comment": (f"lapse-rate tropopause ({source}); NaN: none within the "
                     "data")})}
    if removed:
        total = int(np.count_nonzero(np.isfinite(values)))
        warnings.warn(
            f"detrend: TROPOPAUSE MASK: {removed} values ({100 * removed / total:.1f}"
            f" % of those computed) in {affected} of {alt.shape[0]} profiles lie "
            f"within {half_width:g} km of the tropopause and were set to NaN, "
            "because a vertical filter turns the tropopause into spurious waves. "
            "Narrow it with mask_half_width_km, or keep the values with "
            "mask_tropopause=False (not recommended).",
            TropopauseMaskWarning, stacklevel=3)
        mask = _level_back(near, shape, axis)
        current = np.where(mask, np.nan, current)
        background = np.where(mask, np.nan, background)
    return current, background, extra, note


def _clean_profile(altitude, temperature, resolution_km):
    """One profile's valid points in ascending height, averaged into layers
    of resolution_km if it is finer than that."""
    keep = np.isfinite(altitude) & np.isfinite(temperature)
    z, t = altitude[keep], temperature[keep]
    order = np.argsort(z, kind="stable")
    z, t = z[order], t[order]
    if z.size >= 3 and np.median(np.diff(z)) < resolution_km:
        layer = np.floor((z - z[0]) / resolution_km).astype(int)
        count = np.bincount(layer)
        used = count > 0
        z = np.bincount(layer, weights=z)[used] / count[used]
        t = np.bincount(layer, weights=t)[used] / count[used]
    distinct = np.diff(z, prepend=-np.inf) > 0  # (empty profile: empty)
    return z[distinct], t[distinct]
