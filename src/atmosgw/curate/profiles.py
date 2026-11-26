"""
Operations on sets of profiles stored as [profiles x levels] arrays.

Building blocks for loaders and analyses of profile data (limb sounders,
radio occultation, lidar): join days of data, select profiles, put profiles
on a common height grid, clip to a height range, mask bad points and select
profiles in a latitude/longitude box. They work on a dict of arrays whose
first two dimensions are [profiles x levels] (extra trailing dimensions are
allowed), as produced inside the loaders.

Inputs and outputs
------------------
``data`` is a dict of variables. Each is a NumPy array of shape
[profiles x levels (x ...)], or a ``pandas.DataFrame`` with one row per
profile. ``skip`` is a collection of variable names to leave untouched (e.g.
metadata that is not per profile).

==============================================================  ====================
Function                                                        Returns
==============================================================  ====================
``concat_profiles(days, names, consume=False)``                 dict, profiles of
                                                                all days stacked
``select_profiles(data, keep, skip=())``                        dict, chosen
                                                                profiles only
``interpolate_profiles(alt, values, new_alt, method, extrap)``  array on new_alt
``clip_to_height_range(data, alt, height_range, skip=())``      dict
``mask_points(data, mask, skip=())``                            dict
``profiles_in_box(lat, lon, lat_range, lon_range)``             bool, per profile
==============================================================  ====================

- ``concat_profiles``: ``days`` is a list of dicts; for each name in
  ``names``, the arrays are padded with NaN to the largest size in every
  dimension but the first, converted to float64, and stacked along the first
  dimension (as the MATLAB ``cat_struct``). DataFrames are stacked by row.
  With no days, every variable is a 0 x 0 float64 array. With
  ``consume=True``, each variable is removed from the day dicts once it has
  been stacked, so its memory can be freed during the stacking (for callers
  that do not need the days afterwards); the result is the same.
- ``select_profiles``: ``keep`` is an index array or boolean mask over
  profiles (as the MATLAB ``reduce_struct`` along dimension 1).
- ``interpolate_profiles``: each profile of ``values`` (and each trailing
  slice) is interpolated from its own heights ``alt`` (km, [profiles x
  levels]) to ``new_alt`` (1-D) with ``matlab_compat.interp1``. Only levels
  with a finite height are used, and of repeated heights only the first. A
  profile with fewer than three such levels is all NaN. ``extrap`` is passed
  to ``interp1``: None gives MATLAB's default for the method.
- ``clip_to_height_range``: points with ``alt`` outside ``height_range`` are
  set to NaN in every variable shaped like ``alt`` (all trailing slices), then
  levels with no finite in-range height are removed from every variable.
- ``clip_to_height_range`` and ``mask_points`` take ``copy`` (default
  True). With ``copy=False``, floating-point arrays of the [profiles x
  levels] shape are changed in place (so the input dict's arrays change too)
  instead of copied, which halves the memory needed; the results are the
  same. ``load_limb_sounders`` uses this on its own arrays.
- ``mask_points``: ``mask`` is a boolean [profiles x levels] array; masked
  points become NaN in every variable (all trailing slices). A 2-D variable of
  a different shape is masked through MATLAB-style linear indices (column-
  major), as the MATLAB loader does.
- ``profiles_in_box``: a profile is kept if its latitude and longitude ranges
  overlap the box, i.e. if it could enter the box at any height. Profiles
  with no finite positions are kept.

Origin
------
Ports of steps of the MATLAB ``get_limbsounders`` (uobGravWaves/
MultiInstrumentFunctions, commit 78eef78) and of ``cat_struct`` and
``reduce_struct`` (same repository), the group's own code.

Validation
----------
Through ``tests/test_limb_matlab.py`` (against the MATLAB reference outputs)
and ``tests/test_limb_known.py`` (known answers), and directly in
``tests/test_profiles.py``.

Status
------
Preliminary. AI-ported (Claude) from the MATLAB, to be reviewed by C. Wright.

Usage
-----
>>> import numpy as np
>>> from atmosgw.curate.profiles import interpolate_profiles
>>> alt = np.array([[10.0, 12.0, 14.0, 16.0]])
>>> interpolate_profiles(alt, 2 * alt, np.array([11.0, 15.0]), "linear", None)
array([[22., 30.]])

Change log
----------
2026/11/24 [Claude] initial version

"""

import warnings

import numpy as np
import pandas as pd

from atmosgw.curate.matlab_compat import interp1


def concat_profiles(days, names, consume=False):
    """Stack the days' variables along profiles, padding with NaN."""
    out = {}
    for name in names:
        parts = [day[name] for day in days]
        if not parts:
            out[name] = np.empty((0, 0))
        elif any(isinstance(p, pd.DataFrame) for p in parts):
            out[name] = pd.concat(parts, ignore_index=True)
        else:
            out[name] = _concat_padded([np.asarray(p) for p in parts])
        if consume:
            del parts
            for day in days:
                day.pop(name, None)
    return out


def select_profiles(data, keep, skip=()):
    """Keep only the chosen profiles (index array or boolean mask)."""
    out = {}
    for name, value in data.items():
        if name in skip:
            out[name] = value
        elif isinstance(value, pd.DataFrame):
            out[name] = value.iloc[np.asarray(keep)].reset_index(drop=True)
        else:
            out[name] = np.asarray(value)[keep]
    return out


def interpolate_profiles(alt, values, new_alt, method="linear", extrap=None):
    """Interpolate each profile to new_alt. See the module docstring."""
    alt = np.asarray(alt, dtype=float)
    values = np.asarray(values)
    new_alt = np.asarray(new_alt, dtype=float).ravel()
    if values.ndim > 3:
        raise NotImplementedError("Variables with more than 3 dimensions.")
    n_slices = int(np.prod(values.shape[2:]))  # 1 for 2-D (also when empty)
    as_3d = values.reshape(values.shape[0], values.shape[1], n_slices)
    if method.lower() in _PIECEWISE:
        # all profiles at once; the same arithmetic as interp1, so bit-identical.
        # Levels of values beyond those of alt have no height and are never used
        # (as in the per-profile version and MATLAB; e.g. ALOMAR, Q7)
        levels = as_3d[:, : alt.shape[1]].astype(float)
        out = np.full((values.shape[0], new_alt.size, n_slices), np.nan)
        rows = max(1, int(_CHUNK_ELEMENTS // max(1, alt.shape[1] * new_alt.size)))
        for start in range(0, values.shape[0], rows):
            part = slice(start, start + rows)
            out[part] = _piecewise_rows(alt[part], levels[part], new_alt,
                                        method.lower(), bool(extrap))
    else:
        out = _interpolate_each(alt, as_3d, new_alt, method, extrap)
    return out.reshape((values.shape[0], new_alt.size) + values.shape[2:])


def clip_to_height_range(data, alt, height_range, skip=(), copy=True):
    """NaN outside the height range, then drop levels with nothing in range."""
    alt = np.asarray(alt)
    bad = (alt < min(height_range)) | (alt > max(height_range))
    levels = np.flatnonzero(np.any(np.isfinite(alt) & ~bad, axis=0))
    all_levels = levels.size == bad.shape[1]
    out = {}
    for name, value in data.items():
        if name in skip or isinstance(value, pd.DataFrame):
            out[name] = value
            continue
        value = np.asarray(value)
        if value.shape[:2] == bad.shape:
            value = _set_nan(value, bad, copy)
        out[name] = value if all_levels and not copy else value[:, levels]
    return out


def mask_points(data, mask, skip=(), copy=True):
    """Set masked [profiles x levels] points to NaN in every variable."""
    mask = np.asarray(mask, dtype=bool)
    out = {}
    for name, value in data.items():
        if name in skip or isinstance(value, pd.DataFrame):
            out[name] = value
            continue
        value = np.asarray(value)
        if value.shape[:2] == mask.shape and (value.ndim > 2 or not copy):
            # (2-D of the mask's shape: the same points as the linear indices)
            out[name] = _set_nan(value, mask, copy)
        else:
            out[name] = _set_nan_linear(value, mask)
    return out


def profiles_in_box(lat, lon, lat_range, lon_range):
    """True for profiles that enter the lat/lon box at any height."""
    lat, lon = np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
    with warnings.catch_warnings():
        # all-NaN profiles give NaN quietly, as MATLAB's nanmin/nanmax
        warnings.filterwarnings("ignore", message="All-NaN slice")
        lon_min, lon_max = np.nanmin(lon, axis=1), np.nanmax(lon, axis=1)
        lat_min, lat_max = np.nanmin(lat, axis=1), np.nanmax(lat, axis=1)
    outside = (
        (lon_min > max(lon_range))
        | (lon_max < min(lon_range))
        | (lat_min > max(lat_range))
        | (lat_max < min(lat_range))
    )
    return ~outside


# supporting functions


def _concat_padded(parts):
    """cat_struct along dimension 1: pad with NaN, convert to float64, stack."""
    if any(p.dtype == object or p.dtype.kind in "US" for p in parts):
        return np.concatenate([np.atleast_2d(p) for p in parts], axis=0)
    parts = [np.atleast_2d(p) for p in parts]
    ndim = max(p.ndim for p in parts)
    parts = [p.reshape(p.shape + (1,) * (ndim - p.ndim)) for p in parts]
    other = [max(p.shape[d] for p in parts) for d in range(1, ndim)]
    # filled in place (no padded copy of each part): the same values as padding
    # each part and concatenating
    out = np.full((sum(p.shape[0] for p in parts), *other), np.nan)
    start = 0
    for p in parts:
        out[(slice(start, start + p.shape[0]),)
            + tuple(slice(0, s) for s in p.shape[1:])] = p
        start += p.shape[0]
    return out


_PIECEWISE = ("linear", "nearest", "next", "previous")
_CHUNK_ELEMENTS = 2e7  # profiles x levels x heights per chunk (bounds memory)


def _interpolate_each(alt, values, new_alt, method, extrap):
    """Profile by profile with matlab_compat.interp1 (all methods)."""
    out = np.full((values.shape[0], new_alt.size, values.shape[2]), np.nan)
    for i in range(values.shape[0]):
        good = _usable_levels(alt[i])
        if good.size > 2:
            for k in range(values.shape[2]):
                out[i, :, k] = interp1(
                    alt[i, good], values[i, good, k], new_alt, method, extrap
                )
    return out


def _piecewise_rows(alt, values, q, method, extrapolate):
    """
    Linear, nearest, next or previous interpolation of every profile at once,
    with exactly the arithmetic and rules of matlab_compat.interp1 (so the
    results are identical to _interpolate_each). alt [profiles x levels],
    values [profiles x levels x slices], q 1-D -> [profiles x len(q) x slices].
    """
    n_prof, n_lev, n_slices = values.shape
    out = np.full((n_prof, q.size, n_slices), np.nan)
    if n_prof == 0 or n_lev < 3 or q.size == 0:  # (fewer than 3 levels: NaN)
        return out

    # usable levels (finite heights, first of any repeats), sorted by height
    # and moved to the front of each row; n of them in each profile
    order = np.argsort(alt, axis=1, kind="stable")  # NaN last; repeats in order
    x = np.take_along_axis(alt, order, axis=1)
    v = np.take_along_axis(values, order[:, :, None], axis=1)
    usable = np.isfinite(x)
    usable[:, 1:] &= x[:, 1:] != x[:, :-1]
    front = np.argsort(~usable, axis=1, kind="stable")
    x = np.take_along_axis(x, front, axis=1)
    v = np.take_along_axis(v, front[:, :, None], axis=1)
    n = usable.sum(axis=1)
    x = np.where(np.arange(n_lev) < n[:, None], x, np.inf)
    m = np.maximum(n, 2)[:, None]  # rows with fewer than 3 are set to NaN below

    # interval j holds x[j] <= q < x[j+1]; the last sample goes in the last one
    qq = np.broadcast_to(q, (n_prof, q.size))
    j = (x[:, :, None] <= q[None, None, :]).sum(axis=1) - 1
    j = np.clip(j, 0, m - 2)
    x0 = np.take_along_axis(x, j, axis=1)
    x1 = np.take_along_axis(x, j + 1, axis=1)
    below = qq < x[:, :1]
    above = qq > np.take_along_axis(x, m - 1, axis=1)
    enough = (n > 2)[:, None]

    with np.errstate(invalid="ignore", divide="ignore"):
        for k in range(n_slices):
            vk = v[:, :, k]
            v0 = np.take_along_axis(vk, j, axis=1)
            v1 = np.take_along_axis(vk, j + 1, axis=1)
            first, last = vk[:, :1], np.take_along_axis(vk, m - 1, axis=1)
            if method == "linear":
                t = (qq - x0) / (x1 - x0)
                r = (1 - t) * v0 + t * v1
                r = np.where(v0 == v1, v0, r)
                r = np.where(qq == x0, v0, np.where(qq == x1, v1, r))
            elif method == "nearest":
                r = np.where(qq - x0 < x1 - qq, v0, v1)
                r = np.where(below, first, np.where(above, last, r))
            elif method == "next":
                r = np.where(qq == x0, v0, v1)
                r = np.where(below, first, np.where(above, np.nan, r))
            else:  # previous
                r = np.where(qq == x1, v1, v0)
                r = np.where(below, np.nan, np.where(above, last, r))
            if not extrapolate:
                r = np.where(below | above, np.nan, r)
            r = np.where(np.isnan(qq) | ~enough, np.nan, r)
            out[:, :, k] = r
    return out


def _usable_levels(alt_row):
    """Indices of finite heights, first of any repeats, in index order."""
    _, first = np.unique(alt_row, return_index=True)
    return np.intersect1d(np.flatnonzero(np.isfinite(alt_row)), first)


def _set_nan(value, mask, copy=True):
    """NaN where a [profiles x levels] mask is set, on every trailing slice.
    copy=False changes a floating-point array in place (others are converted
    to float64, i.e. copied, either way)."""
    if value.dtype.kind != "f":
        value = value.astype(float, copy=True)
    elif copy or not value.flags.writeable:
        value = value.copy()
    value[mask] = np.nan
    return value


def _set_nan_linear(value, mask):
    """MATLAB F(find(mask)) = NaN: linear (column-major) indices of the mask."""
    value = value.astype(float, copy=True) if value.dtype.kind != "f" else value.copy()
    index = np.flatnonzero(mask.ravel(order="F"))
    if index.size and index.max() >= value.size:
        raise ValueError("Mask is larger than the variable it is applied to.")
    flat = value.reshape(-1, order="F")
    flat[index] = np.nan
    return flat.reshape(value.shape, order="F")
