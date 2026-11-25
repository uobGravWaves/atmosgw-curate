"""
Operations on swaths stored as [cross-track x along-track (x ...)] arrays.

Purpose
-------
Building blocks for nadir-sounder loaders (AIRS, IASI; later CrIS, IASI-NG):
distances across and along a swath, regridding a swath onto a regular grid in
distance, and putting its levels onto a regular vertical grid. Port of the
regridding in the MATLAB ``prep_airs_3d`` (its ``regularise_airs`` and its
vertical-interpolation step).

Arrays are [cross x along (x level ...)], as MATLAB's ``l1_lat``,
``ret_temp`` (after its permute) and the nadir reader contract
(``docs/design/nadir_loader.md`` section 5); the output Dataset puts
``along_track`` first.

Inputs and outputs
------------------
======================================================  ========================
Function                                                Returns
======================================================  ========================
``arc_distance(lat1, lon1, lat2, lon2)``                ``(arc, azimuth)``:
                                                        great-circle arc in
                                                        degrees, and direction
                                                        (MATLAB ``distance``)
``swath_distances(lat, lon)``                           dict: ``cross`` (km),
                                                        ``along`` (km),
                                                        ``heading`` (degrees),
                                                        ``basis`` (0-based
                                                        middle cross index)
``regular_axis(values, n=None, spacing=None)``          regular 1-D axis over
                                                        the range of values
``regrid(lat, lon, fields, ...)``                       dict (see below)
``regular_levels(z, spacing=3)``                        regular heights, km
``interpolate_levels(z, values, new_z)``                values on new_z
``unwrap_degrees(x)``                                   x with jumps > 180
                                                        degrees removed
======================================================  ========================

``regrid(lat, lon, fields, n_cross=90, n_along=None, dx_cross=None,
dx_along=None, method="linear", extrapolate=None)``:

- ``lat``, ``lon``: [cross x along], degrees.
- ``fields``: dict of arrays [cross x along (x ...)] to regrid (e.g.
  temperature, time as a number). Each trailing slice is regridded
  separately.
- ``n_cross``, ``n_along``: number of grid points across and along track
  (default: as many as the input); ``dx_cross``, ``dx_along``: spacing in km
  instead (override the numbers). As MATLAB's ``NXT``, ``NAT``, ``dXT``,
  ``dAT``.
- ``method``: ``"linear"`` or ``"nearest"``; ``extrapolate``: None (NaN
  outside the input), ``"linear"`` or ``"nearest"``. As MATLAB's
  ``Interpolant`` and ``Extrapolant``.
- Returns a dict: ``lat``, ``lon`` [n_cross x n_along]; ``fields`` (dict, same
  names); ``cross``, ``along`` (km, the regular axes); ``spacing`` ((km
  across, km along); NaN if one point).

Science
-------
As MATLAB, step for step:

1. Distances: on the middle cross-track row (``basis``, MATLAB
   ``round(n_cross/2)``), the cross-track distance of each position from the
   basis point, measured on the first scan (all scans are assumed to have the
   same cross-track geometry), negative on the side of the first position;
   along-track, the cumulative great-circle distance between successive
   basis points, and the heading of each step (the first heading repeated
   for the first point).
2. A regular grid over the range of those distances.
3. Positions of the new grid: along-track, the basis point's latitude,
   longitude and heading interpolated linearly in distance (extrapolated if
   needed); across, a great-circle step of the cross-track distance at
   heading + 90 degrees (``geodesy.reckon``). Longitude and heading are
   unwrapped first (multiples of 360 degrees added at jumps over 180), so a
   swath crossing the dateline, or a heading crossing north, is not
   interpolated across the jump (fixed in the MATLAB on 2026-11-25).
4. Fields: bilinear (or nearest-neighbour) interpolation in distance space,
   with the input's distances as a rectilinear grid (MATLAB
   ``griddedInterpolant``). A NaN spreads to the grid points that use it.

Sphere of radius 6371 km (``geodesy.EARTH_RADIUS_KM``), as MATLAB's defaults.

Vertical (``regular_levels``, ``interpolate_levels``): a grid with the given
spacing over the range of ``z``, starting at the lowest input level (MATLAB
computes it relative to the middle level, but as min:step:max from the
bottom, so despite its comment the grid starts at the lowest level, not the
middle; the same here; for AIRS levels, every 3 km, the two agree);
linear interpolation, NaN outside the input levels (``profiles``
``interpolate_profiles``, which matches MATLAB ``interp1`` exactly).

Origin
------
Port of ``regularise_airs`` and the vertical step of ``prep_airs_3d.m``
(corwin365/MatlabFunctions, commit f1a78c0; C. Wright), the group's own code.
``arc_distance`` is new code matching the documented behaviour of the MATLAB
Mapping Toolbox ``distance`` (spherical default), not a translation of
MathWorks code.

Validation
----------
``tests/test_swath.py``: against MATLAB R2024b output for synthetic swaths
(mid-latitude, crossing the dateline, near the turning latitude) with each
option (``tests/data/matlab_swath_fixtures.mat``); known answers (rotation
in longitude, regular spacing, NaN outside, the vertical grid). Tolerances
stated in the tests.

AI involvement
--------------
AI-ported (Claude) from the MATLAB; to be reviewed by C. Wright.

Usage
-----
>>> import numpy as np
>>> from atmosgw.curate.swath import regular_levels
>>> regular_levels(np.array([0, 3, 6, 9, 12, 15, 20, 25])).tolist()
[0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0]

Status
------
Preliminary.

Change log
----------
2026/11/25 [Claude] initial version

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from atmosgw.curate.geodesy import EARTH_RADIUS_KM, azimuth, km2deg, reckon
from atmosgw.curate.matlab_compat import interp1
from atmosgw.curate.profiles import interpolate_profiles


def arc_distance(lat1, lon1, lat2, lon2):
    """Great-circle arc (degrees) and azimuth from point 1 to point 2 (degrees)."""
    p1, l1, p2, l2 = (np.deg2rad(np.asarray(a, dtype=float))
                      for a in (lat1, lon1, lat2, lon2))
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin((l2 - l1) / 2) ** 2)
    a = np.clip(a, 0, 1)
    arc = np.rad2deg(2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)))
    return arc, azimuth(lat1, lon1, lat2, lon2)


def unwrap_degrees(x):
    """Add multiples of 360 where successive values jump by more than 180."""
    x = np.asarray(x, dtype=float)
    steps = np.diff(x)
    shift = 360.0 * np.cumsum((steps < -180).astype(float) - (steps > 180))
    return x + np.concatenate([[0.0], shift])


def swath_distances(lat, lon):
    """Cross- and along-track distances (km) and headings of a swath."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    basis = _matlab_round(lat.shape[0] / 2) - 1  # MATLAB round(n/2), 1-based

    # cross-track: from the basis point of the first scan, negative on side 1
    arc, _ = arc_distance(lat[basis, 0], lon[basis, 0], lat[:, 0], lon[:, 0])
    arc = np.asarray(arc, dtype=float).copy()
    arc[: basis + 1] = -arc[: basis + 1]
    cross = _deg2km(arc)

    # along-track: cumulative distance and heading between basis points
    seg, head = arc_distance(lat[basis, :-1], lon[basis, :-1],
                             lat[basis, 1:], lon[basis, 1:])
    along = np.concatenate([[0.0], np.cumsum(_deg2km(np.atleast_1d(seg)))])
    head = np.atleast_1d(head)
    heading = np.concatenate([head[:1], head]) if head.size else np.zeros(1)
    return {"cross": cross, "along": along, "heading": heading, "basis": basis}


def regular_axis(values, n=None, spacing=None):
    """Regular axis over [min, max] of values: n points, or a step of spacing."""
    low, high = float(np.nanmin(values)), float(np.nanmax(values))
    if spacing is not None:
        return _matlab_colon(low, abs(float(spacing)), high)
    return _matlab_linspace(low, high, int(n))


def regrid(lat, lon, fields, n_cross=None, n_along=None, dx_cross=None,
           dx_along=None, method="linear", extrapolate=None):
    """Regrid a swath onto a regular grid in distance. See the module docstring."""
    method = method.lower()
    if method not in ("linear", "nearest"):
        raise ValueError("method must be 'linear' or 'nearest'.")
    if extrapolate not in (None, "linear", "nearest"):
        raise ValueError("extrapolate must be None, 'linear' or 'nearest'.")
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    dist = swath_distances(lat, lon)
    basis = dist["basis"]

    new_cross = regular_axis(dist["cross"], n_cross or lat.shape[0], dx_cross)
    new_along = regular_axis(dist["along"], n_along or lat.shape[1], dx_along)

    # positions: along the basis track, then a step across at heading + 90
    def along_track(values):
        return interp1(dist["along"], values, new_along, "linear", True)

    heading = along_track(unwrap_degrees(dist["heading"]))
    lat_mid = along_track(lat[basis, :])
    lon_mid = along_track(unwrap_degrees(lon[basis, :]))
    n_x, n_a = new_cross.size, new_along.size
    new_lat, new_lon = reckon(np.broadcast_to(lat_mid, (n_x, n_a)),
                              np.broadcast_to(lon_mid, (n_x, n_a)),
                              km2deg(np.broadcast_to(new_cross[:, None], (n_x, n_a))),
                              np.broadcast_to(heading + 90, (n_x, n_a)))

    # fields: interpolation in distance space
    query = np.stack(np.meshgrid(new_cross, new_along, indexing="ij"), axis=-1)
    if extrapolate == "nearest":
        query = np.stack([np.clip(query[..., 0], dist["cross"][0], dist["cross"][-1]),
                          np.clip(query[..., 1], dist["along"][0], dist["along"][-1])],
                         axis=-1)
    fill = None if extrapolate == "linear" else np.nan
    new_fields = {}
    for name, value in fields.items():
        value = np.asarray(value, dtype=float)
        flat = value.reshape(value.shape[0], value.shape[1], -1)
        out = np.empty((n_x, n_a, flat.shape[2]))
        for k in range(flat.shape[2]):
            f = RegularGridInterpolator((dist["cross"], dist["along"]), flat[:, :, k],
                                        method=method, bounds_error=False,
                                        fill_value=fill)
            out[:, :, k] = f(query)
        new_fields[name] = out.reshape((n_x, n_a) + value.shape[2:])

    spacing = (abs(new_cross[1] - new_cross[0]) if n_x > 1 else np.nan,
               abs(new_along[1] - new_along[0]) if n_a > 1 else np.nan)
    return {"lat": np.asarray(new_lat), "lon": np.asarray(new_lon),
            "fields": new_fields, "cross": new_cross, "along": new_along,
            "spacing": spacing}


def regular_levels(z, spacing=3.0):
    """Regular heights (km) from the lowest level of z, as MATLAB (see above)."""
    z = np.asarray(z, dtype=float).ravel()
    middle = float(np.median(z))
    if middle not in z:
        middle = float(z[np.argmin(np.abs(z - middle))])
    low, high = z.min() - middle, z.max() - middle
    return _matlab_colon(low, float(spacing), high) + middle


def interpolate_levels(z, values, new_z):
    """Values [..., level] on levels z (km), linearly onto new_z; NaN outside."""
    values = np.asarray(values, dtype=float)
    z = np.asarray(z, dtype=float).ravel()
    flat = values.reshape(-1, values.shape[-1])
    out = interpolate_profiles(np.broadcast_to(z, flat.shape), flat, new_z)
    return out.reshape(values.shape[:-1] + (np.size(new_z),))


# supporting functions


def _deg2km(deg):
    """Degrees of arc to km, as MATLAB deg2km on the default sphere."""
    return np.deg2rad(np.asarray(deg, dtype=float)) * EARTH_RADIUS_KM


def _matlab_round(x):
    """MATLAB round: halves away from zero."""
    return int(np.sign(x) * np.floor(abs(x) + 0.5))


def _matlab_linspace(low, high, n):
    """MATLAB linspace: n points from low to high, both ends exact."""
    if n < 2:
        return np.array([high], dtype=float)
    out = low + np.arange(n) * ((high - low) / (n - 1))
    out[-1] = high
    return out


def _matlab_colon(low, step, high):
    """MATLAB low:step:high for step > 0 (up to high, within rounding).

    The last point is kept if it lies within a few rounding errors of high
    (2 eps relative to the larger end), so a range of 8.99999999999997 steps
    gives 9 points, as in MATLAB. Where the range is within about 1e-13 of a
    whole number of steps, the count can still differ from MATLAB's by one,
    because the distances themselves agree with MATLAB only to about 1e-12 km.
    """
    if step <= 0 or high < low:
        return np.empty(0)
    tol = 2 * np.finfo(float).eps * max(abs(low), abs(high)) / step
    n = int(np.floor((high - low) / step + tol))
    return low + step * np.arange(n + 1)
