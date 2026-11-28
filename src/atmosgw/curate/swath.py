"""
Operations on swaths stored as [cross-track x along-track (x ...)] arrays.

Building blocks for nadir-sounder loaders (AIRS, IASI, CrIS; later IASI-NG):
distances across and along a swath, regridding a swath onto a regular grid in
distance, and putting its levels onto a regular vertical grid. Port of the
regridding in the MATLAB ``prep_airs_3d`` (its ``regularise_airs`` and its
vertical-interpolation step), plus a regridding from footprints at their own
positions (new code, for CrIS).

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
``footprint_track(lat, lon)``                           dict: the track of a
                                                        swath of footprints
``footprint_coordinates(lat, lon)``                     (along, cross) km of
                                                        every footprint
``regrid_footprints(lat, lon, fields, cross_axis,       dict, as ``regrid``
...)``                                                  (see its docstring)
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
- ``flags``: optional data-state flags [cross x along (x ...)] (integers),
  returned in ``flags``: each output point gets the OR of the flags of every
  input point that contributes to it with a non-zero weight.
- ``cross_axis``: a cross-track axis (km) to use instead of one built from
  this swath, so that several segments of a swath share one (overrides
  ``n_cross`` and ``dx_cross``).
- ``cross``: the input's cross-track distance of each scan position (km),
  instead of those of the first scan (MATLAB's assumption).
- ``cross_edge``: cross-track grid points up to this distance (km) beyond
  the input's outermost scan positions take the outermost value instead of
  NaN (returned in ``cross_edge``, a bool per cross-track point); 0 (the
  default, as MATLAB) for none. Only without ``extrapolate``.
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
   ``griddedInterpolant``). A NaN spreads to the grid points that use it
   with a non-zero weight (in MATLAB also to those where its weight is zero,
   e.g. a grid point exactly on a good input point beside a missing one:
   Q14).

Sphere of radius 6371 km (``geodesy.EARTH_RADIUS_KM``), as MATLAB's defaults.

Footprints (``regrid_footprints``; new code, no MATLAB equivalent): for
swaths whose footprints do not lie on straight scan lines. The track is the
point midway between each row's two middle footprints, smoothed over 3 rows
(one CrIS scan: the rows of a scan are unevenly spaced and offset across
track, and an unsmoothed track's heading zigzags with them); each
footprint's cross-track distance and along-track offset from its row's
track point are exact on the sphere (spherical right triangle); values are
interpolated barycentrically in the Delaunay triangles of the footprints in
(along, cross) km, so a field linear in distance is reproduced exactly, and
a NaN vertex with a non-zero weight makes the grid point NaN. Every row
needs positions (a row without them breaks the track's distances).

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
stated in the tests. Footprints: ``tests/test_nadir_cris.py`` (known
distances; a heading that ignores the zigzag within a scan; a linear field
reproduced to 1e-3 K; fields linear in the true along- and cross-track
distance returned as the grid's own distances to 3 km; no bridging of a
missing footprint).

Status
------
Preliminary. AI-ported (Claude) from the MATLAB; to be reviewed by C. Wright.

Usage
-----
>>> import numpy as np
>>> from atmosgw.curate.swath import regular_levels
>>> regular_levels(np.array([0, 3, 6, 9, 12, 15, 20, 25])).tolist()
[0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0]

Change log
----------
2026/11/25 [Claude] initial version
2026/11/27 [Claude] footprint_track, footprint_coordinates, regrid_footprints (CrIS)
2026/11/27 [Claude] review: track smoothed over a scan; vertices option
2026/11/28 [Claude] review: "nearest" is the nearest footprint, not the largest weight

"""

import numpy as np

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
           dx_along=None, method="linear", extrapolate=None, cross_axis=None,
           flags=None, cross=None, cross_edge=0.0):
    """Regrid a swath onto a regular grid in distance. See the module docstring."""
    method = method.lower()
    if method not in ("linear", "nearest"):
        raise ValueError("method must be 'linear' or 'nearest'.")
    if extrapolate not in (None, "linear", "nearest"):
        raise ValueError("extrapolate must be None, 'linear' or 'nearest'.")
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    dist = swath_distances(lat, lon)
    if cross is not None:  # (the input's cross-track distances, instead of scan 1's)
        dist["cross"] = np.asarray(cross, dtype=float)
    basis = dist["basis"]

    if cross_axis is None:
        new_cross = regular_axis(dist["cross"], n_cross or lat.shape[0], dx_cross)
    else:
        new_cross = np.asarray(cross_axis, dtype=float)
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

    # fields: interpolation in distance space. The grid is rectilinear in
    # distance, so cells and weights are found once per axis and applied to
    # every field, level and flag (the arithmetic of scipy's
    # RegularGridInterpolator, which this replaced for speed)
    q_cross, q_along = new_cross, new_along
    # cross-track queries up to cross_edge km beyond the input's edge scan
    # positions take the edge position's value (marked in "cross_edge")
    near = np.zeros(new_cross.size, dtype=bool)
    if cross_edge > 0 and extrapolate is None:
        lo, hi = np.nanmin(dist["cross"]), np.nanmax(dist["cross"])
        near = (((q_cross < lo) & (q_cross >= lo - cross_edge))
                | ((q_cross > hi) & (q_cross <= hi + cross_edge)))
        q_cross = np.where(near, np.clip(q_cross, lo, hi), q_cross)
    if extrapolate == "nearest":
        q_cross = np.clip(q_cross, dist["cross"][0], dist["cross"][-1])
        q_along = np.clip(q_along, dist["along"][0], dist["along"][-1])
    axes = (_axis_weights(dist["cross"], q_cross, method, extrapolate == "linear"),
            _axis_weights(dist["along"], q_along, method, extrapolate == "linear"))
    new_fields = {name: _apply_weights(np.asarray(value, dtype=float), axes)
                  for name, value in fields.items()}
    new_flags = None if flags is None else _apply_flags(np.asarray(flags), axes)

    spacing = (abs(new_cross[1] - new_cross[0]) if n_x > 1 else np.nan,
               abs(new_along[1] - new_along[0]) if n_a > 1 else np.nan)
    return {"lat": np.asarray(new_lat), "lon": np.asarray(new_lon),
            "fields": new_fields, "flags": new_flags, "cross": new_cross,
            "along": new_along, "cross_edge": near,
            "spacing": spacing}


def footprint_track(lat, lon, span=3):
    """
    The track of a swath of footprints [cross x along] (degrees): for each
    along-track row, the point midway (on the sphere) between its two middle
    footprints, smoothed by a centred mean over ``span`` rows (odd; fewer at
    the ends); its cumulative along-track distance (km); and its heading
    (degrees clockwise from north, from the smoothed point ``span`` rows
    before to the one ``span`` rows after, using only points smoothed over
    all ``span`` rows, so fewer at the ends). The default
    span of 3 rows is one CrIS scan of the index layout: the midpoints of a
    scan's 3 rows are 16 km apart but 22 km from the next scan's and offset
    across track by up to a few km, so unsmoothed, the heading from
    neighbouring rows zigzags by about 5 degrees (moving footprints at the
    swath edge up to 90 km along track; found in review, 2019-01-01) and the
    distance summed along the zigzag is too long (0.7% for a 1.5 km offset).
    Returns dict ``lat``, ``lon``, ``along``, ``heading`` (each per row).
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    mid = lat.shape[0] // 2
    xyz = sum(_unit(lat[k], lon[k]) for k in (mid - 1, mid))  # [3 x along]
    n = xyz.shape[1]
    rows = np.arange(n)
    half = np.minimum.reduce([np.full(n, span // 2), rows, n - 1 - rows])
    if span > 1 and n > 2:
        total = np.concatenate([np.zeros((3, 1)), np.cumsum(xyz, axis=1)], axis=1)
        xyz = total[:, rows + half + 1] - total[:, rows - half]
    # (headings between smoothed points over a whole window only, where any)
    full = np.flatnonzero(half == span // 2)
    lo, hi = (full[0], full[-1]) if full.size >= 2 else (0, n - 1)
    t_lat, t_lon = _from_unit(*xyz)
    step, _ = arc_distance(t_lat[:-1], t_lon[:-1], t_lat[1:], t_lon[1:])
    along = np.concatenate([[0.0], np.cumsum(_deg2km(np.atleast_1d(step)))])
    if n > 1:
        before = np.clip(rows - span, lo, hi)
        after = np.clip(rows + span, lo, hi)
        heading = np.atleast_1d(azimuth(t_lat[before], t_lon[before],
                                        t_lat[after], t_lon[after]))
    else:
        heading = np.zeros(n)
    return {"lat": t_lat, "lon": t_lon, "along": along,
            "heading": np.asarray(heading, dtype=float)}


def footprint_coordinates(lat, lon, track=None):
    """
    Along- and cross-track distances (km) of every footprint [cross x along]
    from the track (``footprint_track``): each footprint is measured from
    its own row's track point, exactly on the sphere (spherical right
    triangle with the great circle through the track point at its heading):
    cross = asin(sin d sin A), along offset = atan2(sin d cos A, cos d), with
    d the arc from the track point and A the bearing from the track heading;
    cross is positive to the right of the heading. Returns (along, cross).
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    track = footprint_track(lat, lon) if track is None else track
    arc, bearing = arc_distance(track["lat"][None, :], track["lon"][None, :], lat, lon)
    d = np.deg2rad(np.asarray(arc, dtype=float))
    a = np.deg2rad(np.asarray(bearing, dtype=float) - track["heading"][None, :])
    cross = np.arcsin(np.clip(np.sin(d) * np.sin(a), -1, 1)) * EARTH_RADIUS_KM
    offset = np.arctan2(np.sin(d) * np.cos(a), np.cos(d)) * EARTH_RADIUS_KM
    return track["along"][None, :] + offset, cross


def regrid_footprints(lat, lon, fields, cross_axis, n_along=None, dx_along=None,
                      method="linear", extrapolate=False, flags=None, edge=0.0,
                      max_edge_km=200.0, block=135, vertices=None):
    """
    Regrid footprints at their own positions onto a regular grid in distance,
    for instruments whose footprints do not lie on straight scan lines (CrIS:
    each field of regard's 3 x 3 fields of view rotate towards the swath
    edges).

    - ``lat``, ``lon``: [cross x along] footprint positions (degrees), in any
      layout whose rows follow the track (CrIS: the index layout).
    - ``fields``: dict of [cross x along (x ...)] arrays to regrid.
    - ``cross_axis``: the regular cross-track axis (km, positive right of the
      track); the along-track axis is regular over the track's distances,
      ``n_along`` points or a step of ``dx_along`` km.
    - ``method``: "linear": barycentric interpolation in the Delaunay
      triangles of the footprints in (along, cross) km, which reproduces a
      field linear in distance exactly; "nearest": the nearest footprint in
      (along, cross) km, for grid points inside the triangles.
    - A NaN at a vertex with a non-zero weight makes the grid point NaN (no
      bridging of missing values); ``flags``: each grid point gets the OR of
      the flags of those vertices.
    - Grid points outside the triangles (outside the swath, or in a triangle
      with an edge longer than ``max_edge_km``) are NaN with no flag, unless
      the nearest footprint is within ``edge`` km, or ``extrapolate``: then
      they take its value and flags, and are marked in ``extrapolated``.
    - Triangulated in blocks of ``block`` output rows, each with the
      footprints within 300 km along track of it.
    - ``vertices``: [cross x along] bool, the footprints that may be
      vertices (default: all with positions). Leaving out footprints that
      never hold a value (CrIS fields of regard 1 and 30) lets grid points
      beside them take the edge rule instead of a NaN vertex.

    Returns a dict as ``regrid``: ``lat``, ``lon`` (the grid's own positions,
    from the track and a step across at heading + 90 degrees), ``fields``,
    ``flags``, ``cross``, ``along``, ``spacing``, and ``extrapolated``
    [cross x along] (bool), ``track`` (``footprint_track``).
    """
    from scipy.spatial import Delaunay, cKDTree

    method = method.lower()
    if method not in ("linear", "nearest"):
        raise ValueError("method must be 'linear' or 'nearest'.")
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    track = footprint_track(lat, lon)
    f_along, f_cross = footprint_coordinates(lat, lon, track)
    new_cross = np.asarray(cross_axis, dtype=float)
    new_along = regular_axis(track["along"], n_along or lat.shape[1], dx_along)
    n_x, n_a = new_cross.size, new_along.size

    # the grid's own positions: along the track, then a step across
    def along_track(values):
        return interp1(track["along"], values, new_along, "linear", True)

    heading = along_track(unwrap_degrees(track["heading"]))
    lat_mid = along_track(track["lat"])
    lon_mid = along_track(unwrap_degrees(track["lon"]))
    grid_lat, grid_lon = reckon(np.broadcast_to(lat_mid, (n_x, n_a)),
                                np.broadcast_to(lon_mid, (n_x, n_a)),
                                km2deg(np.broadcast_to(new_cross[:, None], (n_x, n_a))),
                                np.broadcast_to(heading + 90, (n_x, n_a)))

    flat = {name: np.asarray(v, dtype=float).reshape(lat.size, -1)
            for name, v in fields.items()}
    shapes = {name: np.shape(v)[2:] for name, v in fields.items()}
    flat_flags = None if flags is None else np.asarray(flags).reshape(lat.size, -1)
    points = np.column_stack([f_along.ravel(), f_cross.ravel()])
    usable = np.all(np.isfinite(points), axis=1)
    if vertices is not None:
        usable &= np.asarray(vertices, dtype=bool).ravel()
    out = {name: np.full((n_x * n_a, flat[name].shape[1]), np.nan) for name in flat}
    out_flags = (None if flat_flags is None else
                 np.zeros((n_x * n_a, flat_flags.shape[1]), dtype=flat_flags.dtype))
    extrapolated = np.zeros(n_x * n_a, dtype=bool)
    q_along = np.broadcast_to(new_along[None, :], (n_x, n_a)).ravel()
    q_cross = np.broadcast_to(new_cross[:, None], (n_x, n_a)).ravel()
    margin = 300.0  # km along track: more than a field of regard and its triangles
    for start in range(0, n_a, block):
        cols = np.arange(start, min(start + block, n_a))
        a0, a1 = new_along[cols[0]], new_along[cols[-1]]
        near = np.flatnonzero(usable & (points[:, 0] >= a0 - margin)
                              & (points[:, 0] <= a1 + margin))
        query = (np.arange(n_x)[:, None] * n_a + cols[None, :]).ravel()
        q = np.column_stack([q_along[query], q_cross[query]])
        if near.size < 3:
            continue
        tri = Delaunay(points[near])
        tree = cKDTree(points[near])
        simplex = tri.find_simplex(q)
        inside = simplex >= 0
        vertices = np.zeros((q.shape[0], 3), dtype=int)
        weights = np.zeros((q.shape[0], 3))
        if inside.any():
            s = simplex[inside]
            transform = tri.transform[s]
            b = np.einsum("ijk,ik->ij", transform[:, :2], q[inside] - transform[:, 2])
            weights[inside] = np.column_stack([b, 1 - b.sum(axis=1)])
            vertices[inside] = tri.simplices[s]
            # (long-edged triangles bridge more than the footprints' spacing)
            corners = points[near][vertices[inside]]  # [n x 3 x 2]
            sides = np.linalg.norm(corners - np.roll(corners, 1, axis=1), axis=2)
            long_edge = np.max(sides, axis=1) > max_edge_km
            keep = np.flatnonzero(inside)
            inside[keep[long_edge]] = False
        weights[~inside] = 0.0
        if method == "nearest":
            # (the nearest footprint in (along, cross) km, inside the same
            # support; not the triangle's largest-weight vertex, which can be
            # a farther one: found in review, up to 16 km)
            weights[:] = 0.0
            vertices[:] = 0
            if inside.any():
                vertices[inside, 0] = tree.query(q[inside])[1]
                weights[inside, 0] = 1.0
        used = weights > 1e-12
        source = near[vertices]  # [n x 3] footprint indices
        for name, value in flat.items():
            v = value[source]  # [n x 3 x k]
            bad = np.any(np.isnan(v) & used[:, :, None], axis=1)
            total = np.sum(np.where(used[:, :, None], weights[:, :, None] * v, 0.0),
                           axis=1)
            total[bad | ~inside[:, None]] = np.nan
            out[name][query] = total
        if out_flags is not None:
            fv = flat_flags[source]  # [n x 3 x k]
            out_flags[query] = np.bitwise_or.reduce(
                np.where(used[:, :, None], fv, 0), axis=1)
        # outside the triangles: the nearest footprint, if close enough
        outside = np.flatnonzero(~inside)
        if outside.size and (extrapolate or edge > 0):
            dist, nearest = tree.query(q[outside])
            take = outside[dist <= edge] if not extrapolate else outside
            pick = near[nearest[dist <= edge]] if not extrapolate else near[nearest]
            for name, value in flat.items():
                out[name][query[take]] = value[pick]
            if out_flags is not None:
                out_flags[query[take]] = flat_flags[pick]
            extrapolated[query[take]] = True

    new_fields = {name: out[name].reshape((n_x, n_a) + shapes[name])
                  for name in flat}
    new_flags = (None if out_flags is None else
                 out_flags.reshape((n_x, n_a) + np.shape(flags)[2:]))
    spacing = (abs(new_cross[1] - new_cross[0]) if n_x > 1 else np.nan,
               abs(new_along[1] - new_along[0]) if n_a > 1 else np.nan)
    return {"lat": np.asarray(grid_lat), "lon": np.asarray(grid_lon),
            "fields": new_fields, "flags": new_flags, "cross": new_cross,
            "along": new_along, "spacing": spacing,
            "extrapolated": extrapolated.reshape(n_x, n_a), "track": track}


def _unit(lat, lon):
    """(x, y, z) of points on the unit sphere."""
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def _from_unit(x, y, z):
    """(lat, lon) degrees of the directions (x, y, z) (need not be unit)."""
    return (np.rad2deg(np.arctan2(z, np.hypot(x, y))), np.rad2deg(np.arctan2(y, x)))


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
    plan = _level_plan(z, new_z)
    if plan is None:  # heights not finite and distinct: the general route
        z = np.asarray(z, dtype=float).ravel()
        flat = values.reshape(-1, values.shape[-1])
        out = interpolate_profiles(np.broadcast_to(z, flat.shape), flat, new_z)
        return out.reshape(values.shape[:-1] + (np.size(new_z),))
    j, t, q, x0, x1, outside, enough = plan[:7]
    flat = values.reshape(-1, values.shape[-1])
    if len(plan) > 7:
        flat = flat[:, plan[7]]
    out = np.empty((flat.shape[0], q.size))
    block = 200_000  # columns at a time, to bound temporary memory
    for start in range(0, flat.shape[0], block):
        v = flat[start:start + block]
        v0, v1 = v[:, j], v[:, j + 1]
        with np.errstate(invalid="ignore"):
            r = (1 - t) * v0 + t * v1
        r = np.where(v0 == v1, v0, r)
        r = np.where(q == x0, v0, np.where(q == x1, v1, r))
        r[:, outside] = np.nan
        if not enough:
            r[...] = np.nan
        out[start:start + block] = r
    return out.reshape(values.shape[:-1] + (q.size,))


def interpolate_level_flags(z, flags, new_z):
    """
    Data-state flags [..., level] carried to new_z with interpolate_levels:
    each new level takes the flags of the levels that contribute to it with a
    non-zero weight (OR); none outside the input levels.
    """
    flags = np.asarray(flags)
    plan = _level_plan(z, new_z)
    if plan is None:
        raise ValueError("Flags need finite, distinct heights.")
    j, t, q, x0, x1, outside, enough = plan[:7]
    flat = flags.reshape(-1, flags.shape[-1])
    if len(plan) > 7:
        flat = flat[:, plan[7]]
    use0 = (t != 1) & ~outside
    use1 = (t != 0) & ~outside
    out = (np.where(use0, flat[:, j], 0) | np.where(use1, flat[:, j + 1], 0))
    return out.astype(flags.dtype).reshape(flags.shape[:-1] + (q.size,))


def _level_plan(z, new_z):
    """Intervals and weights for interpolating from levels z to new_z, with the
    arithmetic and rules of matlab_compat.interp1 (linear, no extrapolation),
    once for all columns. None if the heights are not finite and distinct."""
    z = np.asarray(z, dtype=float).ravel()
    q = np.asarray(new_z, dtype=float).ravel()
    if not np.all(np.isfinite(z)) or np.unique(z).size != z.size:
        return None
    order = np.argsort(z, kind="stable")
    x = z[order]
    n = x.size
    m = max(n, 2)
    j = np.clip((x[:, None] <= q[None, :]).sum(axis=0) - 1, 0, m - 2)
    x0, x1 = x[j], x[j + 1]
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (q - x0) / (x1 - x0)
    outside = (q < x[0]) | (q > x[-1]) | np.isnan(q)
    plan = (j, t, q, x0, x1, outside, n > 2)
    if not np.array_equal(order, np.arange(n)):
        plan = plan + (order,)
    return plan


def _axis_weights(grid, q, method, extrapolate):
    """Cell index and weights along one axis, as RegularGridInterpolator."""
    i = np.clip(np.searchsorted(grid, q) - 1, 0, grid.size - 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        y = (q - grid[i]) / (grid[i + 1] - grid[i])
    outside = (q < grid[0]) | (q > grid[-1])
    if method == "nearest":
        i = np.where(y <= 0.5, i, i + 1)
        w = [np.ones_like(y)]
        idx = [i]
    else:
        w = [1 - y, y]
        idx = [i, i + 1]
    if not extrapolate:
        w = [np.where(outside, np.nan, wk) for wk in w]
    return idx, w


def _apply_weights(value, axes, block=1024):
    """Interpolate value [cross x along (x ...)] with per-axis cells and weights,
    in blocks of output along-track rows (to bound temporary memory)."""
    (ix, wx), (iy, wy) = axes
    extra = (None,) * (value.ndim - 2)
    n_a = iy[0].size
    out = np.empty((ix[0].size, n_a) + value.shape[2:])
    for start in range(0, n_a, block):
        part = slice(start, start + block)
        acc = None
        for a, wxa in zip(ix, wx, strict=True):
            for b, wyb in zip(iy, wy, strict=True):
                weight = (wxa[:, None] * wyb[None, part])[(...,) + extra]
                # a corner with zero weight adds nothing, even if it is NaN
                # (as the flags, _apply_flags; MATLAB's griddedInterpolant
                # spreads the NaN: Q14)
                with np.errstate(invalid="ignore"):
                    term = np.where(weight == 0, 0.0,
                                    value[np.ix_(a, b[part])] * weight)
                acc = term if acc is None else acc + term
        out[:, part] = acc
    return out


def _apply_flags(flags, axes, block=1024):
    """OR of the flags of every input point with a non-zero weight (in blocks)."""
    (ix, wx), (iy, wy) = axes
    extra = (None,) * (flags.ndim - 2)
    n_a = iy[0].size
    out = np.zeros((ix[0].size, n_a) + flags.shape[2:], dtype=flags.dtype)
    for start in range(0, n_a, block):
        part = slice(start, start + block)
        for a, wxa in zip(ix, wx, strict=True):
            for b, wyb in zip(iy, wy, strict=True):
                used = ((np.nan_to_num(wxa[:, None]) != 0)
                        & (np.nan_to_num(wyb[None, part]) != 0))
                picked = flags[np.ix_(a, b[part])]
                out[:, part] |= np.where(used[(...,) + extra], picked, 0).astype(
                    flags.dtype)
    return out


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
