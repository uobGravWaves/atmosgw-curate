"""
Distances, directions and positions on a spherical Earth.

Purpose
-------
Great-circle calculations used by the instrument loaders: the distance between
two points, the direction from one to the other, and the point reached by
travelling a given distance in a given direction. MATLAB-compatible versions
of the Mapping Toolbox functions ``azimuth``, ``reckon`` and ``km2deg`` (with
their spherical defaults) and of the group's ``nph_haversine``.

Inputs and outputs
------------------
All angles in degrees. Inputs are arrays (or scalars) broadcast together;
outputs are float64 arrays of the broadcast shape (floats for scalars).

==============================================  =================================
Function                                        Returns
==============================================  =================================
``haversine(lat1, lon1, lat2, lon2)``           great-circle distance, km
``azimuth(lat1, lon1, lat2, lon2)``             direction from point 1 to point
                                                2, degrees clockwise from north,
                                                in [0, 360)
``reckon(lat, lon, arc, az)``                   ``(lat2, lon2)`` reached by
                                                travelling ``arc`` degrees of
                                                great circle in direction ``az``;
                                                ``lon2`` in [-180, 180]
``km2deg(km, radius=6371)``                     distance as degrees of arc
==============================================  =================================

``arc`` may be negative (travel in the opposite direction).

MATLAB behaviour reproduced (checked against MATLAB R2024b output):

- ``azimuth`` from the north pole is 180 and from the south pole 0; towards
  the north pole it is 0 and towards the south pole 180; for identical points
  it is 0;
- ``reckon`` longitudes are wrapped to [-180, 180] as ``wrapToPi`` does
  (positive odd multiples of 180 become 180).

Science
-------
Sphere of radius 6371 km (the default of the MATLAB functions, and the value in
``nph_haversine``). Standard spherical-trigonometry formulas, with phi
latitude, lambda longitude, d the arc length and theta the azimuth:

- haversine: a = sin^2(dphi/2) + cos(phi1) cos(phi2) sin^2(dlambda/2),
  distance = R * 2 atan2(sqrt(a), sqrt(1 - a));
- azimuth: theta = atan2(cos(phi2) sin(dlambda),
  cos(phi1) sin(phi2) - sin(phi1) cos(phi2) cos(dlambda));
- reckon: phi2 = asin(sin(phi1) cos(d) + cos(phi1) sin(d) cos(theta)),
  lambda2 = lambda1 + atan2(sin(theta) sin(d) cos(phi1),
  cos(d) - sin(phi1) sin(phi2)).

Using a sphere rather than an ellipsoid gives errors of up to about 0.5% in
distance; this matches the MATLAB code being ported.

Known limitation, kept for agreement with MATLAB (see ``docs/quirks.md``):

- **Q4** ``haversine`` of points within about 1e-7 degrees of exactly
  antipodal can give NaN: rounding makes a slightly greater than 1, so
  sqrt(1 - a) is undefined (MATLAB stops with an error there instead).

Origin
------
``haversine``: port of ``nph_haversine.m`` by Neil Hindley (University of
Bath; corwin365/MatlabFunctions, commit c943399), the group's own code, ported
with his permission; same radius and order of operations. ``azimuth``,
``reckon``, ``km2deg``: new
code written to match the documented behaviour of the MATLAB Mapping Toolbox
functions and checked against their output; not translations of MathWorks
source code.

Validation
----------
``tests/test_geodesy.py``: against MATLAB R2024b Update 6 for random points,
identical and antipodal points, the dateline, poles, cardinal directions,
negative and zero arc lengths, and matrices; known answers (quarter
circumference, cardinal directions, reckon/azimuth/haversine consistency).
Tolerances stated in the test.

AI involvement
--------------
AI-originated (Claude), to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.geodesy import haversine
>>> round(float(haversine(0, 0, 0, 90)))
10008

Status
------
Preliminary.

Change log
----------
2026/11/24 [Claude] initial version

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import numpy as np

from atmosgw.curate.matlab_compat import wrap_to_pi

EARTH_RADIUS_KM = 6371.0


def haversine(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (port of nph_haversine)."""
    lat1, lon1, lat2, lon2 = (np.deg2rad(np.asarray(a, dtype=float))
                              for a in (lat1, lon1, lat2, lon2))
    delta_lat = lat2 - lat1
    delta_lon = lon2 - lon1
    a = (
        np.sin(delta_lat / 2) ** 2
        + np.cos(lat1) * np.cos(lat2) * np.sin(delta_lon / 2) ** 2
    )
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return _scalar(EARTH_RADIUS_KM * c)


def azimuth(lat1, lon1, lat2, lon2):
    """Direction from point 1 to point 2, degrees clockwise from north, [0, 360)."""
    lat1, lon1, lat2, lon2 = np.broadcast_arrays(
        *(np.asarray(a, dtype=float) for a in (lat1, lon1, lat2, lon2))
    )
    p1, p2 = np.deg2rad(lat1), np.deg2rad(lat2)
    dlon = np.deg2rad(lon2) - np.deg2rad(lon1)
    az = np.arctan2(
        np.cos(p2) * np.sin(dlon),
        np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dlon),
    )
    az = np.mod(np.rad2deg(az), 360)
    az = np.where(az == 360, 0.0, az)
    # pole rules, as MATLAB: towards a pole, then from a pole (which takes priority)
    az = np.where(lat2 == 90, 0.0, np.where(lat2 == -90, 180.0, az))
    az = np.where(lat1 == 90, 180.0, np.where(lat1 == -90, 0.0, az))
    return _scalar(az)


def reckon(lat, lon, arc, az):
    """Point reached from (lat, lon) along a great circle; returns (lat2, lon2)."""
    p1 = np.deg2rad(np.asarray(lat, dtype=float))
    l1 = np.deg2rad(np.asarray(lon, dtype=float))
    d = np.deg2rad(np.asarray(arc, dtype=float))
    th = np.deg2rad(np.asarray(az, dtype=float))
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(th))
    l2 = l1 + np.arctan2(
        np.sin(th) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2)
    )
    return _scalar(np.rad2deg(p2)), _scalar(np.rad2deg(wrap_to_pi(l2)))


def km2deg(km, radius=EARTH_RADIUS_KM):
    """Distance in km as degrees of great-circle arc."""
    return _scalar((np.asarray(km, dtype=float) / radius) * (180 / np.pi))


# supporting functions


def _scalar(x):
    """Return a float for 0-d results, otherwise the array."""
    x = np.asarray(x)
    return float(x) if x.ndim == 0 else x
