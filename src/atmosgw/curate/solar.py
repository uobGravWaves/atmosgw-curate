"""
Solar zenith angle, and the AIRS day/night retrieval switch.

Purpose
-------
The solar zenith angle at given times and places, and whether the AIRS
three-dimensional temperature retrieval (Hoffmann and Alexander) used its
day or night set-up there. Port of ``which_airs_retrieval.m``.

Inputs and outputs
------------------
======================================================  ========================
Function                                                Returns
======================================================  ========================
``solar_zenith_angle(time, lat, lon)``                  degrees, [0, 180]
``day_retrieval(time, lat, lon, threshold=96)``         float array: 1 day,
                                                        0 night, NaN where the
                                                        angle is NaN
======================================================  ========================

``time``: ``datetime64`` (any unit; UTC), or numbers as MATLAB datenums with
``datenum=True``. ``lat``, ``lon``: degrees. All broadcast together.

Science
-------
Low-precision solar coordinates, as in L. Hoffmann's C function ``sza``
(quoted in ``which_airs_retrieval.m``): with D the days since
2000-01-01T12:00 UTC, the mean anomaly g = 357.529 + 0.98560028 D degrees,
mean longitude q = 280.459 + 0.98564736 D degrees, ecliptic longitude
L = q + 1.915 sin g + 0.020 sin 2g, obliquity e = 23.439 - 3.6e-7 D;
declination asin(sin e sin L); right ascension atan2(cos e sin L, cos L);
Greenwich mean sidereal time 18.697374558 + 24.06570982441908 D hours; hour
angle from the local sidereal time; then
cos(SZA) = sin(lat) sin(dec) + cos(lat) cos(dec) cos(h).
These are the U.S. Naval Observatory's approximate solar coordinates
(UNVERIFIED: citation not yet checked against the source); accuracy about
0.01-0.1 degrees, ample for the day/night switch.

Day/night: the retrieval uses its day set-up where the angle is below 96
degrees. 96 degrees is where L. Hoffmann's retrieval switches (per
C. Wright, 2026-11-25), not a property of the Sun.

Origin
------
Port of ``which_airs_retrieval.m`` (corwin365/MatlabFunctions, commit
f1a78c0; C. Wright, 2018), itself a translation of L. Hoffmann's C code, used
with the group's permission. The optional ``RowFlag`` majority rule of the
MATLAB is not ported.

Validation
----------
``tests/test_solar.py``: against MATLAB R2024b for 2166 places and times
(``tests/data/matlab_swath_fixtures.mat``); known answers (Sun overhead near
the Tropic of Cancer at noon at the June solstice; the angle and its
supplement on opposite sides of the Earth). Tolerances stated in the test.

AI involvement
--------------
AI-ported (Claude) from the MATLAB; to be reviewed by C. Wright.

Usage
-----
>>> import numpy as np
>>> from atmosgw.curate.solar import solar_zenith_angle
>>> t = np.datetime64("2020-06-21T12:00")
>>> bool(solar_zenith_angle(t, 23.4, 0.0) < 1)
True

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

_J2000_NOON = np.datetime64("2000-01-01T12:00:00", "us")
_J2000_NOON_DATENUM = 730486.5  # MATLAB datenum(2000,1,1,12,0,0)


def solar_zenith_angle(time, lat, lon, datenum=False):
    """Solar zenith angle in degrees. See the module docstring."""
    if datenum:
        d = np.asarray(time, dtype=float) - _J2000_NOON_DATENUM
    else:
        t = np.asarray(time).astype("datetime64[us]")
        d = (t - _J2000_NOON) / np.timedelta64(1, "D")
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)

    g = (357.529 + 0.98560028 * d) * np.pi / 180
    q = 280.459 + 0.98564736 * d
    ecl = (q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g)) * np.pi / 180
    obl = (23.439 - 0.00000036 * d) * np.pi / 180
    dec = np.arcsin(np.sin(obl) * np.sin(ecl))
    ra = np.arctan2(np.cos(obl) * np.sin(ecl), np.cos(ecl))
    gmst = 18.697374558 + 24.06570982441908 * d
    lst = gmst + lon / 15
    h = lst / 12 * np.pi - ra
    lat_r = lat * np.pi / 180
    cos_sza = np.sin(lat_r) * np.sin(dec) + np.cos(lat_r) * np.cos(dec) * np.cos(h)
    return np.arccos(cos_sza) * 180 / np.pi


def day_retrieval(time, lat, lon, threshold=96.0, datenum=False):
    """1 where the day retrieval is used, 0 for night, NaN where undefined."""
    sza = np.asarray(solar_zenith_angle(time, lat, lon, datenum))
    out = np.full(sza.shape, np.nan)
    out[sza < threshold] = 1.0
    out[sza >= threshold] = 0.0
    return out
