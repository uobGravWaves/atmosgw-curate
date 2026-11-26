"""
Statistical outlier screens for nadir-sounder swaths, swappable.

A screen marks values that are statistical outliers for their height level.
Screens are interchangeable components: each is a function with the same
signature, registered by name in ``SCREENS``, so a better method can replace
one later without changing the loader:

    screen(temp, time, lat, lon, **settings) -> bool array

``temp``: [cross x along x level]; ``time``: [cross x along] datenums;
``lat``, ``lon``: [cross x along] degrees. Returns True at outliers, shape of
``temp``. A screen only marks points; the loader records them as the
``statistical_outlier`` state and sets them to NaN (unless
``keep_outliers``).

Choosing a screen: each instrument's registry entry has a default
(``screen``: a name and settings); ``load_nadir_sounders(screen=...)``
overrides it per call with a name, a ``(name, settings)`` pair, a function
with the signature above, or ``"none"``. ``resolve`` turns any of these into
``(function, settings, description, pad_minutes)``; the description is
recorded in the output attribute ``screen``.

Registered screens
------------------
``"rolling"`` (default for IASI; agreed with C. Wright, 2026-11-25):

- scans are grouped into blocks of ``step_minutes`` (1 min) on a fixed UTC
  grid, by the mean time of each along-track row;
- for each block, at each height level, the quartiles are taken over every
  point in the rows within ``window_minutes`` (10 min) centred on the block,
  separately for day and night points if ``split_day_night`` (solar zenith
  angle below or above 96 degrees, where the retrieval switches between its
  day and night set-ups: L. Hoffmann's AIRS retrieval; believed the same for
  IASI by C. Wright, to be confirmed with L. Hoffmann: UNVERIFIED);
- a block's points are outliers if more than ``factor`` (5) interquartile
  ranges beyond the quartiles of their own group (5, agreed 2026-11-25: at 3,
  about 2 % of IASI flags were within 50 K of the level median, most of them
  wave-like; at 5, 3 of 12,473);
- a level is not screened for a group in a window with fewer than
  ``min_points`` (500) finite values of that group, or where the
  interquartile range is zero (a safeguard for quantised data, where both
  quartiles could be equal and every other value would be flagged).

The loader reads extra data beyond each end of a request (``pad_minutes``:
half a window plus half a granule) for this screen, then trims it, so a
scan's flags do not depend on how the data were requested.

``"pooled"``: the MATLAB ``prep_airs_3d`` screen (fixed 2026-11-25): the
quartiles of each height level over every point loaded in one call, ``factor``
(3) interquartile ranges. Kept for comparison with MATLAB; its flags depend
on what is loaded together.

``"none"``: no screen.

Defaults, from measurements on 2026-11-25 (3 h on each of 2007-01-01,
2010-07-01, 2013-01-15 and 2016-06-01, heights 24-57 km). IASI at 3
interquartile ranges: a median 0.05 % (max 1.6 %) of points per granule
flagged, 98 % of them more than 50 K from their level's median (median
127 K), but 301 within 50 K, 181 of those wave-like (at least 3 of 4
neighbours deviating the same way); at 5: 12,473 flagged (97.5 % of those at
3), only 3 within 50 K. Hence 5. AIRS at 3: 0-0.16 % flagged, none more than
50 K off (median 8 K): real variability in AIRS's narrow distributions (the
IASI IQR at these heights is typically 15 K), so AIRS has no screen
(``"none"``).

Quartiles in both: as MATLAB ``prctile`` (NumPy's "hazen": linear
interpolation at (k - 0.5)/n), NaN ignored.

Origin
------
``pooled``: port of the IASI screen of ``prep_airs_3d.m`` (corwin365/
MatlabFunctions, commit f1a78c0). ``rolling``: new.

Status
------
AI-originated (Claude), from the method agreed with C. Wright; to be
reviewed by C. Wright.

Change log
----------
2026/11/25 [Claude] initial version: rolling, pooled, none; swappable
"""

import warnings

import numpy as np

from atmosgw.curate.solar import day_retrieval


def rolling_screen(temp, time, lat, lon, window_minutes=10.0, step_minutes=1.0,
                   factor=5.0, min_points=500, split_day_night=True):
    """Rolling-window screen, day and night separately. See the module docstring."""
    temp = np.asarray(temp, dtype=float)
    time = np.asarray(time, dtype=float)
    out = np.zeros(temp.shape, dtype=bool)
    if temp.size == 0:
        return out
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        row_time = np.nanmean(time, axis=0)
    step = step_minutes / 1440
    half = window_minutes / 2 / 1440
    block = np.floor(row_time / step)
    if split_day_night:
        day = day_retrieval(time, lat, lon, datenum=True)
        groups = [day == 1, day == 0]
    else:
        groups = [np.ones(time.shape, dtype=bool)]

    order = np.argsort(row_time, kind="stable")
    sorted_time = row_time[order]
    for b in np.unique(block[np.isfinite(block)]):
        centre = (b + 0.5) * step
        lo = np.searchsorted(sorted_time, centre - half, side="left")
        hi = np.searchsorted(sorted_time, centre + half, side="left")
        rows_w = order[lo:hi]
        rows_b = np.flatnonzero(block == b)
        window = temp[:, rows_w, :]
        target = temp[:, rows_b, :]
        for group in groups:
            in_w = group[:, rows_w]
            in_b = group[:, rows_b]
            if not in_b.any():
                continue
            values = window[in_w]  # [points x level]
            counts = np.sum(np.isfinite(values), axis=0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                q1, q3 = np.nanpercentile(values, [25, 75], axis=0, method="hazen")
            iqr = q3 - q1
            usable = (counts >= min_points) & (iqr > 0)  # IQR 0: safeguard
            xi, ai = np.nonzero(in_b)
            points = target[xi, ai, :]  # [points x level]
            with np.errstate(invalid="ignore"):
                bad = ((points < q1 - factor * iqr) | (points > q3 + factor * iqr))
            bad &= usable[None, :]
            out[xi, rows_b[ai], :] = bad
    return out


def pooled_screen(temp, time, lat, lon, factor=3.0):
    """MATLAB's screen: quartiles of each level over all points. See above."""
    temp = np.asarray(temp, dtype=float)
    flat = temp.reshape(-1, temp.shape[-1])
    if flat.size == 0:
        return np.zeros(temp.shape, dtype=bool)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q1, q3 = np.nanpercentile(flat, [25, 75], axis=0, method="hazen")
    iqr = q3 - q1
    with np.errstate(invalid="ignore"):
        out = (flat < q1 - factor * iqr) | (flat > q3 + factor * iqr)
    return out.reshape(temp.shape)


SCREENS = {"rolling": rolling_screen, "pooled": pooled_screen, "none": None}

# extra data read beyond each end of a request, per screen (minutes)
PADDING = {"rolling": lambda s: s.get("window_minutes", 10.0) / 2 + 3.0}


def resolve(spec, default):
    """
    (function, settings, description, pad_minutes) for a screen choice.

    spec: None (use default), a name, (name, settings), or a function.
    default: the instrument's (name, settings).
    """
    if spec is None:
        spec = default
    if callable(spec):
        name = getattr(spec, "__name__", "custom")
        return spec, {}, f"custom: {name}", 0.0
    if isinstance(spec, str):
        spec = (spec, {})
    name, settings = spec
    if name not in SCREENS:
        raise ValueError(f"Unknown screen {name!r}; choose from {', '.join(SCREENS)}, "
                         "or pass a function.")
    function = SCREENS[name]
    if function is None:
        return None, {}, "none", 0.0
    text = ", ".join(f"{k}={v}" for k, v in sorted(settings.items()))
    pad = PADDING.get(name, lambda s: 0.0)(settings)
    return function, dict(settings), f"{name}({text})", pad
