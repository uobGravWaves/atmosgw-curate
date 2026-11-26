"""
MATLAB-compatible versions of a few numerical functions.

Reproduce the behaviour of MATLAB's ``interp1``, ``unwrap`` and ``wrapToPi``
where NumPy's nearest equivalents behave differently, so that ports of the
group's MATLAB code give the same results. Use these in ports; new code may
prefer the NumPy/SciPy functions directly.

Inputs and outputs
------------------
``interp1(x, v, xq, method="linear", extrapolate=None)``

===============  ===========  ========  ==========================================
Option           Type         Default   Description
===============  ===========  ========  ==========================================
``x``            1-D array    -         sample points, distinct, any order
``v``            1-D array    -         sample values, same length as ``x``
``xq``           array        -         query points, any shape
``method``       str          "linear"  "linear", "nearest", "next", "previous",
                                        "pchip", "spline" or "makima"
``extrapolate``  bool/None    None      None: MATLAB's default for the method;
                                        True: like ``interp1(..., 'extrap')``;
                                        False: NaN outside the sample range
===============  ===========  ========  ==========================================

Returns a float64 array with the shape of ``xq``.

``unwrap(p)``: 1-D array of finite angles (radians) -> unwrapped copy.

``wrap_to_pi(angle)``: array of angles (radians) -> wrapped to [-pi, pi].

MATLAB behaviour reproduced (each checked against MATLAB R2024b output in
``tests/data/matlab_compat_fixtures.mat``):

- ``interp1``:

  - sample points may be in any order (MATLAB sorts them; ``np.interp``
    silently gives wrong answers for unsorted or descending points);
  - by default, NaN outside the sample range for linear, nearest, next and
    previous, but extrapolation for pchip, spline and makima (``np.interp``
    instead repeats the end values);
  - NaN sample values: for linear, nearest, next and previous they spread
    into the neighbouring intervals, but a query exactly at a sample point
    returns that sample's value; for pchip, spline and makima, NaN samples
    are dropped first (MATLAB warns "Columns of data containing NaN values
    have been ignored");
  - nearest: a query exactly midway between two samples takes the upper one;
  - extrapolation with nearest uses the end values; next extrapolates only
    below the range and previous only above it;
  - linear interpolation is computed as ``(1 - t) * v0 + t * v1``, except
    that between two equal values the result is exactly that value; this
    matches MATLAB bit for bit (NumPy's formula differs in the last bit).
    Checked on 80,000 further MATLAB results, including constant, step-wise
    and nearly constant data; the equal-value rule is checked inside the
    sample range only (with extrapolation it is assumed: UNVERIFIED).

- ``unwrap``: jumps of pi or more between neighbours are removed by adding
  whole multiples of 2*pi (a jump of exactly +-pi is left as it is). The
  correction is accumulated as a whole number of cycles and then multiplied by
  2*pi, which matches MATLAB bit for bit (``np.unwrap`` differs in the last
  bit).

- ``wrap_to_pi``: values already in [-pi, pi] are returned unchanged; positive
  odd multiples of pi become pi and negative odd multiples become -pi.

Not supported: ``interp1`` methods "cubic" and "v5cubic" (in R2024b "cubic"
is cubic convolution, which needs evenly spaced points and falls back to
"spline" otherwise); repeated or NaN sample points; matrices of values;
``unwrap`` with NaN or a tolerance other than pi.

``deg2rad`` and ``rad2deg`` need no copies: ``np.deg2rad`` and ``np.rad2deg``
match MATLAB bit for bit.

Origin
------
New code written to match MATLAB's documented behaviour and checked against
MATLAB output; not a translation of MathWorks source code. Pchip, spline and
makima use ``scipy.interpolate`` (``PchipInterpolator``, ``CubicSpline`` with
not-a-knot end conditions, ``Akima1DInterpolator(method="makima")``).

Validation
----------
``tests/test_matlab_compat.py``: against MATLAB R2024b Update 6 for ascending,
descending and unsorted samples, NaN values and queries, midpoint ties, three
samples and a limb-sounder-like profile, every method with and without
extrapolation; unwrap with random jumps, exact +-pi, 2*pi and 3*pi jumps and
dateline-crossing longitudes; wrap_to_pi at multiples of pi. Linear, nearest,
next, previous, unwrap and wrap_to_pi are bit-identical; the tolerance for
pchip, spline and makima is stated in the test.

Status
------
Preliminary. AI-originated (Claude), to be reviewed by C. Wright.

Usage
-----
>>> import numpy as np
>>> from atmosgw.curate.matlab_compat import interp1
>>> interp1(np.array([3.0, 2.0, 1.0]), np.array([30.0, 20.0, 10.0]),
...         np.array([0.5, 1.5, 2.5]))
array([nan, 15., 25.])

Change log
----------
2026/11/24 [Claude] initial version

"""

import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator

_PIECEWISE = ("linear", "nearest", "next", "previous")
_SMOOTH = ("pchip", "spline", "makima")


def interp1(x, v, xq, method="linear", extrapolate=None):
    """MATLAB-compatible 1-D interpolation. See the module docstring."""
    method = method.lower()
    if method in ("cubic", "v5cubic"):
        raise NotImplementedError(
            f"interp1 method {method!r} is not supported (MATLAB's 'cubic' is cubic "
            "convolution on evenly spaced points); use 'pchip' or 'spline'."
        )
    if method not in _PIECEWISE + _SMOOTH:
        raise ValueError(f"Unknown interp1 method {method!r}.")

    x = np.asarray(x, dtype=float)
    v = np.asarray(v, dtype=float)
    xq = np.asarray(xq, dtype=float)
    if x.ndim != 1 or v.shape != x.shape:
        raise ValueError("x and v must be 1-D arrays of the same length.")
    if np.any(np.isnan(x)):
        raise ValueError("Sample points x must not contain NaN.")

    # pchip, spline and makima ignore NaN samples; the others keep them
    if method in _SMOOTH:
        keep = ~np.isnan(v)
        x, v = x[keep], v[keep]

    # sort the sample points, as MATLAB does
    order = np.argsort(x, kind="stable")
    x, v = x[order], v[order]
    if np.any(np.diff(x) == 0):
        raise ValueError("Sample points x must be distinct.")
    if x.size < 2:
        raise ValueError("At least two sample points are needed.")

    if extrapolate is None:
        extrapolate = method in _SMOOTH

    if method in _SMOOTH:
        if method == "pchip":
            f = PchipInterpolator(x, v, extrapolate=True)
        elif method == "spline":
            f = CubicSpline(x, v, bc_type="not-a-knot", extrapolate=True)
        else:
            f = Akima1DInterpolator(x, v, method="makima", extrapolate=True)
        out = f(xq)
        if not extrapolate:
            out[(xq < x[0]) | (xq > x[-1])] = np.nan
        out[np.isnan(xq)] = np.nan
        return out

    return _interp_piecewise(x, v, xq, method, extrapolate)


def unwrap(p):
    """MATLAB-compatible unwrap of a 1-D array of angles. See the module docstring."""
    p = np.asarray(p, dtype=float)
    if p.ndim != 1:
        raise ValueError("unwrap takes a 1-D array.")
    if not np.all(np.isfinite(p)):
        raise ValueError("unwrap does not support NaN or infinite values.")
    if p.size < 2:
        return p.copy()

    jump = np.diff(p)
    # each jump brought into [-pi, pi); a jump of +pi stays +pi
    wrapped = np.mod(jump + np.pi, 2 * np.pi) - np.pi
    wrapped[(wrapped == -np.pi) & (jump > 0)] = np.pi
    # whole number of cycles to add; none for jumps smaller than pi
    cycles = np.round((wrapped - jump) / (2 * np.pi))
    cycles[np.abs(jump) < np.pi] = 0

    out = p.copy()
    out[1:] = p[1:] + 2 * np.pi * np.cumsum(cycles)
    return out


def unwrap_rows(p):
    """unwrap of the finite values of each row of a 2-D array, all rows at
    once (NaN stay NaN). The same arithmetic as unwrap() on each row's finite
    values, so the results are identical."""
    p = np.asarray(p, dtype=float)
    out = p.copy()
    if p.ndim != 2 or p.shape[1] < 2:
        return out
    finite = np.isfinite(p)
    # each value minus the previous finite value of its row (the jumps unwrap
    # sees between successive finite values); 0 elsewhere
    index = np.where(finite, np.arange(p.shape[1])[None, :], 0)
    np.maximum.accumulate(index, axis=1, out=index)
    previous = np.take_along_axis(p, index, axis=1)
    jump = np.zeros(p.shape)
    with np.errstate(invalid="ignore"):
        jump[:, 1:] = p[:, 1:] - previous[:, :-1]
    jump[:, 1:][~finite[:, 1:] | ~np.isfinite(previous[:, :-1])] = 0.0
    with np.errstate(invalid="ignore"):
        wrapped = np.mod(jump + np.pi, 2 * np.pi) - np.pi
    wrapped[(wrapped == -np.pi) & (jump > 0)] = np.pi
    cycles = np.round((wrapped - jump) / (2 * np.pi))
    cycles[np.abs(jump) < np.pi] = 0
    total = np.cumsum(cycles, axis=1)
    # (as unwrap: every finite value but the first of its row gets the shift,
    # even a zero one, which turns -0.0 into 0.0)
    shift = finite.copy()
    shift[np.arange(p.shape[0]), np.argmax(finite, axis=1)] = False
    out[shift] = p[shift] + 2 * np.pi * total[shift]
    return out


def wrap_to_pi(angle):
    """MATLAB-compatible wrapToPi (radians). See the module docstring."""
    angle = np.array(angle, dtype=float)  # a copy
    outside = (angle < -np.pi) | (angle > np.pi)
    shifted = angle[outside] + np.pi
    wrapped = np.mod(shifted, 2 * np.pi)
    # as wrapTo2Pi: positive multiples of 2*pi map to 2*pi, not 0
    wrapped[(wrapped == 0) & (shifted > 0)] = 2 * np.pi
    angle[outside] = wrapped - np.pi
    return angle[()] if angle.ndim == 0 else angle


# supporting functions


def _interp_piecewise(x, v, xq, method, extrapolate):
    """Linear, nearest, next and previous on sorted, distinct x."""
    out = np.full(xq.shape, np.nan)
    finite = ~np.isnan(xq)
    q = xq[finite]
    below, above = q < x[0], q > x[-1]

    # interval j holds x[j] <= q < x[j+1]; the last sample goes in the last one
    j = np.clip(np.searchsorted(x, q, side="right") - 1, 0, x.size - 2)
    x0, x1, v0, v1 = x[j], x[j + 1], v[j], v[j + 1]

    if method == "linear":
        t = (q - x0) / (x1 - x0)
        r = (1 - t) * v0 + t * v1
        # between equal values, exactly that value (MATLAB does this; the
        # formula above can be one bit off)
        r = np.where(v0 == v1, v0, r)
        # at a sample point, its own value (even if a neighbour is NaN)
        r = np.where(q == x0, v0, np.where(q == x1, v1, r))
    elif method == "nearest":
        r = np.where(q - x0 < x1 - q, v0, v1)  # a tie goes to the upper sample
        r = np.where(below, v[0], np.where(above, v[-1], r))
    elif method == "next":
        r = np.where(q == x0, v0, v1)
        r = np.where(below, v[0], np.where(above, np.nan, r))
    else:  # previous
        r = np.where(q == x1, v1, v0)
        r = np.where(below, np.nan, np.where(above, v[-1], r))

    if not extrapolate:
        r = np.where(below | above, np.nan, r)
    out[finite] = r
    return out
