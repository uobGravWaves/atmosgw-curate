"""
Box and Gaussian smoothing of N-D arrays, as the MATLAB pre-smoothing step.

``smooth_array(values, size, kernel="box", sigma=0.65, skip_nan=False)``
smooths a NumPy array with a separable kernel of odd ``size`` along each
axis. It is the smoothing step of ``analysis.detrend`` (step ``"smooth"``),
where the options are described.

Inputs and outputs
------------------
- ``values``: array of any shape; ``size``: one odd integer >= 1 per axis
  (1 = no smoothing along that axis).
- ``kernel``: ``"box"`` (equal weights) or ``"gaussian"``.
- ``sigma`` (float or one per axis, 0.65): Gaussian only: the standard
  deviation as a fraction of the kernel's half-width, so in grid points it
  is ``sigma * (size - 1) / 2``. The kernel is cut off at ``size`` points
  and normalised to sum to 1.
- Edges: the array is extended by repeating its edge values.
- NaN: a NaN makes every output within the kernel's reach NaN (as MATLAB).
  With ``skip_nan=True`` each output is the weighted mean of the values
  that are present, and is NaN only where the input is NaN.
- Returns a float64 array of the same shape.

Origin
------
New code, written from the documented behaviour of the ``smoothn`` used by
the group's ``prep_airs_3d`` (A. Gannepali, MATLAB File Exchange 725: box or
Gaussian kernel, edges padded by replication); no source was translated
(CONTRIBUTING, "Code from other sources": third-party code).

Validation
----------
``tests/test_detrend_matlab.py``: against MATLAB ``smoothn`` outputs for box
and Gaussian kernels (tolerance stated there); ``tests/test_detrend.py``:
known answers and ``skip_nan``.

Status
------
Matching MATLAB (see Validation). AI-written (Claude); to be reviewed by
C. Wright.

Change log
----------
2026/11/28 [Claude] initial version

"""

import numpy as np
from scipy.ndimage import correlate1d


def smooth_array(values, size, kernel="box", sigma=0.65, skip_nan=False):
    """Smooth values with a box or Gaussian kernel. See the module docstring."""
    values = np.asarray(values, dtype=np.float64)
    size = _per_axis(size, values.ndim, "size")
    if any(int(s) != s or s < 1 or int(s) % 2 == 0 for s in size):
        raise ValueError(f"size must be odd integers >= 1, one per axis; got {size}")
    size = [int(s) for s in size]
    kind = str(kernel).lower()
    if kind not in ("box", "gaussian"):
        raise ValueError(f"kernel must be 'box' or 'gaussian', not {kernel!r}")
    sigma = _per_axis(sigma, values.ndim, "sigma")
    weights = [_weights(n, kind, s) for n, s in zip(size, sigma)]
    if not skip_nan:
        return _apply(values, weights)
    present = np.isfinite(values)
    total = _apply(np.where(present, values, 0.0), weights)
    share = _apply(present.astype(np.float64), weights)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = total / share
    out[~present] = np.nan
    return out


# supporting functions


def _per_axis(value, ndim, name):
    """value as a list with one entry per axis (a scalar is repeated)."""
    if np.ndim(value) == 0:
        return [value] * ndim
    value = list(value)
    if len(value) != ndim:
        raise ValueError(f"{name} needs one value per axis ({ndim}); got "
                         f"{len(value)}")
    return value


def _weights(n, kind, sigma):
    """The 1-D kernel of n points (summing to 1)."""
    half = (n - 1) // 2
    if half == 0:
        return np.ones(1)
    if kind == "box":
        return np.full(n, 1.0 / n)
    if not sigma > 0:
        raise ValueError(f"sigma must be positive; got {sigma}")
    offsets = np.arange(-half, half + 1, dtype=np.float64)
    weights = np.exp(-(offsets * offsets) / (2.0 * (sigma * half) ** 2))
    return weights / weights.sum()


def _apply(values, weights):
    """Apply the 1-D kernels along each axis in turn, edges repeated."""
    out = values
    for axis, w in enumerate(weights):
        if w.size > 1:
            out = correlate1d(out, w, axis=axis, mode="nearest")
    return out
