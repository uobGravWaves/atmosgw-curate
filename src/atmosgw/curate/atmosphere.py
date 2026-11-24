"""
Approximate conversion between height and pressure.

Purpose
-------
Convert altitude to pressure and back with a fixed reference atmosphere, for
datasets that give only one of the two. Faithful ports of the group's MATLAB
``h2p`` and ``p2h``, including their known limitations (Q1-Q3 below), so that
ported loaders match the MATLAB reference outputs. To be replaced by a better
method later; callers should not depend on the limitations.

Inputs and outputs
------------------
``height_to_pressure(altitude)``: altitude in km (any shape) -> pressure in
hPa, same shape, float64.

``pressure_to_height(pressure)``: pressure in hPa (any shape) -> altitude in
km, same shape, float64.

Science
-------
Pressures at seven fixed levels, joined by exponentials:

=========  =========
z (km)     p (hPa)
=========  =========
0          1013.25
11         226.3210
20         54.7489
32         8.6802
47         1.1091
51         0.6694
71         0.0396
=========  =========

Between levels i and i+1, p = p_i exp(-(z - z_i) / H_i), with scale height
H_i = -(z_{i+1} - z_i) / ln(p_{i+1} / p_i). Above 71 km the 51-71 km scale
height is used. ``pressure_to_height`` is the inverse, z = -H ln(p / p_i) + z_i.

The level pressures are based on the U.S. Standard Atmosphere, 1976 [R1]
(confirmed by C. Wright; the values themselves are not yet checked against
the source document: UNVERIFIED).

Known limitations, kept for agreement with MATLAB (see ``docs/quirks.md``):

- **Q1** Every layer is treated as isothermal (a single exponential), whereas
  the 1976 standard has a temperature gradient in most layers; and the level
  heights are used as geometric heights, whereas the standard's layer heights
  are probably geopotential (UNVERIFIED). Both make the conversion
  approximate, increasingly so with height.
- **Q2** ``height_to_pressure`` returns NaN below 0 km.
- **Q3** ``pressure_to_height`` handles pressures above 1013.25 hPa only by
  reusing the layer of the previous finite element, taking elements in MATLAB
  (column-major) order; if there is none, it raises ``ValueError`` (MATLAB
  stops with an error).

Negative pressures are not supported (MATLAB returns complex numbers for
them); NaN in gives NaN out.

Precision, as MATLAB: ``pressure_to_height`` of float32 pressures computes in
single precision, as MATLAB's ``p2h`` does for single input, and returns the
single-precision heights as float64. MATLAB's behaviour was found by
comparison of each intermediate step with MATLAB R2024b, not from
documentation: each single-with-double operation is computed in double and
rounded once to single. The logarithm is computed in double and rounded to
single; MATLAB's own single-precision log is not correctly rounded and
differs from this by one unit in the last place for about 4% of inputs, so
about 2% of heights can differ by one float32 unit (about 1e-7 relative). For
the 55 MLS pressure levels all heights match MATLAB exactly (NumPy's float32
arithmetic alone differs at 7 of them). Used for MLS, whose pressures are
float32. Layer choice compares the float32 values exactly with the double
level pressures; MATLAB's choice for a value exactly at a level is
UNVERIFIED.
``height_to_pressure`` always computes in float64 (its inputs in the loaders
are double).

Origin
------
Ports of ``h2p.m`` and ``p2h.m`` (corwin365/MatlabFunctions, the group's own
code, commit c943399).

Validation
----------
``tests/test_atmosphere.py``: against MATLAB R2024b output at the level
heights and pressures, either side of them, random values, NaN, matrices, and
the Q2 and Q3 cases; known answers (level pressures exactly, round trip).
Tolerance stated in the test.

AI involvement
--------------
AI-ported (Claude) from the MATLAB, to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.atmosphere import height_to_pressure
>>> float(height_to_pressure(11.0))
226.321

Status
------
Preliminary; to be replaced (see Q1).

Change log
----------
2026/11/24 [Claude] initial version (port of h2p.m and p2h.m)
2026/11/24 [Claude] pressure_to_height computes in float32 for float32 input, as MATLAB

How to cite
-----------
Cite the package (see ``CITATION.cff``) and, for the level values, [R1].
"""

import numpy as np

# level heights (km) and pressures (hPa), as in h2p.m and p2h.m
LEVEL_HEIGHTS = np.array([0, 11, 20, 32, 47, 51, 71], dtype=float)
LEVEL_PRESSURES = np.array(
    [1013.25, 226.3210, 54.7489, 8.6802, 1.1091, 0.6694, 0.0396], dtype=float
)
# scale height of each layer (km); the last also used above 71 km
SCALE_HEIGHTS = -np.diff(LEVEL_HEIGHTS) / np.log(
    LEVEL_PRESSURES[1:] / LEVEL_PRESSURES[:-1]
)


def height_to_pressure(altitude):
    """Altitude (km) to pressure (hPa). See the module docstring (port of h2p)."""
    alt = np.asarray(altitude, dtype=float)
    z, p, h = LEVEL_HEIGHTS, LEVEL_PRESSURES, SCALE_HEIGHTS
    out = np.full(alt.shape, np.nan)
    for i in range(z.size - 1):
        idx = (alt >= z[i]) & (alt < z[i + 1])
        out[idx] = p[i] * np.exp(-(alt[idx] - z[i]) / h[i])
    top = alt >= z[-1]
    out[top] = p[-1] * np.exp(-(alt[top] - z[-1]) / h[-1])
    return out[()] if out.ndim == 0 else out


def pressure_to_height(pressure):
    """Pressure (hPa) to altitude (km). See the module docstring (port of p2h)."""
    pres = np.asarray(pressure)
    # MATLAB evaluates single op double in single precision, so p2h of single
    # pressures computes heights in single precision (returned as double)
    work = np.float32 if pres.dtype == np.float32 else np.float64
    pres = pres.astype(work)
    shape = pres.shape
    pres = pres.ravel(order="F")  # p2h.m works through pres(:), column-major
    z, p, h = LEVEL_HEIGHTS, LEVEL_PRESSURES, SCALE_HEIGHTS
    if np.any(pres < 0):
        raise ValueError("Negative pressures are not supported.")

    # layer of each element, as in p2h.m: layer k holds p[k] >= pres > p[k+1];
    # the top layer (index 6, using the 51-71 km scale height) holds pres <= p[6]
    layer = np.full(pres.shape, -1)
    for k in range(z.size - 1):
        layer[(pres <= p[k]) & (pres > p[k + 1])] = k
    layer[pres <= p[-1]] = z.size - 1

    # Q3: above 1013.25 hPa no layer matches, and p2h.m reuses the previous
    # finite element's layer (NaN elements are skipped and leave it unchanged)
    finite = np.isfinite(pres)
    unmatched = finite & (layer < 0)
    if np.any(unmatched):
        index = np.where(finite & (layer >= 0), np.arange(pres.size), -1)
        last_valid = np.maximum.accumulate(index)
        if np.any(last_valid[unmatched] < 0):
            raise ValueError(
                "Pressure above 1013.25 hPa with no earlier element to take a layer "
                "from (the MATLAB p2h stops with an error here)."
            )
        layer[unmatched] = layer[last_valid[unmatched]]

    ref_pres = np.append(p, p[-1])[layer.clip(0)]
    ref_alt = np.append(z, z[-1])[layer.clip(0)]
    scale = np.append(h, h[-1])[layer.clip(0)]
    out = np.full(pres.shape, np.nan)
    f = finite
    if work == np.float64:
        out[f] = (-scale[f] * np.log(pres[f] / ref_pres[f])) + ref_alt[f]
    else:
        # as MATLAB with a single input: each single-with-double operation is
        # computed in double and rounded once to single; log is computed in
        # double and rounded (MATLAB's own single log can differ by 1 unit; see
        # the module docstring). Matches MATLAB for all 55 MLS levels.
        def single(x):
            return np.asarray(x, dtype=np.float64).astype(np.float32)

        ratio = single(pres[f].astype(np.float64) / ref_pres[f])
        log = single(np.log(ratio.astype(np.float64)))
        product = single(-scale[f] * log.astype(np.float64))
        out[f] = single(product.astype(np.float64) + ref_alt[f])
    out = out.reshape(shape, order="F")
    return out[()] if out.ndim == 0 else out
