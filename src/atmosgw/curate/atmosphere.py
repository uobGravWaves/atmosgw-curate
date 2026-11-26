"""
Conversion between height and pressure with the U.S. Standard Atmosphere, 1976.

Convert altitude to pressure and back, for datasets that give only one of the
two. Ports of the group's MATLAB ``h2p`` and ``p2h``, which since 2026-11-25
(both languages) implement the U.S. Standard Atmosphere, 1976 [R1] instead of
the earlier approximation (Q1-Q3 in ``docs/quirks.md``).

Inputs and outputs
------------------
``height_to_pressure(altitude)``: geometric altitude in km (any shape) ->
pressure in hPa, same shape, float64.

``pressure_to_height(pressure)``: pressure in hPa (any shape) -> geometric
altitude in km, same shape, float64.

Valid from -5 km to 149 km geometric (1777.6 to 4.74e-6 hPa), the range of
Table I of [R1]; NaN outside that range, and for NaN or negative
input. Always computed in float64.

Science
-------
All references are to [R1], the U.S. Standard Atmosphere, 1976
(NOAA-S/T 76-1562 / NASA-TM-X-74335; ``docs/theory/relations.md``); pages are
its printed page numbers.

Up to 86 km geometric (84.852 km' geopotential), the standard's defining
equations:

- geometric altitude Z and geopotential altitude H (km'), eqs (18) and (19)
  (p. 8): H = r0 Z / (r0 + Z), Z = r0 H / (r0 - H), with r0 = 6356.766 km
  (table 2, p. 2);
- molecular-scale temperature linear in H in seven layers, eq (23) (p. 10):
  TM = TM,b + LM,b (H - Hb), with base heights Hb and gradients LM,b from
  table 4 (p. 3) and TM,0 = 288.15 K (table 2);
- pressure from the hydrostatic equation, eqs (33a) and (33b) (p. 12):
  P = Pb [TM,b / (TM,b + LM,b (H - Hb))]^(g0' M0 / (R* LM,b)) where
  LM,b is not 0, and P = Pb exp(-g0' M0 (H - Hb) / (R* TM,b)) where it is,
  with g0' = 9.80665 m^2 s^-2 m'^-1, R* = 8.31432e3 N m kmol^-1 K^-1,
  P0 = 1013.25 hPa (table 2) and M0 = 28.9644 kg kmol^-1 (p. 9). The base
  pressures Pb of the higher layers follow from the same equations at H = Hb
  (p. 12). Below sea level eq (33a) is used with the first layer, as for
  Table I, which starts at -5 km geometric (-5.004 km'; p. 51).

Above 86 km the standard's pressure comes from a diffusion model (eq (33c),
p. 12), not a closed formula, so pressures are interpolated in log(pressure)
between the tabulated values of Table I (geometric altitude, p. 68): every
0.5 km from 86 to 100 km and every 1 km from 100 to 149 km. Within a step
log(pressure) is very nearly linear in height (its slope changes by at most
0.0066 per km between steps).

``pressure_to_height`` is the exact inverse: the inverted eqs (33a), (33b)
and (19) below 86 km, and the same table interpolation (Z linear in
log(pressure)) above.

Checks against the document (2026-11-25): the equations reproduce the printed
Table I pressures at -5, 40, 50, 55, 60, 70, 80, 85, 85.5 and 86 km to within
their 5 printed digits (largest difference 7e-5 relative), and the table
values used above 86 km were each checked against the table's independent
P/P0 column (5 digits). At 86 km the equations and the table differ by
1.2e-6 relative (about 0.2 m).

The standard represents mean mid-latitude conditions; the real atmosphere
departs from it with latitude, season and weather.

Origin
------
Ports of ``h2p.m`` and ``p2h.m`` (corwin365/MatlabFunctions), rewritten
2026-11-25 in both languages to implement the standard; the group's own code.

Validation
----------
``tests/test_atmosphere.py``: against MATLAB R2024b output of the rewritten
``h2p.m`` and ``p2h.m``, and known answers from Table I of [R1]; round trip.
Tolerances stated in the test.

Status
------
Implements the standard. Written by Claude from the standard's document
(constants and equations checked against it, as above); to be reviewed by
C. Wright.

Usage
-----
>>> from atmosgw.curate.atmosphere import height_to_pressure
>>> round(float(height_to_pressure(11.0)), 2)
227.0

Change log
----------
2026/11/24 [Claude] initial version (port of h2p.m and p2h.m)
2026/11/24 [Claude] pressure_to_height computes in float32 for float32 input, as MATLAB
2026/11/25 [Claude] U.S. Standard Atmosphere 1976 to 149 km (Q1-Q3), as new h2p.m/p2h.m

"""

import numpy as np

# constants of the standard [R1]: table 2 (p. 2) and p. 9
G0 = 9.80665  # g0', m^2 s^-2 m'^-1
GAS_CONSTANT = 8.31432e3  # R*, N m kmol^-1 K^-1
MOLECULAR_WEIGHT = 28.9644  # M0, kg kmol^-1
EARTH_RADIUS = 6356.766  # r0, km
SURFACE_TEMPERATURE = 288.15  # T0 = TM,0, K
SURFACE_PRESSURE = 1013.25  # P0, hPa

# layers below 86 km [R1] table 4 (p. 3): base geopotential heights (km') and
# molecular-scale temperature gradients (K/km')
BASE_HEIGHTS = np.array([0, 11, 20, 32, 47, 51, 71, 84.8520])
GRADIENTS = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0])
MIN_ALTITUDE = -5.0  # km geometric, the start of Table I (p. 51)

# [R1] Table I, geometric altitude (p. 68): Z (km) and P (hPa = mb), 86-149 km
TABLE_HEIGHTS = np.concatenate([np.arange(86, 100, 0.5), np.arange(100, 150, 1.0)])
TABLE_PRESSURES = np.array([
    3.7338e-3, 3.4163e-3, 3.1259e-3, 2.8602e-3, 2.6173e-3, 2.3951e-3,
    2.1919e-3, 2.0060e-3, 1.8359e-3, 1.6804e-3, 1.5381e-3, 1.4078e-3,
    1.2887e-3, 1.1798e-3, 1.0801e-3, 9.8896e-4, 9.0560e-4, 8.2937e-4,
    7.5966e-4, 6.9592e-4, 6.3765e-4, 5.8439e-4, 5.3571e-4, 4.9122e-4,
    4.5057e-4, 4.1342e-4, 3.7948e-4, 3.4846e-4,  # 86-99.5 km, 0.5 km steps
    3.2011e-4, 2.7192e-4, 2.3144e-4, 1.9742e-4, 1.6882e-4, 1.4477e-4,
    1.2454e-4, 1.0751e-4, 9.3188e-5, 8.1142e-5, 7.1042e-5, 6.2614e-5,
    5.5547e-5, 4.9570e-5, 4.4473e-5, 4.0096e-5, 3.6312e-5, 3.3022e-5,
    3.0144e-5, 2.7615e-5, 2.5382e-5, 2.3401e-5, 2.1635e-5, 2.0055e-5,
    1.8635e-5, 1.7354e-5, 1.6194e-5, 1.5141e-5, 1.4183e-5, 1.3307e-5,
    1.2505e-5, 1.1769e-5, 1.1092e-5, 1.0468e-5, 9.8907e-6, 9.3568e-6,
    8.8617e-6, 8.4018e-6, 7.9739e-6, 7.5751e-6, 7.2028e-6, 6.8550e-6,
    6.5294e-6, 6.2243e-6, 5.9380e-6, 5.6691e-6, 5.4162e-6, 5.1781e-6,
    4.9538e-6, 4.7421e-6,  # 100-149 km, 1 km steps
])

# g0' M0 / R*, in K per km' (the 1e3 converts km' to m')
_HYDROSTATIC = G0 * MOLECULAR_WEIGHT / GAS_CONSTANT * 1e3


def _layer_bases():
    """Base temperatures (K) and pressures (hPa) of the 7 layers, eqs (23), (33)."""
    temps, pres = [SURFACE_TEMPERATURE], [SURFACE_PRESSURE]
    for b in range(GRADIENTS.size - 1):
        dh = BASE_HEIGHTS[b + 1] - BASE_HEIGHTS[b]
        temps.append(temps[b] + GRADIENTS[b] * dh)
        pres.append(_layer_pressure(pres[b], temps[b], GRADIENTS[b], dh))
    return np.array(temps), np.array(pres)


def _layer_pressure(p_base, t_base, gradient, dh):
    """Pressure dh (km') above a layer base: eq (33a), or (33b) if gradient is 0."""
    if np.ndim(gradient) == 0 and gradient == 0:
        return p_base * np.exp(-_HYDROSTATIC * dh / t_base)
    with np.errstate(divide="ignore", invalid="ignore"):
        grad = np.where(gradient == 0, np.nan, gradient)
        power = p_base * (t_base / (t_base + grad * dh)) ** (_HYDROSTATIC / grad)
    return np.where(gradient == 0, p_base * np.exp(-_HYDROSTATIC * dh / t_base), power)


BASE_TEMPERATURES, BASE_PRESSURES = _layer_bases()
TOP_GEOPOTENTIAL = BASE_HEIGHTS[-1]  # 84.852 km' = 86 km geometric


def height_to_pressure(altitude):
    """Geometric altitude (km) to pressure (hPa). See the module docstring."""
    z = np.asarray(altitude, dtype=float)
    out = np.full(z.shape, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        h = EARTH_RADIUS * z / (EARTH_RADIUS + z)  # eq (18)
    low = np.isfinite(h) & (z >= MIN_ALTITUDE) & (h <= TOP_GEOPOTENTIAL)
    layer = np.clip(np.searchsorted(BASE_HEIGHTS, h[low], side="right") - 1,
                    0, GRADIENTS.size - 1)
    out[low] = _layer_pressure(BASE_PRESSURES[layer], BASE_TEMPERATURES[layer],
                               GRADIENTS[layer], h[low] - BASE_HEIGHTS[layer])
    # the equations end at 84.852 km' = 85.99995 km geometric, the table starts
    # at 86 km: in the 5 cm between, the table's 86 km value is used (np.interp
    # clamps at the ends, deliberately)
    high = np.isfinite(z) & ~low & (h > TOP_GEOPOTENTIAL) & (z <= TABLE_HEIGHTS[-1])
    out[high] = np.exp(np.interp(z[high], TABLE_HEIGHTS, np.log(TABLE_PRESSURES)))
    return out[()] if out.ndim == 0 else out


def pressure_to_height(pressure):
    """Pressure (hPa) to geometric altitude (km). See the module docstring."""
    p = np.asarray(pressure, dtype=float)
    out = np.full(p.shape, np.nan)
    p_bottom = float(height_to_pressure(MIN_ALTITUDE))
    p_top86 = _layer_pressure(BASE_PRESSURES[-1], BASE_TEMPERATURES[-1],
                              GRADIENTS[-1], TOP_GEOPOTENTIAL - BASE_HEIGHTS[-2])
    with np.errstate(invalid="ignore"):
        low = np.isfinite(p) & (p > 0) & (p <= p_bottom) & (p >= p_top86)
        high = np.isfinite(p) & (p < p_top86) & (p >= TABLE_PRESSURES[-1])
    # layer: the highest base whose pressure is >= p (base pressures decrease)
    layer = np.clip(np.searchsorted(-BASE_PRESSURES, -p[low], side="right") - 1,
                    0, GRADIENTS.size - 1)
    pb, tb, grad = BASE_PRESSURES[layer], BASE_TEMPERATURES[layer], GRADIENTS[layer]
    with np.errstate(divide="ignore", invalid="ignore"):
        # inverse of eq (33a): H = Hb + (TM,b / L) [(P / Pb)^(-L / k) - 1]
        dh_grad = tb / grad * ((p[low] / pb) ** (-grad / _HYDROSTATIC) - 1)
        # inverse of eq (33b): H = Hb - (TM,b / k) ln(P / Pb)
        dh_iso = -tb / _HYDROSTATIC * np.log(p[low] / pb)
    h = BASE_HEIGHTS[layer] + np.where(grad == 0, dh_iso, dh_grad)
    out[low] = EARTH_RADIUS * h / (EARTH_RADIUS - h)  # eq (19)
    # above 86 km: Z linear in log(pressure) between table values (np.interp
    # needs increasing x, so both are reversed). Pressures between the table's
    # 86 km value and the equations' (1.2e-6 relative apart) are clamped by
    # np.interp to 86 km, deliberately

    out[high] = np.interp(np.log(p[high]), np.log(TABLE_PRESSURES[::-1]),
                          TABLE_HEIGHTS[::-1])
    return out[()] if out.ndim == 0 else out
