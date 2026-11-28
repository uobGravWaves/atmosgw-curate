"""analysis.detrend and analysis.smooth: known answers and behaviour.

Synthetic inputs only (no files), so these run in CI. The comparison with
MATLAB is in tests/test_detrend_matlab.py.

Tolerances (each also stated where used):
- KNOWN = 1e-9 K: float64 known answers (a polynomial of degree <= order
  is its own background; a quadratic with gaps). The least-squares fits of
  values near 250 K round at about 1e-13 K; 1e-9 K is far above that and
  far below any real error (a wave of 0.1 K).
- F32 = 1e-4 K: Dataset-level outputs, which are float32 (spacing 1.5e-5 K
  near 250 K) and, for swaths, computed from float32 inputs (input rounding
  up to 7.6e-6 K, carried through the fit).
"""

import importlib
import json
import warnings

import numpy as np
import pytest
import xarray as xr
from test_data_formats import example_swath, example_time_height

from atmosgw.curate.analysis import (
    TropopauseMaskWarning,
    detrend,
    find_tropopause,
    smooth_array,
)
from atmosgw.curate.data import FORMATS, load_netcdf, states, to_netcdf, validate

DT = importlib.import_module("atmosgw.curate.analysis.detrend")
KNOWN = 1e-9  # K (see the module docstring)
F32 = 1e-4  # K (see the module docstring)
Z = np.arange(0, 60.5, 0.5)  # km: 121 levels at 0.5 km
PL = ("profile", "level")
SG = {"order": 2, "length_km": 18.0, "min_points": None}
POLY = {"order": 4, "dim": "cross_track", "abscissa": "index", "min_points": None}


# ------------------------------------------------------------- helpers


def example_profile(temperature, altitude, **attrs):
    """A valid profile Dataset of temperature [profile x level] (float64)."""
    temperature = np.array(temperature, dtype=np.float64)
    shape = temperature.shape
    altitude = np.broadcast_to(np.asarray(altitude, dtype=np.float64), shape).copy()
    flag = states.new_flags(shape)
    states.set_state(flag, np.isnan(temperature), "missing")
    ds = xr.Dataset({
        "latitude": (PL, np.zeros(shape), {"units": "degrees_north"}),
        "longitude": (PL, np.zeros(shape), {"units": "degrees_east"}),
        "time": (PL, np.full(shape, np.datetime64("2007-01-01T00:00", "us"))),
        "altitude": (PL, altitude, {"units": "km"}),
        "temperature": (PL, temperature, {
            "units": "K", "origin": "retrieved",
            "ancillary_variables": "temperature_flag"}),
        "temperature_flag": (PL, flag, states.flag_attrs()),
        "pressure": (PL, 1013.25 * np.exp(-np.nan_to_num(altitude) / 7.0), {
            "units": "hPa", "origin": "derived",
            "ancillary_variables": "pressure_flag"}),
        "pressure_flag": (PL, states.new_flags(shape), states.flag_attrs()),
    })
    _profile_attrs(ds, ["temperature_uncertainty", "pressure_uncertainty"])
    ds.attrs.update(attrs)
    return ds


def example_pw_profile(residual, waves, altitude):
    """A valid planetary-wave-filtered profile Dataset: temperature_residual
    [profile x level], temperature_pw [profile x level x mode]."""
    residual = np.array(residual, dtype=np.float64)
    waves = np.array(waves, dtype=np.float64)
    ds = example_profile(residual, altitude).drop_vars(["temperature",
                                                        "temperature_flag"])
    for name, values, dims in (("temperature_residual", residual, PL),
                               ("temperature_pw", waves, PL + ("mode",))):
        flag = states.new_flags(values.shape)
        states.set_state(flag, np.isnan(values), "missing")
        ds[name] = (dims, values, {"units": "K", "origin": "derived",
                                   "ancillary_variables": f"{name}_flag"})
        ds[f"{name}_flag"] = (dims, flag, states.flag_attrs())
    _profile_attrs(ds, ["temperature_residual_uncertainty",
                        "temperature_pw_uncertainty", "pressure_uncertainty"])
    ds.attrs["data"] = "planetary_wave_filtered"
    return ds


def _profile_attrs(ds, uncertainties):
    tracked = [n for n, v in FORMATS["profile"]["variables"].items()
               if v.get("tracked") and n not in ds.variables]
    ds.attrs = {a: "" for a in FORMATS["profile"]["attrs"]}
    ds.attrs.update(format="profile", format_version=3, instrument="SYNTH",
                    product="temperature", characterisation="SYNTH-placeholder/0",
                    availability=states.format_availability(
                        {**{n: "not_requested" for n in tracked},
                         **{n: "not_provided" for n in uncertainties}}))


def example_3d_swath(temperature, altitude=None, segment=None):
    """A valid 3-D swath Dataset of temperature [along x cross x level]."""
    temperature = np.asarray(temperature, dtype=np.float32)
    ds = example_swath(*temperature.shape)
    flag = states.new_flags(temperature.shape)
    states.set_state(flag, np.isnan(temperature), "missing")
    dims = ds["temperature"].dims
    ds["temperature"] = (dims, temperature, ds["temperature"].attrs)
    ds["temperature_flag"] = (dims, flag, states.flag_attrs())
    if altitude is not None:
        ds["altitude"] = (("level",), np.asarray(altitude, dtype=np.float64),
                          {"units": "km"})
    if segment is not None:
        ds["segment"] = (("along_track",), np.asarray(segment, dtype=np.float64))
    return ds


def example_2d_swath(bt):
    """A valid 2-D swath Dataset of brightness_temperature [along x cross x
    channel]."""
    bt = np.asarray(bt, dtype=np.float32)
    ds = example_swath(bt.shape[0], bt.shape[1], 1).drop_vars(
        ["temperature", "temperature_flag", "altitude"])
    dims = ("along_track", "cross_track", "channel")
    flag = states.new_flags(bt.shape)
    states.set_state(flag, np.isnan(bt), "missing")
    ds["brightness_temperature"] = (dims, bt, {
        "units": "K", "origin": "measured",
        "ancillary_variables": "brightness_temperature_flag"})
    ds["brightness_temperature_flag"] = (dims, flag, states.flag_attrs())
    listed = states.parse_availability(ds.attrs["availability"])
    listed.pop("temperature_uncertainty")
    listed.update(temperature="not_in_files",
                  brightness_temperature_uncertainty="not_provided")
    ds.attrs.update(product="2d", availability=states.format_availability(listed))
    return ds


def wave_profiles(n_profile=3, altitude=Z):
    """Temperature that rises with height plus a 1 K, 6 km wave: its lapse
    rate never exceeds 2 K/km, so it has no tropopause (no mask)."""
    phase = np.arange(n_profile)[:, None]
    return 220 + 0.5 * altitude + np.sin(2 * np.pi * altitude / 6 + phase)


def kinked_profile(z0=11.0, altitude=Z):
    """6.5 K/km from 288 K at 0 km up to z0, isothermal above."""
    return 288.0 - 6.5 * np.minimum(altitude, z0)


def record(out, var="temperature"):
    return json.loads(out[f"{var}_perturbation"].attrs["detrending"])


def polynomial_rows(shape, order, rng):
    """[along x cross x level] values, each cross-track row a different
    polynomial of degree `order` in the pixel number 1..N (about 250 K)."""
    x = np.arange(1, shape[1] + 1, dtype=np.float64)
    scale = (shape[1] / 2.0) ** np.arange(order + 1)
    coefficients = rng.normal(size=(shape[0], shape[2], order + 1)) / scale
    values = 250 + 5 * np.einsum("alk,ck->acl", coefficients,
                                 x[:, None] ** np.arange(order + 1))
    return values


# ---------------------------------------------------- Savitzky-Golay


@pytest.mark.parametrize("order", [2, 3])
def test_savgol_polynomial_is_its_own_background(order):
    # degree <= order: the fit is exact, edges included (they take the fit
    # to the first and last full frame)
    t = 200 + 3 * Z - 0.05 * Z ** 2 + (1e-3 * Z ** 3 if order == 3 else 0)
    t = np.vstack([t, t[::-1] + 10])
    alt = np.broadcast_to(Z, t.shape)
    bg, info = DT._savgol_step(alt, t, PL, {**SG, "order": order}, False)
    np.testing.assert_allclose(t - bg, 0, atol=KNOWN)  # (KNOWN: rounding)
    assert info["frame_levels"] == [37, 37] and info["regular_grid"]


@pytest.mark.parametrize("length, spacing, frame", [
    (18, 0.5, 37),  # 36, even: 37
    (18, -0.5, 35),  # descending: -36, even: -35 (MATLAB arithmetic)
    (7.5, 1.0, 9),  # 7.5 rounds away from zero to 8, even: 9
    (6.5, 1.0, 7),
    (-6.5, 1.0, 7),  # -6.5 rounds to -7
    (7.4, 1.0, 7),
    (7.6, 1.0, 9),
    (5, 0.5, 11),
])
def test_frame_arithmetic(length, spacing, frame):
    assert DT._frame_levels(length, spacing) == frame


def test_frame_from_the_grid():
    t = wave_profiles(2)
    _, info = DT._savgol_step(np.broadcast_to(Z, t.shape), t, PL, SG, False)
    assert info["frame_levels"] == [37, 37] and info["frame_km"] == 18.5
    down = Z[::-1]
    _, info = DT._savgol_step(np.broadcast_to(down, t.shape), t, PL, SG, False)
    assert info["frame_levels"] == [35, 35] and info["frame_km"] == 17.5


def test_irregular_grid_gets_a_frame_per_profile():
    alt = np.vstack([0.5 * np.arange(121), 0.25 * np.arange(121)])
    t = wave_profiles(2, alt)
    bg, info = DT._savgol_step(alt, t, PL, SG, False)
    assert not info["regular_grid"] and info["frame_levels"] == [37, 73]
    for p in range(2):  # each profile as if alone (then regular)
        alone, one = DT._savgol_step(alt[p:p + 1], t[p:p + 1], PL, SG, False)
        assert one["regular_grid"]
        np.testing.assert_array_equal(bg[p], alone[0])


@pytest.mark.parametrize("n_levels, options", [
    (20, SG),  # 37 levels > 20
    (121, {**SG, "length_km": 1.0, "order": 3}),  # 3 levels, not > order 3
])
def test_unsuitable_frame_raises(n_levels, options):
    t = wave_profiles(1, Z[:n_levels])
    with pytest.raises(ValueError, match="frame"):
        DT._savgol_step(Z[None, :n_levels], t, PL, options, False)


def test_one_nan_removes_half_a_frame_either_side():
    t = wave_profiles(1)
    t[0, 60] = np.nan
    bg, _ = DT._savgol_step(Z[None, :], t, PL, SG, False)
    lost = np.abs(np.arange(121) - 60) <= 18  # half of 37
    np.testing.assert_array_equal(np.isnan(bg[0]), lost)


# ------------------------------------------------------------ skip_nan


def test_savgol_skip_nan_same_where_nothing_missing():
    t = wave_profiles(3)
    alt = np.broadcast_to(Z, t.shape)
    a, _ = DT._savgol_step(alt, t, PL, SG, False)
    b, _ = DT._savgol_step(alt, t, PL, SG, True)
    np.testing.assert_array_equal(a, b)


def test_savgol_skip_nan_only_the_missing_point():
    t = wave_profiles(2)
    t[0, 60] = np.nan
    alt = np.broadcast_to(Z, t.shape)
    bg, _ = DT._savgol_step(alt, t, PL, SG, True)
    expected = np.zeros(t.shape, bool)
    expected[0, 60] = True
    np.testing.assert_array_equal(np.isnan(bg), expected)
    plain, _ = DT._savgol_step(alt, t, PL, SG, False)
    np.testing.assert_array_equal(bg[1], plain[1])  # (the full profile untouched)


def test_savgol_skip_nan_recovers_a_quadratic_with_gaps():
    t = (200 + 3 * Z - 0.05 * Z ** 2)[None, :].copy()
    gaps = [0, 1, 17, 40, 41, 42, 60, 100, 119, 120]
    t[0, gaps] = np.nan
    bg, _ = DT._savgol_step(Z[None, :], t, PL, SG, True)
    np.testing.assert_array_equal(np.isnan(bg), np.isnan(t))
    np.testing.assert_allclose(bg, t, atol=KNOWN)  # (KNOWN: rounding)


def test_savgol_min_points():
    t = wave_profiles(1)
    t[0, 60] = np.nan
    alt = Z[None, :]
    plain, _ = DT._savgol_step(alt, t, PL, SG, False)
    # every frame point required: the same NaN as without skip_nan
    full, _ = DT._savgol_step(alt, t, PL, {**SG, "min_points": 37}, True)
    np.testing.assert_array_equal(np.isnan(full), np.isnan(plain))
    # 36 of 37 enough: only the missing point
    most, _ = DT._savgol_step(alt, t, PL, {**SG, "min_points": 36}, True)
    assert np.flatnonzero(np.isnan(most[0])).tolist() == [60]
    # default, 19 of 37 (just over half): beside a one-sided gap at most 18
    # are missing, so every value present is computed ...
    t[0, 51:70] = np.nan
    gap, _ = DT._savgol_step(alt, t, PL, SG, True)
    np.testing.assert_array_equal(np.isnan(gap), np.isnan(t))
    # ... but a lone value between two gaps (1 of 37 present) is not
    t[0, 30:51] = np.nan
    t[0, 48] = 230.0
    gap, _ = DT._savgol_step(alt, t, PL, SG, True)
    assert np.isnan(gap[0, 48])


def test_polynomial_skip_nan_rows():
    rng = np.random.default_rng(1)
    t = polynomial_rows((3, 11, 2), 4, rng)
    t[0, 4, 0] = np.nan  # one missing pixel
    t[1, :6, 1] = np.nan  # 6 of 11 missing: 5 left, fewer than the 6 needed
    dims = ("along_track", "cross_track", "level")
    plain, _ = DT._polynomial_step(None, t, dims, POLY, False)
    assert np.isnan(plain[0, :, 0]).all() and np.isnan(plain[1, :, 1]).all()
    assert np.isfinite(plain[2]).all()
    gappy, info = DT._polynomial_step(None, t, dims, POLY, True)
    assert np.flatnonzero(np.isnan(gappy[0, :, 0])).tolist() == [4]
    np.testing.assert_allclose(gappy[0, :, 0], t[0, :, 0], atol=KNOWN)
    assert np.isnan(gappy[1, :, 1]).all()  # too few points: stays NaN
    assert info["rows_with_gaps_fitted"] == 1
    finite = np.isfinite(plain)
    np.testing.assert_array_equal(gappy[finite], plain[finite])
    # min_points: 11 required, so the one-gap row stays NaN too
    strict, _ = DT._polynomial_step(None, t, dims, {**POLY, "min_points": 11}, True)
    assert np.isnan(strict[0, :, 0]).all()


# ---------------------------------------------------------- polynomial


def test_polynomial_rows_of_degree_up_to_order_are_their_own_background():
    rng = np.random.default_rng(2)
    dims = ("along_track", "cross_track", "level")
    for degree in range(5):
        t = polynomial_rows((4, 11, 3), degree, rng)
        bg, info = DT._polynomial_step(None, t, dims, POLY, False)
        np.testing.assert_allclose(t - bg, 0, atol=KNOWN)  # (KNOWN: rounding)
    assert info["points"] == 11


def test_polynomial_equals_a_least_squares_fit_per_row():
    rng = np.random.default_rng(3)
    t = 250 + rng.normal(size=(4, 11, 3))
    bg, _ = DT._polynomial_step(None, t, ("along_track", "cross_track", "level"),
                                POLY, False)
    x = np.arange(1, 12, dtype=float)
    for a in range(4):
        for level in range(3):
            c = np.polynomial.polynomial.polyfit(x, t[a, :, level], 4)
            np.testing.assert_allclose(bg[a, :, level],
                                       np.polynomial.polynomial.polyval(x, c),
                                       atol=KNOWN)


def test_polynomial_abscissa_coordinate():
    x = 10 * np.arange(11, dtype=float) ** 1.5  # km, not linear in pixel number
    t = (250 + 0.05 * x - 2e-4 * x ** 2)[None, :, None] * np.ones((2, 1, 3))
    dims = ("along_track", "cross_track", "level")
    coord = xr.Dataset({"distance_cross": (("cross_track",), x)})
    bg, _ = DT._polynomial_step(coord, t, dims,
                                {**POLY, "order": 2, "abscissa": "distance_cross"},
                                False)
    np.testing.assert_allclose(t - bg, 0, atol=KNOWN)
    bg, _ = DT._polynomial_step(coord, t, dims, {**POLY, "order": 2}, False)
    assert np.max(np.abs(t - bg)) > 0.01  # (in pixel number it is not quadratic)
    with pytest.raises(ValueError, match="abscissa"):
        DT._polynomial_step(coord, t, dims, {**POLY, "abscissa": "distance_along"},
                            False)


def test_polynomial_any_dim():
    rng = np.random.default_rng(4)
    t = np.moveaxis(polynomial_rows((3, 9, 2), 4, rng), 1, 0)  # degree 4 along axis 0
    dims = ("along_track", "cross_track", "level")
    bg, _ = DT._polynomial_step(None, t, dims, {**POLY, "dim": "along_track"}, False)
    np.testing.assert_allclose(t - bg, 0, atol=KNOWN)
    # the same as moving that dimension to cross_track
    s = 250 + rng.normal(size=(9, 3, 2))
    a, _ = DT._polynomial_step(None, s, dims, {**POLY, "dim": "along_track"}, False)
    b, _ = DT._polynomial_step(None, np.moveaxis(s, 0, 1), dims, POLY, False)
    np.testing.assert_allclose(a, np.moveaxis(b, 1, 0), atol=KNOWN)


def test_polynomial_per_channel_in_a_2d_swath():
    rng = np.random.default_rng(5)
    bt = polynomial_rows((5, 11, 3), 4, rng)  # [along x cross x channel]
    ds = example_2d_swath(bt)
    out = detrend(ds)
    assert record(out, "brightness_temperature")["steps"][0]["method"] == "polynomial"
    # (F32: float32 input and output)
    np.testing.assert_allclose(out["brightness_temperature_perturbation"], 0,
                               atol=F32)
    validate(out, "swath")


def test_too_few_points_for_the_polynomial():
    ds = example_3d_swath(np.full((2, 4, 3), 250.0))
    with pytest.raises(ValueError, match="needs more than"):
        detrend(ds)


# ----------------------------------------------------------- hindley23


def _pw_inputs():
    rng = np.random.default_rng(6)
    residual = rng.normal(size=(3, 121))
    waves = 220 + rng.normal(size=(3, 121, 4))
    waves[0, 5, :] = np.nan  # every mode missing: background NaN
    waves[1, 7, 1:] = np.nan  # some modes missing: the sum of the others
    residual[2, 9] = np.nan
    return residual, waves


def test_hindley23():
    residual, waves = _pw_inputs()
    ds = example_pw_profile(residual, waves, Z)
    out = detrend(ds)  # (the default for planetary-wave data)
    assert record(out)["steps"][0]["method"] == "hindley23"
    np.testing.assert_array_equal(out["temperature_perturbation"],
                                  residual.astype(np.float32))
    expected = np.nansum(waves, axis=-1)
    expected[0, 5] = np.nan
    np.testing.assert_array_equal(out["temperature_background"],
                                  expected.astype(np.float32))
    flag = out["temperature_background_flag"].values
    assert states.has_state(flag, "not_computed")[0, 5]
    assert states.has_state(flag, "not_computed").sum() == 1
    # the residual is missing at [2, 9] but the waves are not: the background
    # has a value there, so it carries no "no value" state; the perturbation does
    assert np.isfinite(out["temperature_background"].values[2, 9])
    assert not states.no_value_mask(flag)[2, 9]
    assert states.has_state(out["temperature_perturbation_flag"].values,
                            "missing")[2, 9]
    assert "tropopause_altitude" not in out  # (no vertical filter, no mask)
    validate(out, "profile")


def test_hindley23_rules():
    residual, waves = _pw_inputs()
    pw = example_pw_profile(residual, waves, Z)
    raw = example_profile(wave_profiles(), Z)
    with pytest.raises(ValueError, match="hindley_pws"):
        detrend(raw, method="hindley23")
    with pytest.raises(ValueError, match="hindley23"):
        detrend(pw, method="savgol")
    with pytest.raises(ValueError, match="first step"):
        detrend(pw, steps=["savgol", "hindley23"])
    out = detrend(pw, steps=["hindley23", "savgol"], mask_tropopause=False)
    # savgol on the residual; backgrounds add
    bg, _ = DT._savgol_step(np.broadcast_to(Z, residual.shape), residual, PL, SG,
                            False)
    total = np.nansum(waves, axis=-1)
    total[0, 5] = np.nan
    np.testing.assert_allclose(out["temperature_perturbation"], residual - bg,
                               atol=F32)
    np.testing.assert_allclose(out["temperature_background"], total + bg, atol=F32)
    validate(out, "profile")


# ------------------------------------------------------------ smoothing


def test_smooth_constant_unchanged():
    x = np.full((6, 7, 5), 250.0)
    for kernel in ("box", "gaussian"):
        np.testing.assert_allclose(smooth_array(x, (3, 5, 1), kernel), x,
                                   rtol=0, atol=1e-12)  # (summing 1/3s: 1e-13)


def test_box_mean_known_answer():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    # edges repeated: (1 + 1 + 2) / 3 and (4 + 5 + 5) / 3
    np.testing.assert_allclose(smooth_array(x, 3), [4 / 3, 2, 3, 4, 14 / 3],
                               atol=1e-14)
    y = np.arange(12, dtype=float).reshape(3, 4)
    got = smooth_array(y, (3, 1))
    expected = np.vstack([(2 * y[0] + y[1]) / 3, y[1], (y[1] + 2 * y[2]) / 3])
    np.testing.assert_allclose(got, expected, atol=1e-14)
    np.testing.assert_array_equal(smooth_array(y, (1, 1)), y)


@pytest.mark.parametrize("size, sigma", [(5, 0.65), (5, 1.0), (7, 0.65), (3, 0.3)])
def test_gaussian_weights(size, sigma):
    half = (size - 1) // 2
    impulse = np.zeros(2 * size + 1)
    impulse[size] = 1.0
    got = smooth_array(impulse, size, "gaussian", sigma)[size - half:size + half + 1]
    k = np.arange(-half, half + 1)
    expected = np.exp(-k ** 2 / (2 * (sigma * half) ** 2))  # std: sigma * half
    np.testing.assert_allclose(got, expected / expected.sum(), atol=1e-15)
    assert abs(got.sum() - 1) < 1e-14


@pytest.mark.parametrize("size", [2, (3, 4), 0, 1.5])
def test_bad_size_rejected(size):
    with pytest.raises(ValueError, match="odd"):
        smooth_array(np.zeros((5, 5)), size)


def test_smooth_nan():
    x = np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0, 7.0])
    np.testing.assert_array_equal(np.isnan(smooth_array(x, 3)),
                                  [False, True, True, True, False, False, False])
    got = smooth_array(x, 3, skip_nan=True)
    np.testing.assert_array_equal(np.isnan(got), np.isnan(x))
    np.testing.assert_allclose(got[[0, 1, 3, 4]], [4 / 3, 1.5, 4.5, 5], atol=1e-14)


def test_swath_smoothing_stays_within_segments():
    t = np.where(np.arange(6)[:, None, None] < 3, 250.0, 260.0) * np.ones((6, 9, 2))
    steps = [("smooth", {"size": {"along_track": 3}}), "polynomial"]
    out = detrend(example_3d_swath(t, segment=[1, 1, 1, 2, 2, 2]), steps=steps,
                  keep_steps=True)
    np.testing.assert_allclose(out["temperature_perturbation_step1"], t, atol=F32)
    mixed = detrend(example_3d_swath(t), steps=steps, keep_steps=True)
    assert np.max(np.abs(mixed["temperature_perturbation_step1"] - t)) > 3


def test_smooth_needs_size():
    with pytest.raises(TypeError, match="size"):
        detrend(example_3d_swath(np.full((3, 9, 2), 250.0)),
                steps=["smooth", "polynomial"])


# --------------------------------------------------------------- chains


def _random_swath(seed=7, segment=None):
    rng = np.random.default_rng(seed)
    t = (250 + rng.normal(size=(8, 11, 4))).astype(np.float32)
    return example_3d_swath(t, segment=segment), t.astype(np.float64)


def test_chain_equals_steps_by_hand():
    ds, t = _random_swath()
    dims = ds["temperature"].dims
    out = detrend(ds, steps=[("smooth", {"size": {"cross_track": 3,
                                                  "along_track": 3}}),
                             "polynomial"], keep_steps=True)
    smoothed = smooth_array(t, (3, 3, 1))
    bg, _ = DT._polynomial_step(ds, smoothed, dims, POLY, False)
    # (F32: float32 outputs; here the input is the same float32 values)
    np.testing.assert_allclose(out["temperature_perturbation"], smoothed - bg,
                               atol=F32)
    np.testing.assert_allclose(out["temperature_background"], bg, atol=F32)
    np.testing.assert_allclose(out["temperature_perturbation_step1"], smoothed,
                               atol=F32)
    np.testing.assert_array_equal(out["temperature_perturbation_step2"],
                                  out["temperature_perturbation"])
    steps = record(out)["steps"]
    assert [s["method"] for s in steps] == ["smooth", "polynomial"]
    assert steps[0]["size_by_dimension"] == {"along_track": 3, "cross_track": 3,
                                             "level": 1}


def test_chain_backgrounds_add():
    ds, t = _random_swath()
    dims = ds["temperature"].dims
    out = detrend(ds, steps=["polynomial", ("polynomial", {"dim": "along_track"})])
    bg1, _ = DT._polynomial_step(ds, t, dims, POLY, False)
    bg2, _ = DT._polynomial_step(ds, t - bg1, dims, {**POLY, "dim": "along_track"},
                                 False)
    np.testing.assert_allclose(out["temperature_background"], bg1 + bg2, atol=F32)
    np.testing.assert_allclose(out["temperature_perturbation"], t - bg1 - bg2,
                               atol=F32)
    assert "temperature_perturbation_step1" not in out  # (keep_steps False)


def test_chain_errors():
    ds, _ = _random_swath()
    with pytest.raises(ValueError, match="no background"):
        detrend(ds, steps=[("smooth", {"size": 3})])
    with pytest.raises(TypeError, match="not both"):
        detrend(ds, method="polynomial", steps=["polynomial"])
    with pytest.raises(TypeError, match="not both"):
        detrend(ds, steps=["polynomial"], order=3)
    with pytest.raises(TypeError, match="no option"):
        detrend(ds, method="polynomial", frame=3)
    with pytest.raises(TypeError, match="no option"):
        detrend(ds, steps=[("polynomial", {"length_km": 3})])
    with pytest.raises(ValueError, match="Unknown method"):
        detrend(ds, method="fft")
    with pytest.raises(ValueError, match="empty"):
        detrend(ds, steps=[])


def test_aliases():
    ds, t = _random_swath()
    a = detrend(ds, steps=[("boxcar", {"size": {"cross_track": 3}}), "polynomial"])
    b = detrend(ds, steps=[("smooth", {"size": {"cross_track": 3}}), "polynomial"])
    xr.testing.assert_identical(a, b)
    g = detrend(ds, steps=[("gaussian", {"size": {"cross_track": 5}}), "polynomial"])
    assert record(g)["steps"][0]["kernel"] == "gaussian"


# ------------------------------------------------------ defaults by format


def _method(out, var="temperature"):
    return record(out, var)["steps"][0]


def test_defaults_by_format():
    profile = detrend(example_profile(wave_profiles(), Z))
    step = _method(profile)
    assert step["method"] == "savgol" and step["length_km"] == 18.0
    assert step["order"] == 2

    th = example_time_height(n_time=4, n_level=Z.size)
    th = th.assign_coords(altitude=(("level",), Z, {"units": "km"}))
    th["temperature"].values[1:] = wave_profiles(3)
    out = detrend(th)
    assert _method(out)["method"] == "savgol" and _method(out)["length_km"] == 18.0
    validate(out, "time_height")

    swath = detrend(_random_swath()[0])
    step = _method(swath)
    assert step["method"] == "polynomial" and step["order"] == 4
    assert step["dim"] == "cross_track" and step["abscissa"] == "index"

    rng = np.random.default_rng(8)
    two_d = detrend(example_2d_swath(250 + rng.normal(size=(4, 11, 2))))
    assert "brightness_temperature_perturbation" in two_d
    assert _method(two_d, "brightness_temperature")["method"] == "polynomial"

    sonde = detrend(example_profile(wave_profiles(), Z,
                                    converted_from="trajectory"))
    step = _method(sonde)
    assert step["length_km"] == 5.0 and step["frame_levels"] == [11, 11]
    # (an explicit length is kept)
    sonde = detrend(example_profile(wave_profiles(), Z, converted_from="trajectory"),
                    length_km=10.0)
    assert _method(sonde)["length_km"] == 10.0


def test_trajectory_refused():
    ds = xr.Dataset(attrs={"format": "trajectory"})
    with pytest.raises(ValueError, match="to_profiles"):
        detrend(ds)
    with pytest.raises(ValueError, match="format"):
        detrend(xr.Dataset())


def test_hlos_wind_needs_var():
    ds = example_profile(wave_profiles(), Z, product="hlos_wind")
    ds = ds.rename({"temperature": "hlos_wind", "temperature_flag": "hlos_wind_flag"})
    ds["hlos_wind"].attrs.update(units="m s-1",
                                 ancillary_variables="hlos_wind_flag")
    with pytest.raises(ValueError, match="var="):
        detrend(ds)
    out = detrend(ds, var="hlos_wind", mask_tropopause=False)
    assert out["hlos_wind_perturbation"].attrs["units"] == "m s-1"
    with pytest.raises(ValueError, match="not in the Dataset"):
        detrend(ds, var="eastward_wind")


# ---------------------------------------------------------------- outputs


def _flagged_profile():
    t = wave_profiles(3)
    t[0, 60] = np.nan
    ds = example_profile(t, Z)
    states.set_state(ds["temperature_flag"].values, np.eye(3, 121, 10, bool),
                     "a_priori_dominated")
    return ds


def test_outputs():
    ds = _flagged_profile()
    before = ds.copy(deep=True)
    out = detrend(ds)
    xr.testing.assert_identical(ds, before)  # (the input is not altered)
    for name in ("temperature_perturbation", "temperature_background"):
        v = out[name]
        assert v.dims == PL and v.dtype == np.float32
        assert v.attrs["units"] == "K" and v.attrs["origin"] == "derived"
        assert v.attrs["ancillary_variables"] == f"{name}_flag"
        assert out[f"{name}_flag"].dtype == states.FLAG_DTYPE
    total = (out["temperature_perturbation"].values.astype(float)
             + out["temperature_background"].values)
    finite = np.isfinite(total)
    np.testing.assert_allclose(total[finite], ds["temperature"].values[finite],
                               atol=F32)
    info = record(out)
    assert info["variable"] == "temperature" and info["skip_nan"] is False
    step = info["steps"][0]
    assert set(DT.DEFAULTS["savgol"]) <= set(step)
    assert step["frame_levels"] == [37, 37] and step["frame_km"] == 18.5
    assert step["regular_grid"] is True
    assert info["tropopause_mask"]["half_width_km"] == 9.25
    assert info["tropopause_mask"]["values_masked"] == 0
    listed = states.parse_availability(out.attrs["availability"])
    assert listed["temperature_perturbation_uncertainty"] == "not_provided"
    assert listed["temperature_background_uncertainty"] == "not_provided"
    assert "temperature_perturbation" in out.attrs["state_counts"]
    validate(out, "profile")


def test_output_flags():
    ds = _flagged_profile()
    out = detrend(ds)
    given = ds["temperature_flag"].values
    lost = np.zeros((3, 121), bool)
    lost[0] = np.abs(np.arange(121) - 60) <= 18
    for name in ("temperature_perturbation", "temperature_background"):
        values, flag = out[name].values, out[f"{name}_flag"].values
        np.testing.assert_array_equal(np.isnan(values), lost)
        # input states carried
        np.testing.assert_array_equal(flag & given, given)
        np.testing.assert_array_equal(states.has_state(flag, "derived"),
                                      np.isfinite(values))
        np.testing.assert_array_equal(
            states.has_state(flag, "not_computed"),
            np.isnan(values) & ~states.no_value_mask(given))
        assert not states.has_state(flag, "not_computed")[0, 60]  # (missing)
        assert states.has_state(flag, "a_priori_dominated")[1, 11]


def test_validate_each_format():
    validate(detrend(_flagged_profile()), "profile")
    validate(detrend(_random_swath()[0]), "swath")
    th = example_time_height(n_time=4, n_level=Z.size)
    th = th.assign_coords(altitude=(("level",), Z, {"units": "km"}))
    th["temperature"].values[1:] = wave_profiles(3)
    validate(detrend(th), "time_height")


def test_existing_outputs_need_overwrite():
    ds, _ = _random_swath()
    out = detrend(ds)
    with pytest.raises(ValueError, match="overwrite=True"):
        detrend(out)
    again = detrend(out, steps=[("polynomial", {"order": 2})], overwrite=True)
    assert _method(again)["order"] == 2
    assert again.attrs["state_counts"].count("temperature_perturbation:") == 1
    validate(again, "swath")
    # the producer's own perturbation in a 2-D swath
    rng = np.random.default_rng(9)
    two_d = example_2d_swath(250 + rng.normal(size=(4, 11, 2)))
    first = detrend(two_d)
    with pytest.raises(ValueError, match="already in the Dataset"):
        detrend(first)
    validate(detrend(first, overwrite=True), "swath")


def test_netcdf_round_trip(tmp_path):
    out = detrend(_flagged_profile())
    to_netcdf(out, tmp_path / "detrended.nc")
    back = load_netcdf(tmp_path / "detrended.nc")
    for name in ("temperature_perturbation", "temperature_background",
                 "temperature_perturbation_flag", "tropopause_altitude"):
        np.testing.assert_array_equal(back[name].values, out[name].values)
        assert back[name].dtype == out[name].dtype
    assert (back["temperature_perturbation"].attrs["detrending"]
            == out["temperature_perturbation"].attrs["detrending"])


# ------------------------------------------------------------ tropopause


@pytest.mark.parametrize("z0", [11.0, 11.3, 16.7])
def test_find_tropopause_kink(z0):
    t = kinked_profile(z0)
    found = find_tropopause(Z, t)
    assert z0 <= found[0] <= z0 + 0.5  # (within the grid spacing, at or above)
    # descending level order: the same
    assert find_tropopause(Z[::-1], t[::-1])[0] == found[0]


def test_find_tropopause_fine_noisy_profiles():
    rng = np.random.default_rng(10)
    z = np.arange(0, 30.0001, 0.05)
    z0 = np.array([9.0, 12.4, 15.0])
    t = 288 - 6.5 * np.minimum(z[None, :], z0[:, None]) + rng.normal(
        0, 0.05, (3, z.size))
    found = find_tropopause(z, t)
    # (averaged into 0.25 km layers first: within one layer of the kink)
    assert np.all(np.abs(found - z0) <= 0.25 + 0.05), found


def test_find_tropopause_none():
    z = np.arange(15, 40.01, 0.5)  # starts above the tropopause
    assert np.isnan(find_tropopause(z, np.full(z.size, 216.5)))[0]
    z = np.arange(0, 12.01, 0.5)  # top within 2 km of the kink at 11 km
    assert np.isnan(find_tropopause(z, kinked_profile(11.0, z)))[0]
    assert np.isnan(find_tropopause(Z, wave_profiles(1)))[0]
    found = find_tropopause(np.full((2, 3), np.nan), np.full((2, 3), 250.0))
    assert np.isnan(found).all()


def _kinked_dataset():
    t = np.vstack([kinked_profile(11.0), wave_profiles(1)[0]])
    return example_profile(t, Z)


def test_tropopause_mask():
    ds = _kinked_dataset()
    with pytest.warns(TropopauseMaskWarning, match="TROPOPAUSE MASK"):
        out = detrend(ds)
    free = detrend(ds, mask_tropopause=False)
    near = np.zeros((2, 121), bool)
    near[0] = np.abs(Z - 11.0) <= 9.25  # half the 18.5 km frame
    pert, plain = out["temperature_perturbation"].values, free[
        "temperature_perturbation"].values
    np.testing.assert_array_equal(np.isnan(pert), near | np.isnan(plain))
    np.testing.assert_array_equal(pert[~near], plain[~near])
    assert np.isnan(out["temperature_background"].values[near]).all()
    flag = out["temperature_perturbation_flag"].values
    np.testing.assert_array_equal(states.has_state(flag, "not_computed"), near)
    np.testing.assert_array_equal(out["tropopause_altitude"].values, [11.0, np.nan])
    assert out["tropopause_altitude"].dims == ("profile",)
    note = record(out)["tropopause_mask"]
    assert note["values_masked"] == int(near.sum()) and note["profiles_masked"] == 1
    assert "tropopause_altitude" not in free
    validate(out, "profile")


def test_tropopause_options():
    ds = _kinked_dataset()
    with pytest.warns(TropopauseMaskWarning):
        given = detrend(ds, tropopause_km=15.0)
    expected = np.abs(Z - 15.0) <= 9.25
    np.testing.assert_array_equal(np.isnan(given["temperature_perturbation"][0]),
                                  expected)
    np.testing.assert_array_equal(np.isnan(given["temperature_perturbation"][1]),
                                  expected)  # (a number applies to every profile)
    with pytest.warns(TropopauseMaskWarning):
        each = detrend(ds, tropopause_km=[15.0, np.nan])
    assert not np.isnan(each["temperature_perturbation"][1]).any()
    with pytest.warns(TropopauseMaskWarning):
        narrow = detrend(ds, mask_half_width_km=2.0)
    np.testing.assert_array_equal(np.isnan(narrow["temperature_perturbation"][0]),
                                  np.abs(Z - 11.0) <= 2.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error", TropopauseMaskWarning)
        off = detrend(ds, mask_tropopause=False)
    assert np.isfinite(off["temperature_perturbation"]).all()


def test_no_warning_without_a_tropopause():
    ds = example_profile(wave_profiles(), Z)
    with warnings.catch_warnings():
        warnings.simplefilter("error", TropopauseMaskWarning)
        out = detrend(ds)
    assert np.isnan(out["tropopause_altitude"]).all()
    assert np.isfinite(out["temperature_perturbation"]).all()


def test_profile_without_values():
    # (e.g. a padding profile, or one removed entirely by quality control)
    t = np.vstack([kinked_profile(11.0), np.full(Z.size, np.nan)])
    ds = example_profile(t, Z)
    with pytest.warns(TropopauseMaskWarning):
        out = detrend(ds)
    np.testing.assert_array_equal(out["tropopause_altitude"].values, [11.0, np.nan])
    assert np.isnan(out["temperature_perturbation"][1]).all()
    flag = out["temperature_perturbation_flag"].values[1]
    assert not states.has_state(flag, "not_computed").any()
    validate(out, "profile")


def test_tropopause_mask_only_for_vertical_filters():
    with warnings.catch_warnings():
        warnings.simplefilter("error", TropopauseMaskWarning)
        out = detrend(_random_swath()[0])
    assert "tropopause_altitude" not in out


def test_find_tropopause_needs_the_troposphere():
    """Data that start in the stratosphere have no tropopause, even where the
    winter polar stratosphere cools with height (HIRDLS, 70 N, 2007-01-01:
    3.8 K/km at 18-19 km, then nearly isothermal, passed for one before)."""
    z = np.arange(18.0, 60.5, 0.5)
    t = np.where(z < 21, 209.5 - 3.2 * (z - 18.0), 199.9 + 0.2 * (z - 21.0))
    assert np.isnan(find_tropopause(z[None, :], t[None, :])[0])
    # the same shape lower down, reaching the troposphere, is found
    low = np.arange(2.0, 44.5, 0.5)
    assert find_tropopause(low[None, :], t[None, :])[0] == pytest.approx(5.0, abs=0.5)
    # reach_km moves the rule
    assert np.isfinite(find_tropopause(z[None, :], t[None, :], reach_km=18.0,
                                       search_km=(18.0, 22.0))[0])
