"""Direct known-answer tests for atmosgw.curate.profiles.

(These steps are also checked against MATLAB through tests/test_limb_matlab.py.)
"""

import numpy as np
import pandas as pd
import pytest

from atmosgw.curate.profiles import (
    _interpolate_each,
    clip_to_height_range,
    concat_profiles,
    interpolate_profiles,
    mask_points,
    profiles_in_box,
    select_profiles,
)


def test_concat_pads_with_nan_and_makes_float64():
    a = {"x": np.ones((2, 3), dtype=np.float32), "t": pd.DataFrame({"c": [1, 2]})}
    b = {"x": np.full((1, 5, 2), 2.0), "t": pd.DataFrame({"c": [3]})}
    out = concat_profiles([a, b], ["x", "t"])
    assert out["x"].shape == (3, 5, 2) and out["x"].dtype == np.float64
    assert np.all(out["x"][:2, :3, 0] == 1) and np.all(np.isnan(out["x"][:2, 3:, :]))
    assert np.all(np.isnan(out["x"][:2, :3, 1]))
    assert list(out["t"]["c"]) == [1, 2, 3]


def test_concat_no_days():
    out = concat_profiles([], ["x"])
    assert out["x"].shape == (0, 0)


def test_select_profiles_skips_and_tables():
    data = {"x": np.arange(6.0).reshape(3, 2), "meta": np.array([[7.0]]),
            "t": pd.DataFrame({"c": [1, 2, 3]})}
    out = select_profiles(data, np.array([0, 2]), skip=["meta"])
    np.testing.assert_array_equal(out["x"], [[0, 1], [4, 5]])
    assert out["meta"].shape == (1, 1)
    assert list(out["t"]["c"]) == [1, 3]


def test_interpolate_uses_finite_first_heights_and_needs_three():
    alt = np.array([[10.0, 12.0, 12.0, np.nan, 14.0],   # repeat and NaN height
                    [10.0, np.nan, np.nan, np.nan, 12.0]])  # only two usable
    values = np.array([[1.0, 2.0, 99.0, 50.0, 3.0],
                       [1.0, 2.0, 3.0, 4.0, 5.0]])
    out = interpolate_profiles(alt, values, np.array([11.0, 13.0, 15.0]))
    np.testing.assert_array_equal(out[0], [1.5, 2.5, np.nan])
    assert np.all(np.isnan(out[1]))


def test_interpolate_3d_and_descending():
    alt = np.array([[16.0, 14.0, 12.0, 10.0]])
    values = np.stack([alt, 2 * alt], axis=2)
    out = interpolate_profiles(alt, values, np.array([11.0, 15.0]))
    np.testing.assert_allclose(out[0], [[11.0, 22.0], [15.0, 30.0]], rtol=1e-15)


@pytest.mark.parametrize("method", ["linear", "nearest", "next", "previous"])
@pytest.mark.parametrize("extrap", [None, True])
def test_vectorised_interpolation_identical_to_loop(method, extrap):
    # the vectorised path must give exactly what interp1, profile by profile,
    # gives: random heights in any order, repeated and NaN heights, NaN values,
    # constant stretches, profiles with too few levels, 3-D values, queries
    # outside the range and exactly on sample heights
    rng = np.random.default_rng(11)
    n_prof, n_lev = 300, 40
    alt = np.sort(rng.uniform(5, 70, (n_prof, n_lev)), axis=1)
    alt[::3] = alt[::3, ::-1]                       # some descending
    alt[1::7] = rng.permutation(alt[1::7].T).T      # some shuffled
    alt[rng.random((n_prof, n_lev)) < 0.05] = np.nan
    alt[5, 3] = alt[5, 4]                           # repeated height
    alt[9, 3:] = np.nan                             # only three usable
    alt[10, 2:] = np.nan                            # too few
    values = rng.normal(250, 20, (n_prof, n_lev, 2))
    values[rng.random(values.shape) < 0.03] = np.nan
    values[20:30, 10:20, :] = 230.0                 # constant stretches
    new_alt = np.r_[np.arange(0, 80, 0.5), alt[0, 5], alt[1, 7]]
    got = interpolate_profiles(alt, values, new_alt, method, extrap)
    want = _interpolate_each(alt, values, new_alt, method, extrap)
    np.testing.assert_array_equal(got, want)


def test_vectorised_interpolation_in_small_chunks(monkeypatch):
    import atmosgw.curate.profiles as profiles

    rng = np.random.default_rng(12)
    alt = np.sort(rng.uniform(5, 70, (50, 30)), axis=1)
    values = rng.normal(size=(50, 30))
    hs = np.arange(10, 60, 0.5)
    whole = interpolate_profiles(alt, values, hs)
    monkeypatch.setattr(profiles, "_CHUNK_ELEMENTS", 3000)  # a few rows per chunk
    np.testing.assert_array_equal(interpolate_profiles(alt, values, hs), whole)


def test_clip_to_height_range():
    alt = np.array([[10.0, 20.0, 30.0], [11.0, 21.0, 31.0]])
    data = {"Alt": alt, "T": alt + 100, "PW": np.stack([alt, alt], axis=2)}
    out = clip_to_height_range(data, alt, (15, 25))
    np.testing.assert_array_equal(out["T"], [[120.0], [121.0]])
    assert out["PW"].shape == (2, 1, 2)


def test_mask_points_all_slices():
    mask = np.array([[True, False], [False, False]])
    data = {"T": np.ones((2, 2)), "PW": np.ones((2, 2, 3))}
    out = mask_points(data, mask)
    assert np.isnan(out["T"][0, 0]) and np.isnan(out["T"]).sum() == 1
    assert np.all(np.isnan(out["PW"][0, 0])) and np.isnan(out["PW"]).sum() == 3


def test_profiles_in_box_any_height():
    lat = np.array([[-10.0, 5.0], [10.0, 20.0], [np.nan, np.nan]])
    lon = np.array([[0.0, 0.0], [0.0, 0.0], [np.nan, np.nan]])
    np.testing.assert_array_equal(profiles_in_box(lat, lon, (0, 90), (-180, 180)),
                                  [True, True, True])
    np.testing.assert_array_equal(profiles_in_box(lat, lon, (15, 90), (-180, 180)),
                                  [False, True, True])


def test_clip_and_mask_in_place_give_the_same_result():
    # copy=False (used by load_limb_sounders to save memory) must give exactly
    # the results of the default; it may change the input arrays
    from atmosgw.curate.profiles import clip_to_height_range, mask_points

    rng = np.random.default_rng(1)

    def data():
        alt = np.tile(np.arange(0.0, 10.0), (4, 1)) + rng.uniform(0, 0.5, (4, 1))
        return {"Alt": alt, "T": np.where(alt > 8, np.nan, 200 + alt),
                "T32": (200 + alt).astype(np.float32),
                "PW": np.stack([alt, 2 * alt], axis=2),
                "Flag": np.ones((4, 10), dtype=int),
                "Meta": np.array([[3.0]])}

    mask = np.zeros((4, 10), dtype=bool)
    mask[1, 3] = mask[3, 9] = True
    for height_range in ((2.0, 7.0), (-1.0, 99.0)):  # some levels kept / all kept
        rng = np.random.default_rng(1)
        a = clip_to_height_range(data(), data()["Alt"], height_range, skip=["Meta"])
        rng = np.random.default_rng(1)
        d = data()
        b = clip_to_height_range(d, d["Alt"], height_range, skip=["Meta"],
                                 copy=False)
        for name in a:
            np.testing.assert_array_equal(a[name], b[name])
            assert a[name].dtype == b[name].dtype
    rng = np.random.default_rng(1)
    a = mask_points(data(), mask, skip=["Meta"])
    rng = np.random.default_rng(1)
    b = mask_points(data(), mask, skip=["Meta"], copy=False)
    for name in a:
        np.testing.assert_array_equal(a[name], b[name])
        assert a[name].dtype == b[name].dtype
    assert np.isnan(b["T32"][1, 3]) and np.isnan(b["PW"][3, 9]).all()
    assert np.isnan(b["Flag"][1, 3]) and b["Flag"].dtype == np.float64
