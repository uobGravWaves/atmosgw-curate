"""The shared state rules (docs/design/data_formats.md section 6.3), checked
on both loaders and on to_profiles with the same scenarios.

Synthetic files only (written to a temporary folder): limb "Misc" files and
AIRS granules (the writer of test_nadir_known.py). States and counts exact.

Scenarios: a gap inside what was measured (missing); beyond the measured
heights (no_observation); a physically impossible value (outside_physical_range,
NaN, kept with keep_outliers); qualified states only on values; every NaN
explained by exactly one "never there" state or a "removed" state; the
preprocessing attributes; to_profiles keeps each point's states.
"""

import warnings

import numpy as np
import pytest
from scipy.io import savemat
from test_nadir_known import D, seconds, write_granule

from atmosgw.curate.data import states, to_profiles, validate
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_limb_sounders, load_nadir_sounders

LEV = np.arange(10, 52, 2, dtype=float)  # limb file levels, km
GAP_PROF, GAP_LEV = 3, 10                # limb: NaN temperature at 30 km
BAD_PROF, BAD_LEV = 5, 6                 # limb: 50 K at 22 km
NADIR_Z = np.r_[0:61:3, 65:91:5].astype(float)


@pytest.fixture(scope="module")
def limb_files(tmp_path_factory):
    folder = tmp_path_factory.mktemp("conformance_limb")
    n_prof = 8
    alt = np.tile(LEV, (n_prof, 1))
    temp = 200 + 1.5 * alt
    temp[GAP_PROF, GAP_LEV] = np.nan
    temp[BAD_PROF, BAD_LEV] = 50.0
    temp[1, -1] = 50.0     # removed at the top of profile 1 ...
    temp[1, -2] = np.nan   # ... so the gap below it is still inside: missing
    savemat(folder / "conf_2007d001.mat", {
        "Alt": alt, "Temp": temp,
        "Lat": np.tile(np.linspace(-60, 60, n_prof)[:, None], (1, LEV.size)),
        "Lon": np.tile(np.linspace(-150, 150, n_prof)[:, None], (1, LEV.size)),
        "Time": np.tile(datenum(2007, 1, 1) + (np.arange(n_prof)[:, None] + 0.5)
                        / n_prof, (1, LEV.size)),
    })
    return folder


@pytest.fixture(scope="module")
def nadir_files(tmp_path_factory):
    folder = tmp_path_factory.mktemp("conformance_nadir")
    temp = np.broadcast_to(220 + 0.5 * NADIR_Z, (135, 90, 27)).astype(np.float32).copy()
    temp[100, 20, 4] = -999.0  # fill value at a scanned point: missing
    temp[60, 45, 8] = 50.0     # outside the physical range
    write_granule(folder, 53, seconds(53), temp)
    return folder


def limb(folder, **options):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_limb_sounders("2007-01-01", "Misc", misc_path=folder,
                                  misc_identifier="conf", **options)


def nadir(folder, **options):
    return load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=folder,
                               original_grid=True, **options)


def has(ds, name, state):
    return states.has_state(ds[f"{name}_flag"].values, state)


def all_outputs(limb_files, nadir_files):
    """Every output the rules apply to, with its main data variable."""
    swath = nadir(nadir_files)
    return [
        (limb(limb_files), "temperature"),
        (limb(limb_files, original_z=True), "temperature"),
        (limb(limb_files, keep_outliers=True, extrapolate=True), "temperature"),
        (swath, "temperature"),
        (nadir(nadir_files, keep_outliers=True), "temperature"),
        (load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=nadir_files,
                             height_scale=np.arange(0, 100.5, 3)), "temperature"),
        (to_profiles(swath), "temperature"),
    ]


# ------------------------------------------------------------------ rules


def test_gap_inside_measured_is_missing(limb_files, nadir_files):
    ds = limb(limb_files, original_z=True)
    assert has(ds, "temperature", "missing")[GAP_PROF, GAP_LEV]
    assert has(ds, "temperature", "missing")[1, -2]  # next to a removed point
    assert has(ds, "temperature", "missing").sum() == 2
    ds = limb(limb_files)  # interpolated: the heights around 30 km
    gap = has(ds, "temperature", "missing")[GAP_PROF]
    assert gap.any() and not has(ds, "temperature", "missing")[GAP_PROF - 1].any()
    sw = nadir(nadir_files)
    assert has(sw, "temperature", "missing")[100, 20, 4]
    assert has(sw, "temperature", "missing").sum() == 1


def test_beyond_measured_is_no_observation(limb_files, nadir_files):
    ds = limb(limb_files)  # height scale to 60 km, files to 50 km
    above = ds.altitude.values > LEV.max()
    flag = ds.temperature_flag.values
    assert np.all(states.has_state(flag, "no_observation")[above])
    assert not np.any(states.has_state(flag, "missing")[above])
    sw = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=nadir_files,
                             height_scale=np.arange(0, 100.5, 3))
    top = sw.altitude.values > NADIR_Z.max()
    flag = sw.temperature_flag.values[..., top]
    assert np.all(states.has_state(flag, "no_observation"))
    assert not np.any(states.has_state(flag, "missing"))


def test_impossible_value_removed_or_kept(limb_files, nadir_files):
    for removed, kept, where in (
            (limb(limb_files, original_z=True), limb(limb_files, original_z=True,
                                                     keep_outliers=True),
             (BAD_PROF, BAD_LEV)),
            (nadir(nadir_files), nadir(nadir_files, keep_outliers=True), (60, 45, 8))):
        assert has(removed, "temperature", "outside_physical_range")[where]
        assert np.isnan(removed.temperature.values[where])
        assert has(kept, "temperature", "outside_physical_range")[where]
        assert kept.temperature.values[where] == np.float32(50)


def test_removed_point_flagged_in_every_variable(limb_files):
    ds = limb(limb_files, original_z=True)
    for name in ("temperature", "pressure"):
        assert has(ds, name, "outside_physical_range")[BAD_PROF, BAD_LEV]
        assert np.isnan(ds[name].values[BAD_PROF, BAD_LEV])


def test_every_nan_explained_once_and_qualified_only_on_values(limb_files,
                                                                nadir_files):
    for ds, _ in all_outputs(limb_files, nadir_files):
        validate(ds, ds.attrs["format"])
        for name in ds.data_vars:
            if not name.endswith("_flag"):
                continue
            value = ds[name[:-5]].values
            flag = ds[name].values
            nan = np.isnan(value)
            never = (states.has_state(flag, "no_observation")
                     .astype(int) + states.has_state(flag, "missing").astype(int))
            assert np.all(never <= 1), name  # never both
            assert np.all(states.no_value_mask(flag)[nan]), name
            for qualified in ("derived", "extrapolated", "filled"):
                assert not np.any(states.has_state(flag, qualified) & nan), (
                    name, qualified)


def test_preprocessing_attributes(limb_files, nadir_files):
    for ds, _ in all_outputs(limb_files, nadir_files):
        assert ds.attrs["data"] in ("raw", "preprocessed", "planetary-wave-filtered")
        assert (ds.attrs["data"] == "preprocessed") == bool(ds.attrs["preprocessing"])
        assert ds.attrs["heights"] in ("interpolated", "original", "none")
    assert limb(limb_files).attrs["data"] == "raw"   # Misc files
    assert nadir(nadir_files).attrs["data"] == "raw"  # AIRS granule files


def test_nadir_extrapolated_only_outside_a_segments_own_coverage(tmp_path):
    # two segments (granule 54 missing); the second is narrower across track,
    # so the common cross axis (from the first) extends beyond it
    write_granule(tmp_path, 53, seconds(53))
    write_granule(tmp_path, 55, seconds(55), cross_scale=0.8)
    granules = [(D, 53), (D, 55)]
    plain = load_nadir_sounders(None, "AIRS", granules=granules, file_path=tmp_path)
    ext = load_nadir_sounders(None, "AIRS", granules=granules, file_path=tmp_path,
                              extrapolate=True)
    assert not has(plain, "temperature", "extrapolated").any()
    flagged = has(ext, "temperature", "extrapolated")
    second = ext.segment.values == 2
    assert flagged[second].any() and not flagged[~second].any()
    assert np.all(np.isfinite(ext.temperature.values[flagged]))
    # extrapolated points are exactly those NaN without extrapolate
    np.testing.assert_array_equal(flagged, np.isnan(plain.temperature.values)
                                  & np.isfinite(ext.temperature.values))


def test_to_profiles_keeps_each_points_states(nadir_files):
    sw = nadir(nadir_files)
    prof = to_profiles(sw)
    n_al, n_x, n_z = sw.temperature_flag.shape
    np.testing.assert_array_equal(prof.temperature_flag.values,
                                  sw.temperature_flag.values.reshape(n_al * n_x, n_z))
    assert prof.attrs["data"] == sw.attrs["data"]
    assert prof.attrs["preprocessing"] == sw.attrs["preprocessing"]
