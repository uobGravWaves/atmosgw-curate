"""Regridded swath geometry (Q14 in docs/quirks.md): positions of the values,
spacing of partial granules, scans that cannot be placed, segments after a
box, zero-weight corners.

Synthetic AIRS granule files (tests/test_nadir_known.write_granule), so these
run in CI; one real-data check (skipped if absent).
"""

import os
import warnings

import netCDF4
import numpy as np
import pytest
from test_nadir_known import D, seconds, write_granule

from atmosgw.curate import swath
from atmosgw.curate.data import states, validate
from atmosgw.curate.obs import load_nadir_sounders
from atmosgw.curate.paths import DataDirNotSetError, dataset_dir


def _granule(root, g, **options):
    write_granule(root, g, seconds(g), **options)
    return os.path.join(root, "2007", "001", f"airs_2007_001_{g:03d}.nc")


def _set(path, name, function):
    """Replace a variable of a granule file by function(lat, lon, value)."""
    with netCDF4.Dataset(path, "a") as f:
        f.set_auto_maskandscale(False)
        lat, lon = f["l1_lat"][:], f["l1_lon"][:]
        f[name][:] = function(lat, lon, f[name][:])


def test_positions_are_those_of_the_values(tmp_path):
    # temperature = 200 K + latitude (degrees): wherever the regridded value
    # is, its latitude must be (value - 200), also at the swath edges (the
    # regular-grid positions were up to 70 km off there on real data)
    path = _granule(tmp_path, 53)
    _set(path, "ret_temp", lambda lat, lon, t: (200 + lat)[..., None] + 0 * t)
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             original_z=True)
    t = ds.temperature.values[..., 0].astype(float)
    ok = np.isfinite(t)
    assert ok.mean() > 0.95
    # (latitude is interpolated as a direction on the sphere, temperature
    # linearly: up to 1.2e-3 degrees, 130 m, apart at 80 N over 20 km cells)
    np.testing.assert_allclose(ds.latitude.values[ok], t[ok] - 200, atol=2e-3)
    assert "native pixel positions" in ds.attrs["geolocation"]
    validate(ds, "swath")


def test_partial_granule_has_the_spacing_of_whole_ones(tmp_path):
    # granule 54 loses its first 60 scans (fill positions): its segment used
    # to get 135 rows over 75 scans' length (spacing 18 -> 10 km)
    _granule(tmp_path, 53)
    path = _granule(tmp_path, 54)

    def drop_first_60(lat, lon, value):
        value = value.copy()
        value[:60] = -9999
        return value
    _set(path, "l1_lat", drop_first_60)
    _set(path, "l1_lon", drop_first_60)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54)],
                                 file_path=tmp_path)
    native = load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54)],
                                 file_path=tmp_path, original_grid=True)
    for s in np.unique(ds.segment.values):
        step = np.median(np.diff(ds.distance_along.values[ds.segment.values == s]))
        rows = (native.segment.values == s)
        d = swath.swath_distances(native.latitude.values[rows].T,
                                  native.longitude.values[rows].T)["along"]
        np.testing.assert_allclose(step, np.median(np.diff(d)), rtol=0.02)
    assert ds.sizes["along_track"] == 135 + 75


def test_implausibly_narrow_scans_cannot_be_placed(tmp_path):
    # as AIRS 2015-01-07 granule 170: scans 1 % of the usual width, which set
    # the whole swath's cross-track axis (13.6 km wide instead of 1770)
    for g in (53, 54, 55):
        _granule(tmp_path, g)
    _granule(tmp_path, 56, cross_scale=0.01)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = load_nadir_sounders(None, "AIRS", granules=[(D, g) for g in
                                                          (56, 53, 54, 55)],
                                 file_path=tmp_path)
    width = ds.distance_cross.values[-1] - ds.distance_cross.values[0]
    assert 1600 < width < 1900
    assert 56 not in ds.granule.values  # (left out of the regridded swath)
    native = load_nadir_sounders(None, "AIRS", granules=[(D, 56), (D, 53)],
                                 file_path=tmp_path, original_grid=True)
    narrow = native.granule.values == 56
    assert np.all(states.has_state(native.temperature_flag.values[narrow],
                                   "bad_geolocation"))


def test_zero_weight_corner_does_not_remove_good_values():
    lat = np.array([[0.0, 0.2], [0.0, 0.2]])
    lon = np.array([[0.0, 0.0], [0.2, 0.2]])
    t = np.array([[220.0, 220.0], [220.0, np.nan]])
    out = swath.regrid(lat, lon, {"T": t}, n_cross=2, n_along=2)
    np.testing.assert_array_equal(out["fields"]["T"], t)


def _real_airs():
    try:
        folder = dataset_dir("AIRS")
    except (DataDirNotSetError, FileNotFoundError) as error:
        pytest.skip(f"no AIRS data: {error}")
    if not os.path.isfile(os.path.join(folder, "2007", "airs_2007d001.nc")):
        pytest.skip("AIRS 2007-01-01 not found")


def test_real_box_segments_are_unbroken():
    _real_airs()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = load_nadir_sounders(["2007-01-01T00:00", "2007-01-01T06:00"], "AIRS",
                                 lat_range=(-60, -40))
    step = ds.attrs["grid_spacing_km"][0]
    for s in np.unique(ds.segment.values):
        d = np.diff(ds.distance_along.values[ds.segment.values == s])
        assert d.size == 0 or d.max() < 1.5 * step


def test_real_granule_hardly_depends_on_what_is_loaded_with_it():
    # granules 125-126 alone and after 10-11 (0.5.1: 2.2 % of values lost and
    # others changed by up to 2.1 K, from the first segment's axis)
    _real_airs()
    alone = load_nadir_sounders(None, "AIRS", granules=[(D, 125), (D, 126)])
    both = load_nadir_sounders(None, "AIRS", granules=[(D, 10), (D, 11), (D, 125),
                                                       (D, 126)])
    a = alone.temperature.values
    b = both.temperature.values[np.isin(both.granule.values, [125, 126])]
    assert a.shape == b.shape
    lost = np.isfinite(a) & np.isnan(b)
    assert lost.mean() < 0.005
    finite = np.isfinite(a) & np.isfinite(b)
    assert np.percentile(np.abs(a[finite] - b[finite]), 99) < 0.5


def test_extrapolate_marks_only_points_outside_the_segment(tmp_path):
    # (review: with extrapolate=True whole edge columns were marked, from the
    # first scan's distances instead of the segment's own)
    _granule(tmp_path, 53)
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             extrapolate=True)
    edge = states.has_state(ds.temperature_flag.values[:, [0, -1]], "extrapolated")
    assert edge.mean() < 0.05
