"""
Known-answer tests of the CrIS reader (``obs._nadir_cris``) and of the
footprint regridding (``swath.footprint_coordinates``,
``swath.regrid_footprints``), on synthetic files laid out as L. Hoffmann's
daily files: a track along a meridian, 30 fields of regard (FOR) across, and
each FOR's 3 x 3 fields of view (FOV) rotated and spread towards the swath
edges, as in the real files.

Tolerances: positions are built with ``geodesy.reckon`` and measured back
with ``swath.footprint_coordinates`` on the same sphere, so distances agree
to 1e-6 km; a brightness temperature linear in latitude and longitude is
reproduced at the regridded values' own positions to 1e-3 K (barycentric
interpolation in distance is exact for a field linear in distance; the
remainder is the curvature of latitude and longitude against distance over
one triangle, about 20 km).
"""

import netCDF4
import numpy as np
import pytest

from atmosgw.curate import swath
from atmosgw.curate.data import availability, states, validate
from atmosgw.curate.geodesy import km2deg, reckon
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_nadir_sounders
from atmosgw.curate.obs._nadir_cris import index_layout

LAT0, LON0 = -20.0, 10.0  # the track: along the meridian 10 E, heading north
J2000_DAY = datenum(2000, 1, 1)


def field(lat, lon):
    """Brightness temperature linear in position (K)."""
    return 250.0 + 0.5 * (lat - LAT0) + 0.3 * (lon - LON0)


def geometry(n_scan):
    """lat, lon [scan x FOR x FOV] of the synthetic swath."""
    lat = np.empty((n_scan, 30, 9))
    lon = np.empty((n_scan, 30, 9))
    for j in range(30):
        u = (j - 14.5) / 14.5
        theta = np.deg2rad(50.0 * u)  # (rotation of the 3 x 3 pattern)
        spread = 1.0 + 1.5 * u ** 2
        for k in range(9):
            row, col = divmod(k, 3)
            a, c = (row - 1) * 16.0, (col - 1) * 16.0
            a, c = (spread * (a * np.cos(theta) - c * np.sin(theta)),
                    spread * (a * np.sin(theta) + c * np.cos(theta)))
            along = 53.0 * np.arange(n_scan) + a
            cross = 1100.0 * u + c
            t_lat = LAT0 + km2deg(along)
            lat[:, j, k], lon[:, j, k] = reckon(t_lat, np.full(n_scan, LON0),
                                                km2deg(np.full(n_scan, abs(cross))),
                                                np.full(n_scan, 90.0 if cross >= 0
                                                        else 270.0))
    return lat, lon


def write_cris(folder, day="2019-01-01", sat="snpp", n_scan=90, nan_fors=(0, 29),
               nan_footprint=None):
    """A daily file: n_scan scans from 00:00:12, every 8 s."""
    folder.mkdir(parents=True, exist_ok=True)
    d = datenum(int(day[:4]), int(day[5:7]), int(day[8:10]))
    doy = int(d - datenum(int(day[:4]), 1, 1)) + 1
    lat, lon = geometry(n_scan)
    seconds = ((d - J2000_DAY) * 86400 + 12.0 + 8.0 * np.arange(n_scan)[:, None, None]
               + 0.2 * np.arange(30)[None, :, None] + np.zeros((1, 1, 9)))
    bt = field(lat, lon)
    bt[:, list(nan_fors), :] = np.nan
    if nan_footprint is not None:
        bt[nan_footprint] = np.nan
    path = folder / f"pert_cris_{sat}_{day[:4]}_{doy:03d}.nc"
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("NTRACK", n_scan)
        f.createDimension("NXTRACK", 30)
        f.createDimension("NFOV", 9)
        dims = ("NTRACK", "NXTRACK", "NFOV")
        f.createVariable("time", "f8", dims)[:] = seconds
        f.createVariable("lat", "f8", dims)[:] = lat
        f.createVariable("lon", "f8", dims)[:] = lon
        for group, offset in (("4mu", 0.0), ("8mu", 1.0), ("10mu", 2.0),
                              ("15mu_low", 3.0), ("15mu_high", 4.0)):
            f.createVariable(f"bt_{group}", "f4", dims)[:] = bt + offset
            if group in ("4mu", "15mu_low", "15mu_high"):
                f.createVariable(f"bt_{group}_pt", "f4", dims)[:] = bt - 250.0
                f.createVariable(f"bt_{group}_var", "f4", dims)[:] = 1.0
    return path, lat, lon, bt


def test_footprint_coordinates_known_answers():
    # a meridian track, and points straight east and west of its middle rows
    lat_t = LAT0 + km2deg(20.0 * np.arange(5))
    lat = np.empty((4, 5))
    lon = np.empty((4, 5))
    for k, (dist, az) in enumerate(((500.0, 270.0), (1.0, 270.0), (1.0, 90.0),
                                    (500.0, 90.0))):
        lat[k], lon[k] = reckon(lat_t, np.full(5, LON0), km2deg(np.full(5, dist)),
                                np.full(5, az))
    along, cross = swath.footprint_coordinates(lat, lon)
    np.testing.assert_allclose(cross[:, 2], [-500, -1, 1, 500], atol=1e-6)
    np.testing.assert_allclose(along[1:3], np.broadcast_to(20.0 * np.arange(5),
                                                           (2, 5)), atol=1e-6)
    # far across track, the along-track coordinate is that of the row
    np.testing.assert_allclose(along[[0, 3], 2], 40.0, atol=1e-6)


def test_native_grid_is_the_index_layout(tmp_path):
    _, lat, lon, bt = write_cris(tmp_path)
    ds = load_nadir_sounders(None, "CrIS", product="2d", granules=[("2019-01-01", 1)],
                             file_path=tmp_path, original_grid=True)
    validate(ds, "swath")
    assert ds.sizes["along_track"] == 135 and ds.sizes["cross_track"] == 90
    # FOV k (0-based) of FOR j in scan s: cross 3 j + k % 3, along 3 s + k // 3
    for s, j, k in ((0, 0, 0), (3, 14, 4), (10, 29, 8), (44, 7, 5)):
        x, a = 3 * j + k % 3, 3 * s + k // 3
        assert ds.latitude.values[a, x] == lat[s, j, k]
        assert ds.longitude.values[a, x] == lon[s, j, k]
        assert np.array_equal(ds.brightness_temperature.values[a, x, 0],
                              np.float32(bt[s, j, k]), equal_nan=True)
    np.testing.assert_array_equal(index_layout(lat[:45]), ds.latitude.values.T)
    assert ds.attrs["satellite"] == "SNPP" and ds.attrs["grid"] == "native"
    # FORs 1 and 30 are NaN in the files: missing on the native grid
    edge = ds.brightness_temperature_flag.values[:, :3, 0]
    assert np.all(states.has_state(edge, "missing"))


def test_regridded_reproduces_a_linear_field(tmp_path):
    write_cris(tmp_path)
    ds = load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                             granules=[("2019-01-01", 1), ("2019-01-01", 2)],
                             channels=["4mu"])
    validate(ds, "swath")
    assert ds.attrs["grid"] == "regular"
    bt = ds.brightness_temperature.values[..., 0]
    flags = ds.brightness_temperature_flag.values[..., 0]
    good = np.isfinite(bt) & ~states.has_state(flags, "extrapolated")
    assert good.mean() > 0.95
    expected = field(ds.latitude.values, ds.longitude.values)
    assert np.max(np.abs(bt[good] - expected[good])) < 1e-3
    # the grid spans the FORs with values (2-29), at a regular spacing
    width = ds.distance_cross.values[-1] - ds.distance_cross.values[0]
    assert 1700 < width < 2200
    assert np.allclose(np.diff(ds.distance_cross.values),
                       ds.attrs["grid_spacing_km"][1], rtol=1e-9)
    # every NaN has a state; perturbation regridded the same way
    assert not np.any(np.isnan(bt) & ~states.no_value_mask(flags))
    pt = ds.brightness_temperature_perturbation.values[..., 0]
    np.testing.assert_allclose(pt[good], bt[good] - 250.0, atol=1e-3)
    for name in ("temperature", "pressure", "eastward_wind"):
        assert availability(ds, name) == "not_measured"


def test_missing_footprint_is_not_bridged(tmp_path):
    write_cris(tmp_path / "a")
    write_cris(tmp_path / "b", nan_footprint=(20, 15, 4))
    full, holed = (load_nadir_sounders(None, "CrIS", product="2d", channels=["4mu"],
                                       granules=[("2019-01-01", 1)],
                                       file_path=tmp_path / name)
                   for name in ("a", "b"))
    a = full.brightness_temperature.values[..., 0]
    b = holed.brightness_temperature.values[..., 0]
    lost = np.isfinite(a) & np.isnan(b)
    assert 0 < lost.sum() <= 6  # (only the grid points around the footprint)
    flags = holed.brightness_temperature_flag.values[..., 0]
    assert np.all(states.has_state(flags[lost], "missing"))
    same = np.isfinite(a) & np.isfinite(b)
    np.testing.assert_array_equal(a[same], b[same])


def test_nearest_takes_a_footprint_value(tmp_path):
    _, _, _, bt = write_cris(tmp_path)
    ds = load_nadir_sounders(None, "CrIS", product="2d", channels=["4mu"],
                             granules=[("2019-01-01", 1)], file_path=tmp_path,
                             interp_type="nearest")
    values = ds.brightness_temperature.values[..., 0]
    finite = values[np.isfinite(values)]
    assert np.all(np.isin(finite, bt.astype(np.float32)[np.isfinite(bt)]))


def test_granules_satellites_and_records(tmp_path):
    write_cris(tmp_path, sat="snpp")
    write_cris(tmp_path, sat="j1")
    # the first six minutes: granule 1 (scans 1-45, mean 00:03:08)
    ds = load_nadir_sounders(["2019-01-01T00:00", "2019-01-01T00:06"], "CrIS",
                             product="2d", file_path=tmp_path)
    assert set(ds.granule.values.tolist()) == {1}
    assert ds.attrs["files_not_found"] == ""  # (the day before is not reported)
    j1 = load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                             granules=[("2019-01-01", 2)], satellite="NOAA-20",
                             file_source=True)
    assert j1.attrs["satellite"] == "NOAA-20"
    assert list(j1.original_files.values) == ["pert_cris_j1_2019_001.nc"]
    assert set(j1.granule.values.tolist()) == {2}
    with pytest.warns(UserWarning, match="were found|not found"):
        ds = load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                                 granules=[("2019-01-01", 3), ("2019-01-02", 1)])
    assert "2019-01-01 granule 3" in ds.attrs["files_not_found"]
    assert "2019-01-02 granule 1" in ds.attrs["files_not_found"]
    with pytest.raises(ValueError, match="satellite"):
        load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                            granules=[("2019-01-01", 1)], satellite="Aqua")
    with pytest.raises(TypeError, match="CrIS only"):
        load_nadir_sounders(None, "AIRS", file_path=tmp_path,
                            granules=[("2019-01-01", 1)], satellite="SNPP")


def test_channels_and_perturbations(tmp_path):
    write_cris(tmp_path)
    ds = load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                             granules=[("2019-01-01", 1)], original_grid=True,
                             channels=["10mu", "15mu_high"])
    assert list(ds.channel_name.values) == ["10mu", "15mu_high"]
    bt = ds.brightness_temperature.values
    good = np.isfinite(bt[..., 0])
    np.testing.assert_allclose(bt[..., 1][good] - bt[..., 0][good], 2.0, atol=1e-4)
    pt = ds.brightness_temperature_perturbation.values
    assert np.all(np.isnan(pt[..., 0]))  # (no perturbation for 10.4 micron)
    assert np.any(np.isfinite(pt[..., 1]))
    with pytest.raises(ValueError, match="Unknown CrIS channel"):
        load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                            granules=[("2019-01-01", 1)], channels=["9mu"])
