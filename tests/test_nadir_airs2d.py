"""Tests for the AIRS 2-D reader (L. Hoffmann's granule files), on synthetic
files written to a temporary folder in the files' layout. Exact checks."""

import netCDF4
import numpy as np
import pytest

from atmosgw.curate.data import states, validate
from atmosgw.curate.geodesy import km2deg, reckon
from atmosgw.curate.obs import load_nadir_sounders

DAY = "2020-01-31"
GROUPS = ("4mu", "8mu", "15mu_low", "15mu_high")


def write(root, granule, year_folder=False, bt_value=None):
    """One granule file: 135 x 90, a curved track, constant brightness temperatures."""
    folder = root / ("2020" if year_folder else "") / "031"
    folder.mkdir(parents=True, exist_ok=True)
    n_al, n_x = 135, 90
    head = 10 + np.linspace(0, 4, n_al)
    lat, lon = np.empty(n_al), np.empty(n_al)
    lat[0], lon[0] = -40.0 + granule * 0.1, 20.0
    for k in range(1, n_al):
        lat[k], lon[k] = reckon(lat[k - 1], lon[k - 1], km2deg(18), head[k - 1])
    la, lo = reckon(np.broadcast_to(lat[:, None], (n_al, n_x)),
                    np.broadcast_to(lon[:, None], (n_al, n_x)),
                    km2deg(np.broadcast_to(np.linspace(-800, 800, n_x), (n_al, n_x))),
                    np.broadcast_to(head[:, None] + 90, (n_al, n_x)))
    start = ((np.datetime64(DAY) - np.datetime64("2000-01-01")) / np.timedelta64(1, "s")
             + 328 + (granule - 1) * 360)
    path = folder / f"airs_2020_031_{granule:03d}.nc"
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("NTRACK", None)
        f.createDimension("NXTRACK", n_x)
        dims = ("NTRACK", "NXTRACK")
        f.createVariable("time", "f8", dims)[:] = start + np.broadcast_to(
            np.arange(n_al)[:, None] * 8 / 3, (n_al, n_x))
        f.createVariable("lat", "f8", dims)[:] = la
        f.createVariable("lon", "f8", dims)[:] = lo
        for k, g in enumerate(GROUPS):
            bt = np.full((n_al, n_x), 230.0 + k, dtype=np.float32)
            if bt_value is not None and g == "4mu":
                bt[60, 45] = bt_value
            f.createVariable(f"bt_{g}", "f4", dims)[:] = bt
            if g != "8mu":
                f.createVariable(f"bt_{g}_pt", "f4", dims)[:] = np.float32(0.5 + k)
                f.createVariable(f"bt_{g}_var", "f4", dims)[:] = np.float32(2.0)
    return path


def test_channels_derived_variables_and_format(tmp_path):
    write(tmp_path, 121)
    ds = load_nadir_sounders(None, "AIRS", product="2d", file_path=tmp_path,
                             granules=[(DAY, 121)])
    validate(ds, "swath")
    assert list(ds.channel_name.values) == list(GROUPS)
    np.testing.assert_allclose(ds.channel_wavenumber.values,
                               [1e4 / 4.3, 1e4 / 8.1, 1e4 / 15, 1e4 / 15])
    for k in range(4):  # constant fields stay exactly constant through regridding
        assert np.all(ds.brightness_temperature.values[..., k] == np.float32(230 + k))
    pt = ds.brightness_temperature_perturbation.values
    assert np.all(pt[..., 0] == np.float32(0.5)) and np.isnan(pt[..., 1]).all()
    flag = ds.brightness_temperature_perturbation_flag.values
    assert states.has_state(flag[..., [0, 2, 3]], "derived").all()
    assert states.has_state(flag[..., 1], "no_observation").all()
    assert ds.brightness_temperature_perturbation.attrs["origin"] == "derived"
    assert ds.brightness_temperature_variance.attrs["units"] == "K2"


def test_channel_choice_and_year_folders(tmp_path):
    write(tmp_path, 121, year_folder=True)
    ds = load_nadir_sounders(None, "AIRS", product="2d", file_path=tmp_path,
                             granules=[(DAY, 121)], channels=["15mu_high", "4mu"])
    validate(ds, "swath")
    assert list(ds.channel_name.values) == ["15mu_high", "4mu"]
    assert np.all(ds.brightness_temperature.values[..., 0] == np.float32(233))
    with pytest.raises(ValueError):
        load_nadir_sounders(None, "AIRS", product="2d", file_path=tmp_path,
                            granules=[(DAY, 121)], channels=["9mu"])


def test_bad_value_removed_everywhere_it_matters(tmp_path):
    write(tmp_path, 121, bt_value=500.0)  # outside the physical range
    ds = load_nadir_sounders(None, "AIRS", product="2d", file_path=tmp_path,
                             granules=[(DAY, 121)], original_grid=True)
    validate(ds, "swath")
    flag = ds.brightness_temperature_flag.values
    assert states.has_state(flag[60, 45, 0], "outside_physical_range")
    assert np.isnan(ds.brightness_temperature.values[60, 45, 0])
    assert np.isnan(ds.brightness_temperature_perturbation.values[60, 45, 0])


def test_missing_granule_recorded(tmp_path):
    write(tmp_path, 121)
    ds = load_nadir_sounders(None, "AIRS", product="2d", file_path=tmp_path,
                             granules=[(DAY, 121), (DAY, 122)])
    assert "granule 122" in ds.attrs["files_not_found"]
    assert set(ds.granule.values) == {121}
