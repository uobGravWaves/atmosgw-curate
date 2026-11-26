"""Failures are loud: set-up faults raise, missing or empty data warn, bad
inputs are refused with a message that says how to recover (obs/_common.py,
the loaders' _prepare, the readers).

Synthetic files only (Misc MAT-files, MLS HDF5, netCDF), so every test runs
in CI.
"""

import os
import warnings

import h5py
import netCDF4
import numpy as np
import pytest
from scipy.io import savemat

from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_limb_sounders, load_nadir_sounders, load_time_height
from atmosgw.curate.readers import read_netcdf


def _misc(folder, day, temp=230.0, n_lev=16, **override):
    alt = np.tile(np.arange(10.0, 10.0 + 2 * n_lev, 2.0), (2, 1))
    content = {"Alt": alt, "Lat": np.full_like(alt, 30.0),
               "Lon": np.full_like(alt, 10.0), "Temp": np.full_like(alt, temp),
               "Time": np.full_like(alt, datenum(2007, 1, day) + 0.5)}
    content.update(override)
    savemat(os.path.join(folder, f"mine_2007d{day:03d}.mat"), content)


def _load(folder, time_range="2007-01-01", **options):
    return load_limb_sounders(time_range, "Misc", misc_path=folder,
                              misc_identifier="mine", **options)


@pytest.mark.parametrize("loader, args", [
    (load_limb_sounders, ("2007-01-01", "SABER")),
    (load_nadir_sounders, ("2007-01-01T05:00", "AIRS")),
    (load_time_height, ("2007-01-23", "ALOMAR")),
])
def test_missing_data_folder_raises(tmp_path, monkeypatch, loader, args):
    monkeypatch.setenv("ATMOSGW_DATA", str(tmp_path))
    monkeypatch.setenv("ATMOSGW_CONFIG", str(tmp_path / "none.toml"))
    with pytest.raises(FileNotFoundError, match="does not exist"):
        loader(*args)
    with pytest.raises(FileNotFoundError, match="from file_path"):
        loader(*args, file_path=tmp_path / "typo")


def test_missing_days_warned(tmp_path):
    _misc(tmp_path, 1)
    _misc(tmp_path, 3)
    with pytest.warns(UserWarning, match="1 files looked for were not found"):
        ds = _load(tmp_path, ["2007-01-01", "2007-01-03"])
    assert ds.attrs["files_not_found"] == "2007-01-02"
    with pytest.warns(UserWarning, match="None of the Misc files"):
        _load(tmp_path, "2007-01-05")
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # (date_warning=False silences it)
        _load(tmp_path, ["2007-01-01", "2007-01-03"], date_warning=False)


def test_file_without_data_warned(tmp_path):
    _misc(tmp_path, 1, temp=np.nan)
    with pytest.warns(UserWarning, match="hold no values"):
        ds = _load(tmp_path)
    assert ds.attrs["files_without_data"] == "mine_2007d001.mat"


def test_variables_of_the_wrong_shape_refused(tmp_path):
    _misc(tmp_path, 1, Temp=np.full((2, 13), 230.0))
    with pytest.raises(ValueError, match="mine_2007d001.mat: Temp has shape"):
        _load(tmp_path)


def test_unreadable_file_named(tmp_path):
    (tmp_path / "mine_2007d001.mat").write_bytes(b"not a MAT-file")
    with pytest.raises(OSError, match="mine_2007d001.mat"):
        _load(tmp_path)


def test_missing_variable_named(tmp_path):
    _misc(tmp_path, 1)
    with pytest.raises(KeyError, match="Extra not found in .*mine_2007d001.mat"):
        _load(tmp_path, additional_vars=["Extra"])


def test_height_range_outside_the_height_scale(tmp_path):
    _misc(tmp_path, 1)
    with pytest.raises(ValueError, match="leaves no heights of the height scale"):
        _load(tmp_path, height_range=(100, 120))


def test_one_or_two_levels_give_nan_not_an_error(tmp_path):
    for n_lev in (1, 2):
        _misc(tmp_path, 1, n_lev=n_lev)
        ds = _load(tmp_path)
        assert ds.sizes["profile"] == 2 and np.all(np.isnan(ds.temperature.values))


@pytest.mark.parametrize("value, message", [
    ("20070101", "write it as '2007-01-01'"),
    (None, "empty"),
    ([], "empty"),
])
def test_bad_time_inputs(tmp_path, value, message):
    _misc(tmp_path, 1)
    with pytest.raises(ValueError, match=message):
        _load(tmp_path, value)


def test_datetime64_arrays_of_any_unit(tmp_path):
    _misc(tmp_path, 1)
    times = np.array(["2007-01-01T00:00", "2007-01-01T23:00"], dtype="datetime64[ns]")
    assert _load(tmp_path, times).sizes["profile"] == 2


@pytest.mark.parametrize("options, message", [
    ({"granules": [("2007-01-01", 0)]}, "1 to 240"),
    ({"granules": [("2007-01-01", 241)]}, "1 to 240"),
    ({"granules": [("2007-01-01", 5.5)]}, "1 to 240"),
    ({"granules": [("2007-01-01", 5)], "n_cross": 0}, "n_cross"),
    ({"granules": [("2007-01-01", 5)], "dx_along": -1}, "dx_along"),
    ({"granules": [("2007-01-01", 5)], "lat_range": (-100, 100)}, "lat_range"),
])
def test_nadir_options_checked(tmp_path, options, message):
    with pytest.raises(ValueError, match=message):
        load_nadir_sounders(None, "AIRS", file_path=tmp_path, **options)


def test_repeated_granule_loaded_once(tmp_path):
    with pytest.warns(UserWarning, match="None of the AIRS"):
        ds = load_nadir_sounders(None, "AIRS", file_path=tmp_path,
                                 granules=[("2007-01-01", 5), ("2007-01-01", 5)])
    assert ds.attrs["files_not_found"] == "2007-01-01 granule 5"


def test_mls_file_without_temperature_refused(tmp_path):
    (tmp_path / "2026").mkdir()
    path = tmp_path / "2026" / "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d001.he5"
    with h5py.File(path, "w") as f:
        geo = "HDFEOS/SWATHS/Temperature/Geolocation Fields"
        for name in ("Latitude", "Longitude", "Time"):
            f[f"{geo}/{name}"] = np.zeros(3)
        f[f"{geo}/Pressure"] = np.array([100.0, 10, 1])
    with pytest.raises(KeyError, match="L2gpValue not found in .*2026d001.he5"):
        load_limb_sounders("2026-01-01", "MLS", file_path=tmp_path)


def test_every_missing_value_masked(tmp_path):
    path = tmp_path / "x.nc"
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("n", 3)
        v = f.createVariable("t", "f8", ("n",))
        v.missing_value = np.array([-999.0, 200.0])
        v[:] = [220.0, -999.0, 200.0]
    np.testing.assert_array_equal(read_netcdf(path, ["t"])["t"],
                                  [220.0, np.nan, np.nan])


def test_one_row_of_heights_is_per_level(tmp_path):
    # (a Misc file may give its heights as one row, as GNSS gives MSL_alt)
    alt = np.arange(10.0, 41.0, 2.0)[None, :]
    savemat(tmp_path / "mine_2007d001.mat", {
        "Alt": alt, "Lat": np.full((2, 16), 30.0), "Lon": np.full((2, 16), 10.0),
        "Temp": np.full((2, 16), 230.0),
        "Time": np.full((2, 16), datenum(2007, 1, 1) + 0.5)})
    ds = _load(tmp_path, original_z=True)
    np.testing.assert_array_equal(ds.altitude.values, np.tile(alt, (2, 1)))


def test_padding_and_neighbours_not_reported(tmp_path):
    # IASI's default screen reads beyond the request: the days before and
    # after are not reported missing (nothing is in the folder, so only the
    # day asked for is)
    with pytest.warns(UserWarning, match="2010-07-01"):
        ds = load_nadir_sounders(["2010-07-01T00:00", "2010-07-01T00:20"], "IASI",
                                 file_path=tmp_path)
    assert ds.attrs["files_not_found"] == "2010-07-01"
