"""Formats: saved files read back as valid Datasets, extras' own dimensions
are accepted, empty planetary-wave results keep their mode dimension, and
to_profiles keeps a supplied pressure and uncertainty and does not share
memory with its input."""

import numpy as np
import xarray as xr
from scipy.io import savemat
from test_nadir_known import D, seconds, write_granule
from test_time_height_known import sessions_2007, write_year

from atmosgw.curate.data import problems, to_profiles, validate
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_limb_sounders, load_nadir_sounders, load_time_height


def small_swath(tmp_path):
    write_granule(tmp_path, 53, seconds(53))
    return load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                               n_cross=12, n_along=9)


def small_time_height(tmp_path):
    write_year(tmp_path, 2007, sessions_2007())
    return load_time_height(["2007-01-24", "2007-01-25"], "ALOMAR",
                            file_path=tmp_path)


def _misc(tmp_path, **extra):
    alt = np.tile(np.arange(10.0, 41.0, 2.0), (2, 1))
    savemat(tmp_path / "mine_2007d001.mat", {
        "Alt": alt, "Lat": np.full_like(alt, 30.0), "Lon": np.full_like(alt, 10.0),
        "Temp": np.full_like(alt, 230.0),
        "Time": np.full_like(alt, datenum(2007, 1, 1) + 0.5), **extra})
    return dict(misc_path=tmp_path, misc_identifier="mine")


def test_saved_file_reads_back_valid(tmp_path):
    ds = load_limb_sounders("2007-01-01", "Misc", original_z=True, **_misc(tmp_path))
    ds.to_netcdf(tmp_path / "out.nc")
    back = xr.load_dataset(tmp_path / "out.nc")
    assert back.time.dtype == np.dtype("datetime64[ns]")  # (xarray's default)
    validate(back, "profile")
    np.testing.assert_array_equal(back.time.values.astype("datetime64[us]"),
                                  ds.time.values)


def test_extras_own_dimensions_are_valid(tmp_path):
    opts = _misc(tmp_path, Meta3=np.zeros((2, 16, 4)))
    ds = load_limb_sounders("2007-01-01", "Misc", additional_vars=["Meta3"], **opts)
    assert "Meta3_dim2" in ds.dims
    assert problems(ds, "profile") == []


def test_to_profiles_keeps_supplied_pressure_and_uncertainty(tmp_path):
    swath = small_swath(tmp_path)
    shape = swath.temperature.shape
    swath["pressure"] = (swath.temperature.dims, np.full(shape, 123.0, np.float32),
                         {"units": "hPa", "origin": "retrieved",
                          "ancillary_variables": "pressure_flag"})
    swath["pressure_flag"] = (swath.temperature.dims,
                              np.zeros(shape, swath.temperature_flag.dtype),
                              swath.temperature_flag.attrs)
    swath["temperature_uncertainty"] = (swath.temperature.dims,
                                        np.full(shape, 0.5, np.float32),
                                        {"units": "K", "comment": "1-sigma"})
    prof = to_profiles(swath)
    assert np.all(prof.pressure.values == 123.0)
    assert prof.pressure.attrs["origin"] == "retrieved"
    assert np.all(prof.temperature_uncertainty.values == 0.5)


def test_to_profiles_does_not_share_memory_with_its_input(tmp_path):
    th = small_time_height(tmp_path)
    prof = to_profiles(th)
    prof.temperature.values[...] = -1
    assert not np.any(th.temperature.values == -1)
