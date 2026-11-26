"""data.to_netcdf and data.load_netcdf: every format read back identical
(unpacked) or within half the stated resolution (packed); the integer type
chosen by range; invalid Datasets, existing files and failed writes leave no
file; a packed file re-saved unpacked; empty Datasets. Synthetic inputs only
(written by the other tests' helpers), so these run in CI."""

import warnings

import netCDF4
import numpy as np
import pytest
import xarray as xr
from test_formats_round_trip import _misc, small_swath, small_time_height
from test_nadir_known import D, seconds, write_granule

from atmosgw.curate.data import FormatError, load_netcdf, to_netcdf, to_profiles
from atmosgw.curate.data import netcdf as nc
from atmosgw.curate.obs import load_limb_sounders, load_nadir_sounders


def _profile(tmp_path):
    return load_limb_sounders("2007-01-01", "Misc", original_z=True, file_source=True,
                              dx=True, **_misc(tmp_path))


def _profile_interpolated(tmp_path):
    return load_limb_sounders("2007-01-01", "Misc", **_misc(tmp_path))


def _swath(tmp_path):
    return small_swath(tmp_path)


def _swath_options(tmp_path):
    write_granule(tmp_path, 53, seconds(53))
    return load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                               n_cross=12, n_along=9, file_source=True,
                               keep_native_geolocation=True, day_night=True)


def _time_height(tmp_path):
    return small_time_height(tmp_path)


def _swath_profiles(tmp_path):
    return to_profiles(small_swath(tmp_path))


BUILDERS = [_profile, _profile_interpolated, _swath, _swath_options, _time_height,
            _swath_profiles]


def _stored_types(path):
    """{variable: type in the file}."""
    with netCDF4.Dataset(path) as f:
        return {name: v.dtype for name, v in f.variables.items()}


def _same_attrs(a, b):
    assert set(a) == set(b)
    for key in a:
        if isinstance(a[key], (list, tuple, np.ndarray)):
            np.testing.assert_array_equal(np.asarray(a[key]), np.asarray(b[key]))
        else:
            assert a[key] == b[key], key


def _assert_identical(ds, back):
    assert set(ds.variables) == set(back.variables)
    assert dict(ds.sizes) == dict(back.sizes)
    _same_attrs(ds.attrs, back.attrs)
    for name in ds.variables:
        a, b = ds[name], back[name]
        assert a.dims == b.dims, name
        text = a.dtype.kind in "US" and b.dtype.kind in "US"  # (str length may differ)
        assert a.dtype == b.dtype or text, name
        np.testing.assert_array_equal(a.values, b.values, err_msg=name)
        _same_attrs(a.attrs, b.attrs)


@pytest.mark.parametrize("build", BUILDERS, ids=lambda f: f.__name__)
def test_every_format_reads_back_identical(tmp_path, build):
    ds = build(tmp_path)
    path = tmp_path / "out.nc"
    to_netcdf(ds, path)
    back = load_netcdf(path)
    _assert_identical(ds, back)
    stored = _stored_types(path)
    assert stored["time"] == np.int64  # (exact microseconds)
    assert all(t != np.int16 for t in stored.values() if not isinstance(t, type))


def _packed_close(a, b, resolution):
    """Within half the resolution, plus the rounding of the variable's type."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    np.testing.assert_array_equal(np.isnan(a), np.isnan(b))
    rounding = (np.finfo(np.float32).eps * np.nanmax(np.abs(a)) * 4
                if np.asarray(a).dtype == np.float32 else 1e-9 * resolution)
    assert np.nanmax(np.abs(a - b)) <= resolution / 2 + rounding


def test_packed_swath(tmp_path):
    ds = _swath_options(tmp_path)
    path = tmp_path / "packed.nc"
    to_netcdf(ds, path, pack=True)
    back = load_netcdf(path)
    stored = _stored_types(path)
    assert stored["temperature"] == np.int16
    assert stored["latitude"] == np.int32 and stored["longitude"] == np.int32
    assert stored["native_latitude"] == np.int32
    # grid axes, times, flags and indices as they are
    assert stored["altitude"] == np.float64 and stored["distance_along"] == np.float64
    assert stored["time"] == np.int64 and stored["temperature_flag"] == np.uint16
    assert back.temperature.dtype == np.float32
    assert back.latitude.dtype == np.float64
    # (0.01 K: values within 0.005 K, plus float32 rounding near 300 K)
    _packed_close(ds.temperature.values, back.temperature.values, 0.01)
    _packed_close(ds.latitude.values, back.latitude.values, 1e-7)
    _packed_close(ds.longitude.values, back.longitude.values, 1e-7)
    for name in ("time", "temperature_flag", "granule", "distance_along",
                 "distance_cross", "altitude", "day_retrieval"):
        np.testing.assert_array_equal(ds[name].values, back[name].values)


def test_packed_profiles_keep_float64_and_leave_pressure(tmp_path):
    ds = _profile(tmp_path)
    # (computed distances, not float32 values, so that they are packed)
    dx = np.where(np.isnan(ds.dx.values), np.nan, 123.456789012 + ds.altitude.values)
    ds["dx"] = (ds.dx.dims, dx, ds.dx.attrs)
    path = tmp_path / "packed.nc"
    to_netcdf(ds, path, pack=True)
    back = load_netcdf(path)
    stored = _stored_types(path)
    assert stored["temperature"] == np.int16
    assert stored["altitude"] == np.float64  # (the file's own values: float32-exact)
    assert stored["dx"] in (np.int16, np.int32)  # (computed; type by its range)
    assert stored["pressure"] == np.float64
    assert stored["source_profile"] == np.float64
    assert back.temperature.dtype == np.float64
    _packed_close(ds.temperature.values, back.temperature.values, 0.01)
    _packed_close(ds.dx.values, back.dx.values, 1e-6)
    np.testing.assert_array_equal(ds.altitude.values, back.altitude.values)
    np.testing.assert_array_equal(ds.pressure.values, back.pressure.values)


def test_packed_time_height_uncertainty_and_site(tmp_path):
    ds = _time_height(tmp_path)
    path = tmp_path / "packed.nc"
    to_netcdf(ds, path, pack=True)
    back = load_netcdf(path)
    stored = _stored_types(path)
    assert stored["temperature_uncertainty"] in (np.int16, np.int32)
    assert stored["altitude"] == np.float64  # (a grid axis)
    assert stored["latitude"] == np.float64  # (the site: no dimensions)
    _packed_close(ds.temperature_uncertainty.values,
                  back.temperature_uncertainty.values, 0.001)
    np.testing.assert_array_equal(ds.altitude.values, back.altitude.values)


@pytest.mark.parametrize("values, expected", [
    ([100.0, 400.0], "int16"),
    ([0.0, 655.33], "int16"),  # (65533 steps: the int16 range less the fill)
    ([0.0, 700.0], "int32"),
    ([0.0, 1e8], None),  # (too wide for int32 at 0.01 K: stored unpacked)
    ([np.nan, np.nan], None),
    ([250.0, np.inf], None),
])
def test_integer_type_chosen_by_range(values, expected):
    variable = xr.Variable(("a", "b"), np.array([values], dtype=np.float64))
    packing = nc._packing(variable, 0.01)
    assert packing.get("dtype") == expected
    if expected:  # (every value representable: none clipped)
        steps = np.round((np.asarray(values) - packing["add_offset"]) / 0.01)
        info = np.iinfo(expected)
        assert np.all(steps > info.min) and np.all(steps <= info.max)


def test_extreme_values_round_trip_without_clipping(tmp_path):
    ds = _profile(tmp_path)
    temp = ds.temperature.values.copy()
    temp[0, 0], temp[0, 1] = 100.0, 760.0  # (span 660 K: needs int32 at 0.01 K)
    ds["temperature"] = (ds.temperature.dims, temp, ds.temperature.attrs)
    path = tmp_path / "wide.nc"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        to_netcdf(ds, path, pack=True)
    assert _stored_types(path)["temperature"] == np.int32
    _packed_close(temp, load_netcdf(path).temperature.values, 0.01)


def _leftovers(folder):
    return sorted(p.name for p in folder.iterdir() if p.suffix in (".nc", ".partial"))


def test_invalid_dataset_is_refused_and_nothing_written(tmp_path):
    ds = _profile(tmp_path)
    bad = ds.copy()
    bad["temperature"].attrs["units"] = "degC"
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(FormatError, match="units"):
        to_netcdf(bad, out / "bad.nc")
    plain = xr.Dataset({"x": ("a", np.arange(3.0))})
    with pytest.raises(FormatError, match="format"):
        to_netcdf(plain, out / "plain.nc")
    assert _leftovers(out) == []


def test_existing_file_needs_clobber(tmp_path):
    ds = _profile(tmp_path)
    path = tmp_path / "out.nc"
    to_netcdf(ds, path)
    before = path.read_bytes()
    changed = ds.copy()
    changed.attrs["options"] = "changed"
    with pytest.raises(FileExistsError, match="clobber"):
        to_netcdf(changed, path)
    assert path.read_bytes() == before
    to_netcdf(changed, path, clobber=True)
    assert load_netcdf(path).attrs["options"] == "changed"
    assert _leftovers(tmp_path) == ["out.nc"]


def test_failed_write_leaves_no_file(tmp_path, monkeypatch):
    ds = _profile(tmp_path)
    original = xr.Dataset.to_netcdf

    def write_then_fail(self, path, *args, **kwargs):  # (a file exists, then error)
        original(self, path, *args, **kwargs)
        raise RuntimeError("disk full")

    monkeypatch.setattr(xr.Dataset, "to_netcdf", write_then_fail)
    with pytest.raises(RuntimeError, match="disk full"):
        to_netcdf(ds, tmp_path / "new.nc")
    assert _leftovers(tmp_path) == []
    # an existing file is left as it was
    monkeypatch.setattr(xr.Dataset, "to_netcdf", original)
    to_netcdf(ds, tmp_path / "old.nc")
    before = (tmp_path / "old.nc").read_bytes()
    monkeypatch.setattr(xr.Dataset, "to_netcdf", write_then_fail)
    with pytest.raises(RuntimeError):
        to_netcdf(ds, tmp_path / "old.nc", clobber=True)
    assert (tmp_path / "old.nc").read_bytes() == before
    assert _leftovers(tmp_path) == ["old.nc"]


def test_packed_file_resaved_unpacked(tmp_path):
    ds = _swath(tmp_path)
    to_netcdf(ds, tmp_path / "packed.nc", pack=True)
    back = load_netcdf(tmp_path / "packed.nc")
    to_netcdf(back, tmp_path / "again.nc")  # (the file's packing is not reused)
    stored = _stored_types(tmp_path / "again.nc")
    assert stored["temperature"] == np.float32 and stored["latitude"] == np.float64
    _assert_identical(back, load_netcdf(tmp_path / "again.nc"))


def test_empty_datasets(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        limb = load_limb_sounders("2007-01-05", "Misc", **_misc(tmp_path))
        swath = load_nadir_sounders(["2007-01-05T00:00", "2007-01-05T01:00"], "AIRS",
                                    file_path=tmp_path)
    assert limb.sizes["profile"] == 0 and swath.sizes["along_track"] == 0
    for k, ds in enumerate((limb, swath)):
        for pack in (False, True):
            path = tmp_path / f"empty_{k}_{pack}.nc"
            to_netcdf(ds, path, pack=pack)
            _assert_identical(ds, load_netcdf(path))


def test_load_netcdf_refuses_files_not_in_a_format(tmp_path):
    xr.Dataset({"x": ("a", np.arange(3.0))}).to_netcdf(tmp_path / "plain.nc")
    with pytest.raises(FormatError, match="xarray.load_dataset"):
        load_netcdf(tmp_path / "plain.nc")


def test_options_must_be_true_or_false(tmp_path):
    ds = _profile(tmp_path)
    with pytest.raises(TypeError):
        to_netcdf(ds, tmp_path / "x.nc", pack="yes")
    with pytest.raises(TypeError):
        to_netcdf(ds, tmp_path / "x.nc", clobber=1)
    with pytest.raises(FileNotFoundError, match="does not exist"):
        to_netcdf(ds, tmp_path / "no_such_folder" / "x.nc")


def test_float32_exact_positions_are_not_packed(tmp_path):
    ds = _profile(tmp_path)  # (positions and heights read from the file: exact)
    lat = ds.latitude.values
    assert np.array_equal(lat.astype(np.float32).astype(float), lat, equal_nan=True)
    path = tmp_path / "packed.nc"
    to_netcdf(ds, path, pack=True)
    stored = _stored_types(path)
    assert stored["latitude"] == np.float64 and stored["altitude"] == np.float64
    np.testing.assert_array_equal(load_netcdf(path).latitude.values, lat)
    computed = ds.copy()
    computed["latitude"] = (ds.latitude.dims, lat + 1e-9, ds.latitude.attrs)
    to_netcdf(computed, tmp_path / "computed.nc", pack=True)
    assert _stored_types(tmp_path / "computed.nc")["latitude"] in (np.int16, np.int32)


def test_missing_text_reads_back_empty(tmp_path):
    ds = _profile(tmp_path)
    ds["note"] = ("profile", np.array(["a", np.nan], dtype=object))
    to_netcdf(ds, tmp_path / "text.nc")
    assert load_netcdf(tmp_path / "text.nc").note.values.tolist() == ["a", ""]
