"""2-D products in load_nadir_sounders, through a stand-in instrument.

No 2-D reader exists yet: these tests register a stand-in reader that returns
synthetic brightness temperatures [cross x along x channel], to check that the
core handles 2-D data (no levels; one or several channels) and that its output
follows the swath format. Exact checks only.
"""

import types

import numpy as np
import pytest

from atmosgw.curate.data import states, validate
from atmosgw.curate.geodesy import km2deg, reckon
from atmosgw.curate.obs import nadir
from atmosgw.curate.obs._limb_tools import record

DAY = 733043.0  # 2007-01-01


def make_reader(n_channels, with_list=True, fill=None):
    """A stand-in reader: one granule of 135 x 90 synthetic points."""
    def load(settings, inst_info, var_names):
        record(inst_info)
        n_al, n_x = 135, 90
        head = 350 + np.linspace(0, 4, n_al)
        lat, lon = np.empty(n_al), np.empty(n_al)
        lat[0], lon[0] = 10.0, 30.0
        for k in range(1, n_al):
            lat[k], lon[k] = reckon(lat[k - 1], lon[k - 1], km2deg(18), head[k - 1])
        off = np.linspace(-800, 800, n_x)
        la, lo = reckon(np.broadcast_to(lat, (n_x, n_al)),
                        np.broadcast_to(lon, (n_x, n_al)),
                        km2deg(np.broadcast_to(off[:, None], (n_x, n_al))),
                        np.broadcast_to(head + 90, (n_x, n_al)))
        bt = 230 + np.arange(n_channels)[None, None, :] + 0 * la[..., None]
        if fill is not None:
            bt[fill] = np.nan
        data = {"Lat": np.asarray(la), "Lon": np.asarray(lo),
                "Time": DAY + 0.2 + np.broadcast_to(np.arange(n_al) * 8 / 3 / 86400,
                                                    (n_x, n_al)),
                "Temp": bt, "Granule": np.full(n_al, 53.0), "Day": np.full(n_al, DAY),
                "SourceFile": np.ones(n_al), "SourceRow": np.arange(1.0, n_al + 1)}
        if with_list:
            data["Channel"] = 600.0 + 10 * np.arange(n_channels)
        return data, ["synthetic.nc"]
    return types.SimpleNamespace(load=load)


@pytest.fixture
def stand_in(monkeypatch):
    def register(reader):
        entry = {**nadir.INSTRUMENTS["AIRS"], "reader": reader, "products": ("2d",),
                 "characterisation": {"2d": "TEST-2D-placeholder/0"}}
        monkeypatch.setitem(nadir.INSTRUMENTS, "TEST2D", entry)
    return register


def load(**options):
    return nadir.load_nadir_sounders(None, "TEST2D", product="2d",
                                     granules=[("2007-01-01", 53)],
                                     file_path="/unused", **options)


def test_several_channels(stand_in):
    stand_in(make_reader(3))
    ds = load()
    validate(ds, "swath")
    assert ds.brightness_temperature.dims == ("along_track", "cross_track", "channel")
    assert "altitude" not in ds and "level" not in ds.dims
    np.testing.assert_array_equal(ds.channel_wavenumber.values, [600, 610, 620])
    for c in range(3):  # a constant field stays exactly constant through regridding
        assert np.all(ds.brightness_temperature.values[..., c] == np.float32(230 + c))
    assert ds.attrs["spatial_dims"] == "along_track cross_track"
    assert np.isnan(ds.attrs["grid_spacing_km"][2])
    assert ds.brightness_temperature.attrs["origin"] == "measured"


def test_one_channel_without_list_has_no_channel_dimension(stand_in):
    stand_in(make_reader(1, with_list=False))
    ds = load()
    validate(ds, "swath")
    assert ds.brightness_temperature.dims == ("along_track", "cross_track")
    assert ds.brightness_temperature_flag.dims == ("along_track", "cross_track")


def test_fill_values_and_native_grid(stand_in):
    stand_in(make_reader(2, fill=(slice(10, 12), slice(50, 52), 1)))
    ds = load(original_grid=True)
    validate(ds, "swath")
    flag = ds.brightness_temperature_flag.values
    assert states.has_state(flag, "missing").sum() == 4
    assert np.isnan(ds.brightness_temperature.values[50:52, 10:12, 1]).all()
    assert ds.attrs["grid"] == "native"


def test_3d_only_options_are_ignored_for_2d(stand_in):
    stand_in(make_reader(2))
    ds = load(height_scale=[10, 20], original_z=True)
    validate(ds, "swath")
    assert "level" not in ds.dims
