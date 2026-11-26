"""Real-data check of load_time_height("ALOMAR") against the MATLAB-parity
path, load_limb_sounders(..., "ALOMAR", original_z=True) (itself checked
against MATLAB, tests/test_limb_matlab.py). Skipped without the ALOMAR files
(a skip is not a pass).

Both read the same files; rows are matched by (file, source_profile) and
levels by their 0.15 km lattice index. Tolerances:

- temperature, its NaN positions, its state flags and additional
  variables (T_ecmwf): identical (the same float64 values, only placed
  differently), with outliers removed (default) and kept;
- the same rows, none lost or added;
- heights: 1e-9 km (the time_height heights are the lattice values; the
  file values differ from them by at most about 1e-11 km);
- times: 1 ms (the limb loader rounds to the millisecond);
- pressure of to_profiles against the limb loader's estimated pressure:
  1e-11 relative (both from atmosphere.height_to_pressure, on heights
  differing by <= ~1e-11 km: dp/p = dz/H, H ~ 6-7 km, so ~2e-12).

Dates: 2007-01-23 (the MATLAB reference case); 2007-05-13 to 2007-05-14
(the height grid changes); 2007-07-02 (inside the session that started on
2007-06-30); 2007-06-04 (two sessions overlapping in time); 2005-12-30 to
2006-01-14 (two yearly files). A wider run (34 dates: a session of every
year 2003-2026, the first of each of the 6 height grids, an 8-day session)
was made once, 2026-11-26.
"""

import numpy as np
import pytest

from atmosgw.curate.data import to_profiles, validate
from atmosgw.curate.obs import load_limb_sounders, load_time_height
from atmosgw.curate.paths import DataDirNotSetError, dataset_dir

STEP = 0.15


def _data_or_skip():
    try:
        if dataset_dir("ALOMAR", must_exist=False).is_dir():
            return
    except DataDirNotSetError:
        pass
    pytest.skip("no ALOMAR data")


@pytest.mark.parametrize("keep_outliers", [True, False])
@pytest.mark.parametrize("time_range", ["2007-01-23", ["2007-05-13", "2007-05-14"],
                                        "2007-07-02", "2007-06-04",
                                        ["2005-12-30", "2006-01-14"]])
def test_same_values_as_limb_loader(time_range, keep_outliers):
    _data_or_skip()
    options = dict(keep_outliers=keep_outliers, file_source=True,
                   additional_vars=["T_ecmwf"])
    th = load_time_height(time_range, "ALOMAR", **options)
    validate(th, "time_height")
    with pytest.warns(UserWarning, match="load_time_height"):
        limb = load_limb_sounders(time_range, "ALOMAR", original_z=True, **options)
    assert th.sizes["time"] > 0

    # rows: the same (file, row) pairs
    def keys(ds, source_profile, source_file):
        files = ds.original_files.values
        return [(files[f - 1], int(p)) for f, p in zip(source_file, source_profile)]

    th_keys = keys(th, th.source_profile.values, th.source_file.values)
    limb_keys = keys(limb, limb.source_profile.values[:, 0].astype(int),
                     limb.source_file.values[:, 0].astype(int))
    assert sorted(th_keys) == sorted(limb_keys)
    assert len(set(th_keys)) == len(th_keys)
    row = {k: i for i, k in enumerate(th_keys)}
    r = np.array([row[k] for k in limb_keys])  # th row of each limb profile

    # levels: by lattice index
    th_k = np.round(th.altitude.values / STEP).astype(int)
    alt = limb.altitude.values
    has = np.isfinite(alt)
    k = np.round(np.where(has, alt, 0) / STEP).astype(int)
    col = np.searchsorted(th_k, k)
    assert np.all(th_k[np.minimum(col, th_k.size - 1)][has] == k[has])
    col = np.minimum(col, th_k.size - 1)
    rows = np.broadcast_to(r[:, None], alt.shape)

    np.testing.assert_allclose(th.altitude.values[col][has], alt[has], rtol=0,
                               atol=1e-9)
    for name in ("temperature", "temperature_flag", "T_ecmwf"):
        np.testing.assert_array_equal(th[name].values[rows, col][has],
                                      limb[name].values[has], err_msg=name)
    # every value of a row is at one of the limb row's heights
    placed = np.zeros(th.temperature.shape, dtype=bool)
    placed[rows[has], col[has]] = True
    assert np.all(np.isnan(th.temperature.values[~placed]))

    dt = np.abs(th.time.values[r][:, None] - limb.time.values).astype(np.int64)
    assert dt.max() <= 1000  # microseconds

    prof = to_profiles(th)
    validate(prof, "profile")
    np.testing.assert_allclose(prof.pressure.values[rows, col][has],
                               limb.pressure.values[has], rtol=1e-11, atol=0)
