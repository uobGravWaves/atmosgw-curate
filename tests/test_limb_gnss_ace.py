"""ACE and GNSS: pressure units and the GNSS corrections (Q9-Q11 in
docs/quirks.md; obs/_limb_ace_gnss.py).

Synthetic: an ACE file written with scipy (pressure in atm); GNSS files as
read_mat would return them (the MATLAB table MetaData cannot be written from
Python), given to the reader by replacing read_mat. Known answers: the
conversion factors, MetaData positions at every level, the occultation id's
time for hour-rounded times, and which copy of a repeated occultation is kept.

Real data (skipped if absent; a skip is not a pass): pressure within a factor
of 3 of the U.S. Standard Atmosphere at 12-40 km for every limb instrument,
which fails for a unit error of 100 (Pa) or 1013 (atm).
"""

import os
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy.io import savemat

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.data import validate
from atmosgw.curate.matlab_time import datenum, datenum_to_datetime
from atmosgw.curate.obs import _limb_ace_gnss, load_limb_sounders
from atmosgw.curate.paths import DataDirNotSetError, dataset_dir

D1 = "2007-01-01"
T0 = datenum(2007, 1, 1)


def _ace(tmp_path):
    folder = tmp_path / "2007"
    folder.mkdir()
    alt = np.tile(np.arange(10.0, 41.0, 2.0), (2, 1))
    pres_atm = height_to_pressure(alt) / 1013.25
    savemat(folder / "ace_v52_2007d001.mat", {
        "Alt": alt, "Lat": np.full_like(alt, 30.0), "Lon": np.full_like(alt, 10.0),
        "Time": np.full_like(alt, T0 + 0.5), "Temp": np.full_like(alt, 230.0),
        "Pres": pres_atm})
    return pres_atm


def test_ace_pressure_converted_from_atmospheres(tmp_path):
    pres_atm = _ace(tmp_path)
    ds = load_limb_sounders(D1, "ACE", file_path=tmp_path, original_z=True)
    np.testing.assert_allclose(ds.pressure.values, pres_atm * 1013.25, rtol=1e-15)
    assert ds.pressure.attrs["units"] == "hPa"
    assert "atmospheres" in ds.pressure.attrs["comment"]
    validate(ds, "profile")


# GNSS: six profiles of 4 levels; the occultation ids and centres are
#   0 ucar  A (positions, good time, PTRes 0.1)   1 jpl A (no positions, 0.4)
#   2 jpl   B (no positions, hour time, 0.4)       3 romsaf B (positions, 1.2)
#   4 jpl   C (alone: no positions, time 03:00 but id 03:37)
#   5 romsaf D (fill value in Temp and Pres at level 2)
MSL = np.array([[10.0, 20.0, 30.0, 40.0]])
IDS = ["cosmic1c1-G01-200701010110", "cosmic1c1-G01-200701010110",
       "champ-G02-200701010205", "champ-G02-200701010205",
       "champ-G04-200701010337", "metopa-G05-200701010412"]
CENTRES = ["ucar", "jpl", "jpl", "romsaf", "jpl", "romsaf"]
FILE_TIMES = [T0 + (1 + 10.5 / 60) / 24, T0 + 1 / 24, T0 + 2 / 24,
              T0 + (2 + 5.4 / 60) / 24, T0 + 3 / 24, T0 + (4 + 12.3 / 60) / 24]
META_LAT = [10.0, 10.5, 20.0, 20.3, -40.0, 60.0]
META_LON = [100.0, 100.2, -170.0, -170.1, 5.0, 30.0]
RES = [0.1, 0.4, 0.4, 1.2, 0.4, 1.2]
FILL = np.float32(9.96921e36)


def _gnss_content():
    n, k = len(IDS), MSL.shape[1]
    lat = np.tile(np.array(META_LAT)[:, None], (1, k)).astype(np.float32)
    lon = np.tile(np.array(META_LON)[:, None], (1, k)).astype(np.float32)
    for row in (1, 2, 4):  # (JPL: no per-level positions)
        lat[row], lon[row] = np.nan, np.nan
    temp = np.tile(np.array([220.0, 215.0, 225.0, 240.0], dtype=np.float32), (n, 1))
    temp = temp + np.arange(n, dtype=np.float32)[:, None]
    pres_pa = (height_to_pressure(np.tile(MSL, (n, 1))) * 100).astype(np.float32)
    temp[5, 2], pres_pa[5, 2] = FILL, FILL
    meta = pd.DataFrame({"occid": IDS, "latitude": META_LAT, "longitude": META_LON,
                         "time": FILE_TIMES, "agency": CENTRES})
    return {"Temp": temp, "Pres": pres_pa, "Lat": lat, "Lon": lon, "MSL_alt": MSL,
            "PTRes": np.tile(np.array(RES, dtype=np.float32)[:, None], (1, k)),
            "MetaData": meta}


@pytest.fixture
def gnss(tmp_path, monkeypatch):
    """A GNSS 'file' for 2007-01-01 under tmp_path, served by a fake read_mat."""
    folder = tmp_path / "2007"
    folder.mkdir()
    (folder / "merged_ro_2007_001.mat").write_bytes(b"")
    monkeypatch.setattr(_limb_ace_gnss, "read_mat", lambda path: _gnss_content())
    return tmp_path


def _load(gnss, **options):
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # (no warning expected)
        return load_limb_sounders(D1, "GNSS", file_path=gnss, original_z=True,
                                  file_source=True, **options)


def test_gnss_one_copy_per_occultation(gnss):
    ds = _load(gnss)
    # A: ucar over jpl; B: romsaf (positions, good time) over the finer jpl;
    # C: the only copy, kept although jpl; D: kept
    assert ds.source_profile.values[:, 0].tolist() == [1, 4, 5, 6]
    assert "2 repeated copies" in ds.temperature.attrs["comment"]
    validate(ds, "profile")


def test_gnss_keep_duplicate_occultations(gnss):
    ds = _load(gnss, keep_duplicate_occultations=True)
    assert ds.source_profile.values[:, 0].tolist() == [1, 2, 3, 4, 5, 6]


def test_gnss_positions_from_metadata(gnss):
    ds = _load(gnss, keep_duplicate_occultations=True)
    for row in (1, 2, 4):
        np.testing.assert_array_equal(ds.latitude.values[row], META_LAT[row])
        np.testing.assert_array_equal(ds.longitude.values[row], META_LON[row])
    # profiles with positions in the file keep them
    np.testing.assert_array_equal(ds.latitude.values[0],
                                  np.float32(META_LAT[0]).astype(float))
    assert ds.latitude.attrs["comment"].startswith("3 profiles with no position")


def test_gnss_hour_rounded_times_take_the_id_time(gnss):
    ds = _load(gnss, keep_duplicate_occultations=True)
    got = ds.time.values[:, 0]
    want = [datenum_to_datetime(FILE_TIMES[0], "ms"),  # within 10 min of its id
            np.datetime64("2007-01-01T01:10"),         # 01:00 in the file
            np.datetime64("2007-01-01T02:05"),         # 02:00, 5 min, on the hour
            datenum_to_datetime(FILE_TIMES[3], "ms"),
            np.datetime64("2007-01-01T03:37"),
            datenum_to_datetime(FILE_TIMES[5], "ms")]
    np.testing.assert_array_equal(got, np.array(want, dtype="datetime64[us]"))
    assert "occultation id" in ds.time.attrs["comment"]


def test_gnss_pressure_in_hpa_and_fill_values_nan(gnss):
    ds = _load(gnss, keep_duplicate_occultations=True)
    want = height_to_pressure(np.tile(MSL, (6, 1)))
    want[5, 2] = np.nan
    np.testing.assert_allclose(ds.pressure.values, want, rtol=1e-6)  # (float32)
    assert np.isnan(ds.temperature.values[5, 2])
    assert ds.temperature_flag.values[5, 2] != 0  # missing, not physical range
    assert "Pa / 100" in ds.pressure.attrs["comment"]


def test_gnss_fixes_unit():
    content = _gnss_content()
    rows = _limb_ace_gnss.gnss_fixes(content, {"inst": "GNSS"})
    assert rows.tolist() == [0, 3, 4, 5]
    assert len(content["MetaData"]) == 4 and content["Lat"].shape == (4, 4)
    assert np.isnan(content["Temp"][3, 2]) and np.isnan(content["Pres"][3, 2])


# real data: pressure against the standard atmosphere, every limb instrument

REAL = [("ACE", D1), ("GNSS", D1), ("GNSS", "2020-01-10"), ("HIRDLS", D1),
        ("MLS", D1), ("SABER", D1), ("MIPAS", "2007-01-02"), ("SOFIE", "2007-06-01")]


@pytest.mark.parametrize("instrument, day", REAL)
def test_real_pressure_near_standard_atmosphere(instrument, day):
    try:
        folder = dataset_dir(instrument)
    except (DataDirNotSetError, FileNotFoundError, KeyError) as error:
        pytest.skip(f"no {instrument} data: {error}")
    if not os.path.isdir(folder):
        pytest.skip(f"no {instrument} data")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = load_limb_sounders(day, instrument, original_z=True)
    assert ds.sizes["profile"] > 0, "no profiles: data missing?"
    z, p = ds.altitude.values, ds.pressure.values
    band = (z > 12) & (z < 40) & np.isfinite(p)
    ratio = np.median(p[band] / height_to_pressure(z[band]))
    assert 1 / 3 < ratio < 3, f"{instrument}: pressure / standard atmosphere {ratio}"
