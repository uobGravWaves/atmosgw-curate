"""
Known-answer tests for the fixes from the external review of 2026-11-28 (of
b6ba9fd), outside the nadir regridding (those are in test_nadir_cris.py and
test_nadir_iasi.py): IAGOS times and altitudes, trajectory working names,
FAAM at the dateline, Aeolus azimuths and copies, DWD files out of order,
GNSS without pressure, MLS with a quality field missing, time-height
profiles, netCDF publication and 2-D nadir extras.

Each test builds the case the review reproduced, on synthetic files.
Tolerances: exact, except angles (Aeolus 1e-6 degree; the FAAM dateline
mean 1e-3 degree, as the samples are float32) and pytest.approx defaults
(relative 1e-6) for converted values.
"""

import shutil
import warnings

import netCDF4
import numpy as np
import pytest

from atmosgw.curate.data import states, to_profiles, validate
from atmosgw.curate.data.netcdf import _publish_new
from atmosgw.curate.obs import (
    load_limb_sounders,
    load_nadir_sounders,
    load_trajectories,
)
from atmosgw.curate.obs._traj_dwd_sondes import _parse, read_rows

IAGOS_NEW = "IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4"


# ----------------------------------------------------------------- IAGOS


def test_iagos_time_fill_is_not_a_time(tmp_path):
    # an old-format file's UTC_time of -9999 (its declared fill) became
    # 21:13:21 the day before departure, with a good temperature
    from test_new_loaders import iagos_file
    path = tmp_path / "IAGOS_timeseries_2019031816400001.nc"
    iagos_file(path, new=False)
    with netCDF4.Dataset(path, "a") as f:
        f.set_auto_maskandscale(False)
        f["UTC_time"][2] = -9999.0
    ds = load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path)
    t = ds.time.values[0, :50]
    assert np.isnat(t[2])
    assert np.all(t[~np.isnat(t)] >= np.datetime64("2019-03-18T16:40"))


def test_iagos_flagged_altitude_is_nan(tmp_path):
    # the flight's altitude source (GPS) is chosen per flight; a sample whose
    # GPS altitude is flagged 7 is NaN, not 10.1 km
    from test_new_loaders import iagos_file
    path = tmp_path / IAGOS_NEW
    iagos_file(path)
    with netCDF4.Dataset(path, "a") as f:
        flag = f.createVariable("gps_alt_AC_validity_flag", "i4", ("UTC_time",))
        values = np.zeros(50, dtype=int)
        values[10] = 7
        flag[:] = values
    ds = load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path)
    alt = ds.altitude.values[0]
    assert ds.altitude_source.values[0] == "gps_alt_AC"
    assert np.isnan(alt[10]) and alt[9] == pytest.approx(10.1)
    assert "flagged 3 or 7 are NaN" in ds.altitude.attrs["comment"]


def test_trajectory_working_mask_names_refused_as_extras(tmp_path):
    # an extra named TempBad replaced the temperature screening mask: every
    # temperature was removed
    from test_new_loaders import iagos_file
    iagos_file(tmp_path / IAGOS_NEW)
    with pytest.raises(ValueError, match="TempBad"):
        load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path,
                          additional_vars=["TempBad"])


# ------------------------------------------------------------------ FAAM


def test_faam_one_second_mean_across_the_dateline(tmp_path):
    from test_new_loaders import faam_file
    path = tmp_path / "core_faam_20250709_v005_r0_c999.nc"
    faam_file(path)
    with netCDF4.Dataset(path, "a") as f:
        lon = f["LON_GIN"][:]
        lon[0, :16], lon[0, 16:] = 179.9, -179.9
        f["LON_GIN"][:] = lon
    ds = load_trajectories("2025-07-09", "FAAM", file_path=tmp_path)
    lon = ds.longitude.values[0]
    assert abs(abs(lon[0]) - 180.0) < 1e-3  # (the arithmetic mean was 0)
    assert lon[1] == pytest.approx(-1.0)  # (other seconds: the plain mean)


# ---------------------------------------------------------------- Aeolus


def test_aeolus_azimuth_interpolated_across_north(tmp_path):
    from test_new_loaders import write_aeolus
    write_aeolus(tmp_path)
    path = tmp_path / "Aeolus_L2B_monthly_2019_01.nc"
    with netCDF4.Dataset(path, "a") as f:
        az = f["ray_az"][:]
        az[0], az[1] = 1.0, 359.0  # (bins 2.5 and 1.5 km of the first column)
        f["ray_az"][:] = az
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path,
                                height_scale=np.array([1.5, 2.0, 2.5]))
    az = ds.line_of_sight_azimuth.values[0]
    assert az[0] == pytest.approx(359.0, abs=1e-6)
    assert az[2] == pytest.approx(1.0, abs=1e-6)
    assert min(az[1], 360 - az[1]) < 1e-6  # (north, not 180)
    # a profile that does not cross north is untouched (where it has values:
    # its invalid observation is removed from every variable)
    other = ds.line_of_sight_azimuth.values[1]
    assert np.isfinite(other).any() and np.all(other[np.isfinite(other)] == 260.0)


def test_aeolus_copies_across_monthly_files_read_once(tmp_path):
    from test_new_loaders import write_aeolus
    write_aeolus(tmp_path)
    shutil.copy(tmp_path / "Aeolus_L2B_monthly_2019_01.nc",
                tmp_path / "Aeolus_L2B_monthly_2018_12.nc")
    ds = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path,
                            original_z=True)
    assert ds.sizes["profile"] == 2  # (4 before: each file's copies)


def test_aeolus_valid_copy_not_dropped_for_an_invalid_one(tmp_path):
    # two copies equal in every field but validity, the first invalid: the
    # valid copy is read, in its own column's profile (no extra profile)
    from test_new_loaders import write_aeolus
    write_aeolus(tmp_path, duplicate=True)
    path = tmp_path / "Aeolus_L2B_monthly_2019_01.nc"
    with netCDF4.Dataset(path, "a") as f:
        valid = f["ray_valid"][:]
        valid[0] = 0
        f["ray_valid"][:] = valid
    ds = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path,
                            original_z=True)
    hlos = ds.hlos_wind.values
    assert ds.sizes["profile"] == 2  # (as without copies)
    assert hlos[0, 2] == pytest.approx(10.0)  # (the valid copy's 1000 cm/s)
    assert not states.has_state(ds.hlos_wind_flag.values[0, 2], "instrument_flag")


# ------------------------------------------------------------------- DWD


def test_dwd_ascent_moved_to_the_end_of_a_file_is_found(tmp_path):
    # the lines of one ascent (here 60) moved after every other row: the
    # forward scan stopped at the first later row and missed them
    from test_new_loaders import write_sondes
    path = write_sondes(tmp_path)
    lines = path.read_bytes().splitlines(keepends=True)
    moved = [ln for ln in lines[1:] if b";2020011005;" in ln]
    rest = [ln for ln in lines[1:] if b";2020011005;" not in ln]
    path.write_bytes(lines[0] + b"".join(rest + moved))
    rows = read_rows(path, 2020011005, 2020011005)
    assert rows.shape[0] == len(moved) == 60
    np.testing.assert_array_equal(rows, _parse(moved))
    ds = load_trajectories("2020-01-10", "DWDsondes", file_path=tmp_path)
    assert ds.sizes["flight"] == 2


# ------------------------------------------------------------------ GNSS


def test_gnss_without_pressure_estimates_it_from_msl_alt(tmp_path, monkeypatch):
    from test_limb_gnss_ace import _gnss_content

    from atmosgw.curate.obs import _limb_ace_gnss
    folder = tmp_path / "2007"
    folder.mkdir()
    (folder / "merged_ro_2007_001.mat").write_bytes(b"")
    content = _gnss_content()
    content.pop("Pres")
    monkeypatch.setattr(_limb_ace_gnss, "read_mat", lambda path: content)
    ds = load_limb_sounders("2007-01-01", "GNSS", file_path=tmp_path, original_z=True)
    validate(ds, "profile")
    assert ds.pressure.attrs["origin"] == "derived"
    assert np.all(np.isfinite(ds.pressure.values))


# ------------------------------------------------------------------- MLS


def _drop_field(path, name):
    import h5py
    with h5py.File(path, "a") as f:
        del f[f"HDFEOS/SWATHS/Temperature/Data Fields/{name}"]


def test_mls_missing_field_screens_the_rest_and_marks_incomplete(tmp_path):
    # Convergence missing: Status, Quality, precision and pressure rules are
    # still applied; every value of the file is qc_incomplete (kept)
    from test_limb_selection import _mls
    _mls(tmp_path, status=(1, 0, 0))
    _drop_field(tmp_path / "2026" / "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d001.he5",
                "Convergence")
    ds = load_limb_sounders("2026-01-01", "MLS", file_path=tmp_path, original_z=True)
    flag = ds.temperature_flag.values
    assert states.has_state(flag[0], "instrument_flag").all()  # (Status odd)
    assert not states.has_state(flag[1, 2:6], "instrument_flag").any()
    # (a qualified state: on every value that exists, none removed for it)
    kept = np.isfinite(ds.temperature.values)
    assert kept.any() and states.has_state(flag[kept], "qc_incomplete").all()
    assert np.all(ds.temperature.values[1, 2:6] == 230)
    assert "missing from a file (Convergence)" in ds.temperature.attrs["comment"]


@pytest.mark.parametrize("missing", ["Status", "Quality", "Convergence",
                                     "L2gpPrecision"])
def test_mls_each_remaining_rule_applied(tmp_path, missing):
    # one field missing; each of the other profile or level rules still marks
    # its failure: profile 0 odd Status, 1 Quality 0.5 (bad at 100 hPa and
    # below), 2 Convergence 1.04; precision -1 at profile 1, 10 hPa (where
    # nothing else fails)
    from test_limb_selection import _mls
    precision = np.ones((3, 7), np.float32)
    precision[1, 4] = -1.0
    _mls(tmp_path, status=(1, 0, 0), quality=(1.5, 0.5, 1.5),
         convergence=(1.0, 1.0, 1.04), precision=precision)
    _drop_field(tmp_path / "2026" / "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d001.he5",
                missing)
    ds = load_limb_sounders("2026-01-01", "MLS", file_path=tmp_path, original_z=True)
    bad = states.has_state(ds.temperature_flag.values, "instrument_flag")
    # (levels 2-5: 100, 82.5, 10, 0.001 hPa, inside the useful range)
    expected = {"Status": bad[0, 2:6].all(), "Quality": bad[1, 2],
                "Convergence": bad[2, 2:6].all(), "L2gpPrecision": bad[1, 4]}
    for rule, applied in expected.items():
        assert applied == (rule != missing), rule
    assert not bad[1, [3, 5]].any()  # (Quality 0.5 is good above 83 hPa)


def test_mls_incomplete_marked_only_for_the_incomplete_file(tmp_path):
    import h5py
    from test_limb_selection import _mls
    _mls(tmp_path)
    one = tmp_path / "2026" / "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d001.he5"
    two = tmp_path / "2026" / "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d002.he5"
    shutil.copy(one, two)
    with h5py.File(two, "a") as f:
        times = f["HDFEOS/SWATHS/Temperature/Geolocation Fields/Time"]
        times[...] = times[()] + 86400
    _drop_field(two, "Quality")
    ds = load_limb_sounders(["2026-01-01", "2026-01-02T23:59"], "MLS",
                            file_path=tmp_path, original_z=True)
    incomplete = states.has_state(ds.temperature_flag.values, "qc_incomplete")
    kept = np.isfinite(ds.temperature.values)
    # (the day from level 2: 1000 hPa is outside the useful range, removed in
    # every variable, time too)
    day = np.broadcast_to(ds.time.values[:, 2:3].astype("datetime64[D]"), kept.shape)
    second = day == np.datetime64("2026-01-02")
    assert (kept & second).any() and incomplete[kept & second].all()
    assert not incomplete[~second].any()


# ------------------------------------------------------ time-height profiles


def test_time_height_profiles_keep_versions_and_copy(tmp_path):
    from test_time_height_known import sessions_2007, write_year

    from atmosgw.curate.obs import load_time_height
    write_year(tmp_path, 2007, sessions_2007())
    series = load_time_height(["2007-01-24", "2007-01-25"], "ALOMAR",
                              file_path=tmp_path, file_source=True)
    p = to_profiles(series)
    validate(p, "profile")
    assert list(p.original_file_version.values) == list(
        series.original_file_version.values)
    before = series.session.values[0]
    p["session"].values[0] = before + 99
    assert series.session.values[0] == before  # (not shared)
    for name in ("temperature", "pressure"):
        refs = p[name].attrs.get("ancillary_variables", "").split()
        assert all(r in p for r in refs)


# ---------------------------------------------------------------- netCDF


def test_file_saved_before_a_new_state_still_valid(tmp_path):
    # flag variables carrying the registry without its last state (as files
    # saved before qc_incomplete was added) pass validate; a wrong one fails
    from test_nadir_cris import write_cris

    from atmosgw.curate.data import FormatError
    write_cris(tmp_path)
    ds = load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                             granules=[("2019-01-01", 1)], original_grid=True)
    attrs = states.flag_attrs()
    old = {"flag_masks": np.asarray(attrs["flag_masks"])[:-1],
           "flag_meanings": " ".join(attrs["flag_meanings"].split()[:-1])}
    flags = [n for n in ds.data_vars if n.endswith("_flag")]
    for name in flags:
        ds[name].attrs.update(old)
    validate(ds, "swath")
    ds[flags[0]].attrs["flag_meanings"] = "x " + old["flag_meanings"]
    with pytest.raises(FormatError, match="registry"):
        validate(ds, "swath")


def test_new_file_not_replaced_by_a_second_writer(tmp_path):
    # the second of two writers racing for one new path gets FileExistsError
    # and the first's file stays (before: both succeeded, the last won)
    path = tmp_path / "out.nc"
    first, second = tmp_path / ".a.partial", tmp_path / ".b.partial"
    first.write_bytes(b"first")
    second.write_bytes(b"second")
    _publish_new(first, path)
    with pytest.raises(FileExistsError):
        _publish_new(second, path)
    assert path.read_bytes() == b"first" and not first.exists()


# ------------------------------------------- file longitudes (C. Wright)


def test_file_longitudes_0_360_accepted_outside_bad_in_every_loader(tmp_path):
    # file longitudes must lie in -180..360: 180-360 are wrapped and kept,
    # beyond is bad_geolocation (Q13, limb; nadir and trajectories had
    # marked 180-360 bad). The requested lon_range takes any value (unchanged).
    from test_nadir_iasi import DAY_SECONDS, write_orbit
    from test_new_loaders import faam_file
    path = write_orbit(tmp_path / "iasi", DAY_SECONDS, 90)
    with netCDF4.Dataset(path, "a") as f:
        f["l1_lon"][:] = f["l1_lon"][:] + 180.0  # (200-218 E)
    ds = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 1)],
                             file_path=tmp_path / "iasi", screen="none")
    lon = ds.longitude.values[np.isfinite(ds.temperature.values[..., 0])]
    assert lon.size and np.all((lon > -161) & (lon < -141))
    boxed = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 1)],
                                file_path=tmp_path / "iasi", screen="none",
                                lon_range=(190, 230))  # (asked for in 0-360)
    assert boxed.sizes["along_track"] > 0
    with netCDF4.Dataset(path, "a") as f:  # (one scan beyond 360: bad)
        lon = f["l1_lon"][:]
        lon[10] = lon[10] + 200.0
        f["l1_lon"][:] = lon
    native = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 1)],
                                 file_path=tmp_path / "iasi", screen="none",
                                 original_grid=True)
    flag = native.temperature_flag.values[..., 0]
    assert states.has_state(flag[10], "bad_geolocation").all()
    assert not states.has_state(flag[11], "bad_geolocation").any()
    faam = tmp_path / "faam" / "core_faam_20250709_v005_r0_c999.nc"
    faam.parent.mkdir()
    faam_file(faam)
    with netCDF4.Dataset(faam, "a") as f:
        lon = f["LON_GIN"][:]
        lon[0], lon[1] = 359.0, 400.0  # (second 0: -1 E; second 1: out of range)
        f["LON_GIN"][:] = lon
    ds = load_trajectories("2025-07-09", "FAAM", file_path=faam.parent)
    assert ds.longitude.values[0, 0] == pytest.approx(-1.0)
    assert states.has_state(ds.temperature_flag.values[0, 1], "bad_geolocation")
    # one bad sample among good ones: the second is bad, not a plausible mean
    with netCDF4.Dataset(faam, "a") as f:
        lon = f["LON_GIN"][:]
        lon[2, :] = 10.0
        lon[2, 5] = 400.0
        f["LON_GIN"][:] = lon
    ds = load_trajectories("2025-07-09", "FAAM", file_path=faam.parent)
    assert states.has_state(ds.temperature_flag.values[0, 2], "bad_geolocation")


# ------------------------------------------------------------ nadir 2-D


def test_nadir_2d_extras_refused_clearly(tmp_path):
    from test_nadir_cris import write_cris
    write_cris(tmp_path)
    with pytest.raises(TypeError, match="additional_vars"):
        load_nadir_sounders(None, "CrIS", product="2d", file_path=tmp_path,
                            granules=[("2019-01-01", 1)], additional_vars=["bt_4mu"])
