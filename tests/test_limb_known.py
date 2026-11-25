"""Known-answer tests for load_limb_sounders, on synthetic files.

Port of the known-answer checks in the MATLAB test_get_limbsounders.m (section
1), run through the public function. The synthetic files are written here by
Python, with the same recipe as the MATLAB test (so they need no data):
2007-01-01 and 2007-01-02, 8 profiles x 21 levels per day; profile k (1-8) has
heights 10:2:50 km + 0.3k, Temp = 200 + 1.5 * Alt (linear, so interpolation
is exact), Lat linspace(-60, 60, 8), Lon linspace(-157.5, 157.5, 8), and time
(k - 0.5)/8 of the day.

Differences from the MATLAB checks: the out-of-range date check uses ACE (the
MATLAB one uses HIRDLS, whose reader is not ported yet); the check of the
MATLAB helper list_non_modal_size has no Python counterpart.

Tolerances: 1e-9 (K, km, degrees) for values that are exact in exact
arithmetic, as in the MATLAB test; times to 1 ms (output times are rounded to
the millisecond).
"""

import warnings

import numpy as np
import pytest
from scipy.io import savemat

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_limb_sounders

HS = np.arange(18, 60.5, 0.5)
D1, D2 = "2007-01-01", "2007-01-02"
TOP = 50 + 0.3 * np.arange(1, 9)  # top of each profile, km
EXP_TEMP = 200 + 1.5 * np.tile(HS, (8, 1))
EXP_FINITE = np.tile(HS, (8, 1)) <= TOP[:, None]
LEV = np.arange(10, 52, 2, dtype=float)


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    """Write the synthetic files; return their folder."""
    folder = tmp_path_factory.mktemp("synthetic")
    (folder / "sampler").mkdir()
    n_prof, n_lev = 8, LEV.size
    d0 = datenum(2007, 1, 1)
    for day in (0, 1):
        tag = f"2007d{day + 1:03d}"
        alt = 0.3 * np.arange(1, n_prof + 1)[:, None] + LEV
        s = {
            "Alt": alt,
            "Lat": np.tile(np.linspace(-60, 60, n_prof)[:, None], (1, n_lev)),
            "Lon": np.tile(np.linspace(-157.5, 157.5, n_prof)[:, None], (1, n_lev)),
            "Time": np.tile(d0 + day + (np.arange(1, n_prof + 1)[:, None] - 0.5)
                            / n_prof, (1, n_lev)),
            "Temp": 200 + 1.5 * alt,
        }
        savemat(folder / f"synlin_{tag}.mat", s)

        b = {k: v.copy() for k, v in s.items()}
        if day == 0:
            b["Lat"][2, 4] = 95
            b["Temp"][5, 6] = 50
        savemat(folder / f"synbad_{tag}.mat", b)

        savemat(folder / f"syn3d_{tag}.mat",
                {**s, "Extra3": s["Temp"][:, :, None] * np.arange(1, 4)})

        p = {k: v for k, v in s.items() if k != "Temp"}
        p["Temp_Residual"] = 0.1 * alt
        p["Temp_PW"] = alt[:, :, None] * np.arange(1, 5)
        if day == 0:
            p["Lat"] = p["Lat"].copy()
            p["Lat"][2, 4] = 95
        p["Meta"] = np.array([[7.0]])
        savemat(folder / f"synpw_{tag}.mat", p)

        x = dict(s)
        x["PerPointT"] = s["Temp"].T
        x["PerProf"] = 10.0 * np.arange(1, n_prof + 1)[:, None]
        x["PerLev"] = (LEV + 0.5)[None, :]
        x["Meta"] = np.array([[7.0]])
        x["TempMeta"] = np.array([[3.0]])
        x["FlagVar"] = np.tile(np.mod(np.arange(n_lev), 3).astype(float), (n_prof, 1))
        savemat(folder / f"synext_{tag}.mat", x)

        w = dict(s)
        w["Lon"] = np.tile(np.array([190, -190, 350, -350, 0, 10, 20, 30.0])[:, None],
                           (1, n_lev))
        savemat(folder / f"synlon_{tag}.mat", w)

        if day == 0:
            sampled = {"Alt": alt.T, "Lat": s["Lat"].T, "Lon": s["Lon"].T,
                       "Time": s["Time"].T, "T": s["Temp"].T,
                       "Prs": height_to_pressure(alt).T,
                       "SourceProf": np.tile(np.arange(1, n_prof + 1.0), (n_lev, 1))}
            savemat(folder / "sampler" / f"sampled_{int(d0)}.mat",
                    {"Sampled_Data": sampled})
    return folder


def misc(syn, identifier, time_range=D1, **options):
    return load_limb_sounders(time_range, "Misc", misc_path=syn,
                              misc_identifier=identifier, **options)


def n_profiles(ds):
    return ds.sizes["profile"]


# linear profiles


def test_linear_temperature_exact_and_nan_above_top(syn):
    ds = misc(syn, "synlin")
    temp = ds.temperature.values
    assert temp.shape == (8, HS.size)
    np.testing.assert_array_equal(np.isfinite(temp), EXP_FINITE)
    assert np.max(np.abs(temp[EXP_FINITE] - EXP_TEMP[EXP_FINITE])) < 1e-9


def test_output_altitude_is_height_scale(syn):
    np.testing.assert_array_equal(misc(syn, "synlin").altitude.values,
                                  np.tile(HS, (8, 1)))


def test_position_and_time_constant_within_profiles(syn):
    ds = misc(syn, "synlin")
    lat, lon = ds.latitude.values, ds.longitude.values
    assert np.all(np.nanmax(lat, 1) - np.nanmin(lat, 1) == 0)
    assert np.all(np.nanmax(lon, 1) - np.nanmin(lon, 1) < 1e-12)
    for row in ds.time.values:  # NaT above the top of each profile
        valid = row[~np.isnat(row)]
        assert valid.max() - valid.min() <= np.timedelta64(1, "ms")


def test_column_height_scale_and_lower_case_name(syn):
    # identical apart from the 'options' attribute, which records the call as given
    def without_options(ds):
        ds.attrs.pop("options")
        return ds

    ds = without_options(misc(syn, "synlin"))
    assert without_options(misc(syn, "synlin", height_scale=HS[:, None])).identical(ds)
    assert without_options(load_limb_sounders(D1, "misc", misc_path=syn,
                                              misc_identifier="synlin")).identical(ds)


def test_q6_pressure_interpolated_in_log_pressure(syn):
    # Q6 (fixed): pressure is interpolated in log(pressure), so between two
    # levels it is the exponential through them. Tolerance 1e-14 relative:
    # a few roundings in log, the interpolation and exp
    ds = misc(syn, "synlin")
    alt = 0.3 * 1 + LEV  # profile 1's original heights
    p = height_to_pressure(alt)
    h = 19.0  # between profile 1's levels at 18.3 and 20.3 km
    j = np.searchsorted(alt, h) - 1
    t = (h - alt[j]) / (alt[j + 1] - alt[j])
    log_interp = np.exp((1 - t) * np.log(p[j]) + t * np.log(p[j + 1]))
    got = ds.pressure.values[0, np.flatnonzero(HS == h)[0]]
    assert got == pytest.approx(log_interp, rel=1e-14)
    linear = (1 - t) * p[j] + t * p[j + 1]
    assert linear > got * (1 + 5e-4)  # the old (linear) value was ~0.1% higher


def test_q6_nearest_keeps_pressure_values(syn):
    # with 'nearest', pressure is not put through log/exp: values are exactly
    # those of the data
    ds = misc(syn, "synlin", interp_type="nearest")
    p = height_to_pressure(0.3 * 1 + LEV)
    got = ds.pressure.values[0]
    assert np.all(np.isin(got[np.isfinite(got)], p))


# time handling (profile k of each day is at (k - 0.5)/8 of the day)

TIME_CASES = [
    (D1, "fuzzy", 8),
    ("2007-01-01T12:00", "fuzzy", 8),
    ("2007-01-01T12:00", "generous", 8),
    ([ "2007-01-01T06:00", "2007-01-01T18:00"], "strict", 4),
    (["2007-01-01T12:00", "2007-01-02T02:24"], "generous", 16),
    ([D1, D2], "fuzzy", 16),
    ([D2, D1], "fuzzy", 16),
    ([D1, D2], "strict", 8),
]


@pytest.mark.parametrize("time_range,mode,expected", TIME_CASES)
def test_time_handling(syn, time_range, mode, expected):
    assert n_profiles(misc(syn, "synlin", time_range, time_handling=mode)) == expected


def test_strict_single_time_warns_and_returns_nothing(syn):
    with pytest.warns(UserWarning, match="exactly that instant"):
        ds = misc(syn, "synlin", time_handling="strict")
    assert n_profiles(ds) == 0


# lat/lon selection (lats linspace(-60, 60, 8), lons linspace(-157.5, 157.5, 8))


def test_lat_lon_box(syn):
    assert n_profiles(misc(syn, "synlin", lat_range=(0, 90))) == 4
    assert n_profiles(misc(syn, "synlin", lon_range=(-180, 0))) == 4
    assert n_profiles(misc(syn, "synlin", lat_range=(0, 90), lon_range=(-180, 0))) == 0


# provenance


def test_source_profile_file_and_names_over_two_days(syn):
    ds = misc(syn, "synlin", [D1, D2], file_source=True)
    np.testing.assert_array_equal(ds.source_profile.values[:, 0],
                                  np.r_[1:9, 1:9].astype(float))
    np.testing.assert_array_equal(ds.source_file.values[:, 0], [1.0] * 8 + [2.0] * 8)
    assert list(ds.original_files.values) == ["synlin_2007d001.mat",
                                              "synlin_2007d002.mat"]


# outliers (day 1 of synbad: Lat[2, 4] = 95, Temp[5, 6] = 50 K; 0-based)


def test_outliers_nan_in_every_field(syn):
    ds = misc(syn, "synbad", original_z=True)
    for name in ("latitude", "longitude", "temperature", "pressure", "altitude"):
        v = ds[name].values
        assert np.isnan(v[2, 4]) and np.isnan(v[5, 6]) and np.isnan(v).sum() == 2, name
    t = ds.time.values
    assert np.isnat(t[2, 4]) and np.isnat(t[5, 6]) and np.isnat(t).sum() == 2


def test_keep_outliers(syn):
    ds = misc(syn, "synbad", original_z=True, keep_outliers=True)
    assert ds.latitude.values[2, 4] == 95 and ds.temperature.values[5, 6] == 50


# original heights clipped to a height range


def test_original_z_height_range(syn):
    ds = misc(syn, "synlin", original_z=True, height_range=(20, 40))
    alt_full = 0.3 * np.arange(1, 9)[:, None] + LEV
    alt = ds.altitude.values
    in_range = (alt >= 20) & (alt <= 40)
    assert alt.shape[1] == np.sum(np.any((alt_full >= 20) & (alt_full <= 40), axis=0))
    for name in ("latitude", "longitude", "temperature", "pressure"):
        np.testing.assert_array_equal(np.isfinite(ds[name].values), in_range)
    np.testing.assert_array_equal(~np.isnat(ds.time.values), in_range)
    temp = ds.temperature.values
    assert np.max(np.abs(temp[in_range] - (200 + 1.5 * alt[in_range]))) < 1e-9


# planetary-wave data, with a 3-D temperature_pw (mode k is k * Alt)


def pw(syn, time_range=D1, **options):
    return misc(syn, "synpw", time_range, hindley_pws=True, **options)


def test_pw_outlier_in_all_modes(syn):
    pw_values = pw(syn, original_z=True).temperature_pw.values
    assert np.all(np.isnan(pw_values[2, 4, :])) and np.isnan(pw_values).sum() == 4


def test_pw_original_z_height_range_masks_all_modes(syn):
    ds = pw(syn, original_z=True, height_range=(20, 40))
    alt = ds.altitude.values
    in_range = (alt >= 20) & (alt <= 40) & np.isfinite(ds.latitude.values)
    for k in range(4):
        np.testing.assert_array_equal(np.isfinite(ds.temperature_pw.values[:, :, k]),
                                      in_range)


def test_pw_interpolated_exactly(syn):
    ds = pw(syn, D2)
    assert ds.temperature_pw.dims == ("profile", "level", "mode")
    assert ds.temperature_pw.shape == (8, HS.size, 4)
    expected = np.tile(HS, (8, 1))[:, :, None] * np.arange(1, 5)
    err = np.abs(ds.temperature_pw.values - expected)
    assert np.nanmax(err) < 1e-9


# additional variables


def test_3d_additional_var_over_two_days(syn):
    ds = misc(syn, "syn3d", [D1, D2], additional_vars=["Extra3"])
    extra = ds.Extra3.values
    assert extra.shape == (16, HS.size, 3)
    expected = ds.temperature.values[:, :, None] * np.arange(1, 4)
    np.testing.assert_array_equal(np.isfinite(extra), np.isfinite(expected))
    assert np.nanmax(np.abs(extra - expected)) < 1e-9


def test_additional_var_shapes(syn):
    ds = misc(syn, "synext", [D1, D2], additional_vars=["PerPointT", "PerProf",
                                                         "PerLev"])
    temp = ds.temperature.values
    k = np.tile(np.r_[1:9, 1:9][:, None], (1, HS.size))
    np.testing.assert_array_equal(ds.PerPointT.values, temp)  # transposed in file
    for name, expected in (("PerProf", 10.0 * k),
                           ("PerLev", np.tile(HS, (16, 1)) - 0.3 * k + 0.5)):
        v = ds[name].values
        np.testing.assert_array_equal(np.isfinite(v), np.isfinite(temp))
        assert np.nanmax(np.abs(v - expected)) < 1e-9, name


def test_scalar_additional_var_is_refused(syn):
    with pytest.raises(ValueError, match="ignore_non_modal"):
        misc(syn, "synext", additional_vars=["Meta"])


def test_missing_additional_var_is_refused(syn):
    with pytest.raises(KeyError):
        misc(syn, "synext", additional_vars=["Nope"])


def test_ignore_non_modal_passes_metadata_through(syn):
    ds = misc(syn, "synext", additional_vars=["Meta", "TempMeta"],
              ignore_non_modal=True)
    assert ds.Meta.values == 7 and ds.TempMeta.values == 3
    assert ds.Meta.dims == ()
    np.testing.assert_array_equal(ds.temperature.values,
                                  misc(syn, "synlin").temperature.values)


def test_ignore_non_modal_with_pw(syn):
    ds = pw(syn, D2, additional_vars=["Meta"], ignore_non_modal=True)
    assert ds.altitude.shape == (8, HS.size)
    np.testing.assert_array_equal(ds.temperature_pw.values,
                                  pw(syn, D2).temperature_pw.values)
    assert ds.Meta.values == 7


def test_integer_flag_interpolated_with_nearest(syn):
    ds = misc(syn, "synext", additional_vars=["FlagVar"])
    flag = ds.FlagVar.values
    np.testing.assert_array_equal(np.isfinite(flag), np.isfinite(ds.temperature.values))
    alt_k = 0.3 * np.arange(1, 9)[:, None] + LEV
    flag_k = np.tile(np.mod(np.arange(LEV.size), 3), (8, 1))
    compared = 0
    for p in range(8):
        for h in np.flatnonzero(np.isfinite(flag[p])):
            dist = np.abs(alt_k[p] - HS[h])
            order = np.argsort(dist, kind="stable")
            if dist[order[1]] - dist[order[0]] > 1e-6:  # not a near-tie
                assert flag[p, h] == flag_k[p, order[0]]
                compared += 1
    assert compared > 100
    assert set(np.unique(flag[np.isfinite(flag)])) <= {0.0, 1.0, 2.0}


# longitudes


def test_longitudes_wrapped_not_removed(syn):
    expected = np.array([-170, 170, -10, 10, 0, 10, 20, 30.0])
    ds = misc(syn, "synlon", original_z=True)
    np.testing.assert_array_equal(ds.longitude.values,
                                  np.tile(expected[:, None], (1, ds.sizes["level"])))
    assert not np.isnan(ds.temperature.values).any()
    ds = misc(syn, "synlon")
    assert np.nanmax(np.abs(ds.longitude.values - expected[:, None])) < 1e-9


# sampler output


def test_sampler(syn):
    ds = load_limb_sounders(D1, "Sampler", sample_path=syn / "sampler")
    temp = ds.temperature.values
    np.testing.assert_array_equal(np.isfinite(temp), EXP_FINITE)
    assert np.max(np.abs(temp[EXP_FINITE] - EXP_TEMP[EXP_FINITE])) < 1e-9


# dates outside the instrument's record


def test_out_of_range_date_warns_and_returns_nothing(syn):
    with pytest.warns(UserWarning, match="only available"):
        ds = load_limb_sounders("2001-01-01", "ACE", file_path=syn)
    assert n_profiles(ds) == 0


def test_no_data_root_needed_with_explicit_paths(syn, monkeypatch, tmp_path):
    # as on a machine (e.g. CI) with no data root set: file_path, misc_path and
    # sample_path must work without one
    monkeypatch.delenv("ATMOSGW_DATA", raising=False)
    monkeypatch.setenv("ATMOSGW_CONFIG", str(tmp_path / "no_such_config.toml"))
    with pytest.warns(UserWarning, match="only available"):
        assert n_profiles(load_limb_sounders("2001-01-01", "ACE", file_path=syn)) == 0
    assert n_profiles(misc(syn, "synlin")) == 8
    assert n_profiles(load_limb_sounders(D1, "Sampler",
                                         sample_path=syn / "sampler")) == 8


# output format and arguments


def test_output_format(syn):
    ds = misc(syn, "synlin", [D1, D2], file_source=True, dx=True)
    assert ds.temperature.dims == ("profile", "level")
    assert ds.time.dtype == np.dtype("datetime64[us]")
    assert ds.time.values[0, 0] == np.datetime64("2007-01-01T01:30")
    assert ds.temperature.attrs["units"] == "K"
    assert ds.pressure.attrs["units"] == "hPa"
    assert ds.attrs["instrument"] == "Misc"
    assert np.all(np.isnan(ds.dx.values[0]))
    # neighbouring profiles: 17.14 degrees of latitude and 45 of longitude apart
    assert 1000 < np.nanmin(ds.dx.values[1:8]) < 6000


def test_no_source_fields_unless_asked(syn):
    ds = misc(syn, "synlin")
    assert "source_profile" not in ds and "original_files" not in ds


def test_arguments_checked(syn):
    with pytest.raises(TypeError, match="datenum"):
        misc(syn, "synlin", 733043.0)
    with pytest.raises(ValueError, match="AIRS"):
        load_limb_sounders(D1, "AIRS")
    with pytest.raises(ValueError, match="Unknown instrument"):
        load_limb_sounders(D1, "SAB")  # no abbreviations
    with pytest.raises(TypeError, match="original_z"):
        misc(syn, "synlin", original_z=1)
    with pytest.raises(TypeError, match="additional_vars"):
        misc(syn, "synlin", additional_vars="Meta")
    with pytest.raises(ValueError, match="time_handling"):
        misc(syn, "synlin", time_handling=3)
    with pytest.raises(ValueError, match="1979"):
        misc(syn, "synlin", "1970-01-01")
    with pytest.raises(ValueError, match="misc_path"):
        load_limb_sounders(D1, "Misc")


# real data (2007-01-01, as the MATLAB test); skipped if absent


def _hirdls_or_skip():
    from atmosgw.curate.paths import DataDirNotSetError, dataset_dir

    try:
        if dataset_dir("HIRDLS", must_exist=False).is_dir():
            return
    except DataDirNotSetError:
        pass
    pytest.skip("no HIRDLS data")


def test_hirdls_extras_shapes():
    _hirdls_or_skip()
    ds = load_limb_sounders(D1, "HIRDLS", lat_range=(-20, 20),
                            additional_vars=["TemperaturePrecision",
                                             "OrbitAscendingFlag"])
    assert ds.TemperaturePrecision.shape == ds.temperature.shape
    flag = ds.OrbitAscendingFlag.values
    assert flag.shape == ds.temperature.shape
    spread = np.nanmax(flag, axis=1) - np.nanmin(flag, axis=1)
    assert np.all(spread[np.isfinite(spread)] == 0)  # one value per profile


def test_hirdls_no_fill_values():
    _hirdls_or_skip()
    ds = load_limb_sounders(D1, "HIRDLS", original_z=True, keep_outliers=True,
                            lat_range=(-20, 20),
                            additional_vars=["TemperaturePrecision"])
    assert not np.any(ds.temperature.values == -999)
    assert not np.any(ds.TemperaturePrecision.values == -999)
    assert not np.any(ds.altitude.values < 0)


def _data_or_skip(name):
    from atmosgw.curate.paths import DataDirNotSetError, dataset_dir

    try:
        if dataset_dir(name, must_exist=False).is_dir():
            return
    except DataDirNotSetError:
        pass
    pytest.skip(f"no {name} data")


@pytest.mark.parametrize("instrument,date,point_var,profile_var", [
    ("MIPAS", D2, "TEM_SD", "sza"),
    ("MLS", D1, "L2gpPrecision", "Quality"),
])
def test_extras_shapes(instrument, date, point_var, profile_var):
    _data_or_skip(instrument)
    ds = load_limb_sounders(date, instrument, lat_range=(-20, 20),
                            additional_vars=[point_var, profile_var])
    assert ds[point_var].shape == ds.temperature.shape
    value = ds[profile_var].values
    assert value.shape == ds.temperature.shape
    spread = np.nanmax(value, axis=1) - np.nanmin(value, axis=1)
    assert np.all(spread[np.isfinite(spread)] == 0)  # one value per profile


def test_mipas_missing_additional_var_is_refused():
    _data_or_skip("MIPAS")
    with pytest.raises(KeyError, match="Nope"):
        load_limb_sounders(D2, "MIPAS", additional_vars=["Nope"])


def test_alomar_file_names_and_year_without_file():
    _data_or_skip("ALOMAR")
    ds = load_limb_sounders("2007-01-23", "ALOMAR", file_source=True)
    assert list(ds.original_files.values) == ["2007_V2026-07-27_file4juli.h5"]
    ds = load_limb_sounders("2022-06-01", "ALOMAR")
    assert ds.temperature.shape == (0, HS.size)


def test_q7_alomar_each_day_has_its_own_heights():
    # Q7 (fixed): each day's own heights are used. In 2007 the grid starts at
    # 5.25 km (765 levels) until 2007-05-14, then at 5.10 km (766 levels); the
    # first day's heights used to be applied to every day. Checked against the
    # fixed MATLAB module (MultiInstrumentFunctions e13e247), 2026-11-24.
    _data_or_skip("ALOMAR")
    early = load_limb_sounders("2007-01-23", "ALOMAR", original_z=True,
                               keep_outliers=True)
    late = load_limb_sounders("2007-06-01", "ALOMAR", original_z=True,
                              keep_outliers=True)
    assert early.altitude.values[0, 0] == pytest.approx(5.25, abs=1e-9)
    assert late.altitude.values[0, 0] == pytest.approx(5.10, abs=1e-9)
    assert np.isfinite(late.altitude.values[0]).sum() == 766


def test_airs_and_iasi_not_ported():
    for name in ("AIRS", "IASI"):
        with pytest.raises(ValueError, match="AIRS and IASI are not ported"):
            load_limb_sounders(D1, name)


def test_verbose_prints(syn, capsys):
    misc(syn, "synlin", verbose=True)
    out = capsys.readouterr().out
    assert "Loading data" in out
    assert "------> Temp (" in out  # one line per variable interpolated
    misc(syn, "synlin")
    assert capsys.readouterr().out == ""


def test_no_warnings_in_normal_use(syn):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        misc(syn, "synlin", [D1, D2], file_source=True, dx=True)


# which file is used when several match a day (Q8), and the output's record
# of the load (attributes). New checks, no MATLAB counterpart.


def test_newest_mls_version():
    from atmosgw.curate.obs._limb_tools import newest_mls_version

    def name(version):
        return f"/x/MLS-Aura_L2GP-Temperature_{version}_2026d001.he5"

    assert newest_mls_version([name("v05-04-c01"), name("v06-04-c01")]) == \
        name("v06-04-c01")
    assert newest_mls_version([name("v06-04-c01"), name("v06-03-c02")]) == \
        name("v06-04-c01")
    assert newest_mls_version([name("v06-01-c01"), name("v06-01-c02")]) == \
        name("v06-01-c02")
    assert newest_mls_version(["/x/odd_2026d001.he5", name("v05-04-c01")]) == \
        name("v05-04-c01")  # no version in the name: counts as oldest
    assert newest_mls_version(["/x/a_2026d001.he5", "/x/b_2026d001.he5"]) == \
        "/x/a_2026d001.he5"  # equal: the first, as MATLAB


def test_choose_file_records_several_files_not_folders(tmp_path):
    from atmosgw.curate.obs._limb_tools import choose_file

    (tmp_path / "x_2007d001").mkdir()
    for name in ("a_2007d001.mat", "b_2007d001.mat"):
        (tmp_path / name).write_bytes(b"")
    folder = str(tmp_path / "x_2007d001") + "/"
    a, b = str(tmp_path / "a_2007d001.mat"), str(tmp_path / "b_2007d001.mat")
    info = {}
    assert choose_file([], info, "day") is None
    assert choose_file([folder, a], info, "day") == folder  # first, as MATLAB
    assert "record" not in info  # only one file: nothing recorded
    assert choose_file([a, b], info, "2007-01-01") == a
    assert info["record"]["multiple_matches"] == [
        "2007-01-01: a_2007d001.mat, b_2007d001.mat -> used a_2007d001.mat"]


def _write_mls(path, temperature):
    """A minimal MLS Level 2 temperature file: 3 profiles x 5 pressure levels."""
    import h5py

    group = "HDFEOS/SWATHS/Temperature"
    with h5py.File(path, "w") as f:
        f[f"{group}/Data Fields/L2gpValue"] = np.full((3, 5), temperature,
                                                      dtype=np.float32)
        geo = f"{group}/Geolocation Fields"
        f[f"{geo}/Latitude"] = np.array([-10, 0, 10], dtype=np.float32)
        f[f"{geo}/Longitude"] = np.array([0, 1, 2], dtype=np.float32)
        f[f"{geo}/Pressure"] = np.array([100, 30, 10, 3, 1], dtype=np.float32)
        # TAI seconds since 1993-01-01 for 2026-01-01 00:00-02:00 (no leap seconds)
        start = (np.datetime64("2026-01-01") - np.datetime64("1993-01-01")) \
            / np.timedelta64(1, "s")
        f[f"{geo}/Time"] = start + np.array([0.0, 3600, 7200])


def test_q8_mls_newest_version_used_and_recorded(tmp_path):
    # Q8 (fixed 2026-11-25, MATLAB and Python): with v05-04 and v06-04 files for
    # the same day, the first in name order (v05-04) was used; now the newest
    (tmp_path / "2026").mkdir()
    old = "MLS-Aura_L2GP-Temperature_v05-04-c01_2026d001.he5"
    new = "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d001.he5"
    _write_mls(tmp_path / "2026" / old, 210.0)
    _write_mls(tmp_path / "2026" / new, 230.0)
    with pytest.warns(UserWarning, match="more than one matching MLS file"):
        ds = load_limb_sounders("2026-01-01", "MLS", file_path=tmp_path,
                                original_z=True)
    assert np.all(ds.temperature.values == 230)  # exact: copied, not computed
    assert ds.attrs["files_read"] == new
    assert ds.attrs["multiple_matches"] == f"2026-01-01: {old}, {new} -> used {new}"


def test_attributes_record_the_load(syn):
    import json

    from atmosgw.curate import __version__

    # D1 to D2 + 1 day: no file for 2007-01-03
    ds = misc(syn, "synlin", [D1, "2007-01-03"], lat_range=(-90, 90))
    assert ds.attrs["atmosgw_curate_version"] == __version__
    assert ds.attrs["files_read"] == "synlin_2007d001.mat, synlin_2007d002.mat"
    assert ds.attrs["files_not_found"] == "2007-01-03"
    assert ds.attrs["multiple_matches"] == ""
    options = json.loads(ds.attrs["options"])
    assert options["time_range"] == [D1, "2007-01-03"]
    assert options["instrument"] == "Misc"
    assert options["height_range"] == [0, "inf"]
    assert options["lat_range"] == [-90, 90]
    assert options["original_z"] is False and options["height_scale"] is None
    # per variable: where it came from in the file, and how it was obtained
    assert ds.temperature.attrs["source_variable"] == "Temp"
    assert "estimated from height" in ds.pressure.attrs["comment"]  # no Pres in file
    assert ds.altitude.attrs["comment"].startswith("the height scale")
    assert "comment" not in ds.temperature.attrs


def test_several_misc_files_for_a_day_warn(syn, tmp_path):
    import shutil

    for name in ("dup_a_2007d001.mat", "dup_b_2007d001.mat"):
        shutil.copy(syn / "synlin_2007d001.mat", tmp_path / name)
    with pytest.warns(UserWarning, match="1 day"):
        ds = misc(tmp_path, "dup")
    assert ds.attrs["files_read"] == "dup_a_2007d001.mat"  # the first, as MATLAB
    assert "dup_b_2007d001.mat" in ds.attrs["multiple_matches"]


def test_attributes_written_to_netcdf(syn, tmp_path):
    import xarray as xr

    ds = misc(syn, "synlin")
    ds.to_netcdf(tmp_path / "out.nc")
    with xr.open_dataset(tmp_path / "out.nc") as back:
        assert back.attrs["options"] == ds.attrs["options"]
        assert back.pressure.attrs["comment"] == ds.pressure.attrs["comment"]
