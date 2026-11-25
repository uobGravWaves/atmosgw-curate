"""Tests for atmosgw.curate.data (formats, states, validate).

No MATLAB equivalent: known-answer checks only. A small valid example of each
format must pass ``validate``; each kind of fault must be reported.
"""

import numpy as np
import pytest
import xarray as xr

from atmosgw.curate.data import (
    FORMATS,
    FormatError,
    availability,
    problems,
    states,
    validate,
)

# ---------------------------------------------------------------- states


def test_registry_bits_unique_powers_of_two():
    bits = [s.bit for s in states.STATES]
    assert len(set(bits)) == len(bits)
    assert all(b > 0 and (b & (b - 1)) == 0 for b in bits)
    assert max(bits) <= np.iinfo(states.FLAG_DTYPE).max
    assert len({s.name for s in states.STATES}) == len(bits)


def test_set_has_and_count():
    f = states.new_flags((2, 3))
    states.set_state(f, [[True, False, False], [False, False, True]], "missing")
    states.set_state(f, [[True, True, False], [False, False, False]], "derived")
    assert states.has_state(f, "missing").tolist() == [[True, False, False],
                                                       [False, False, True]]
    assert states.no_value_mask(f).tolist() == [[True, False, False],
                                                [False, False, True]]
    assert states.state_counts(f) == {"missing": 2, "derived": 2}
    with pytest.raises(KeyError):
        states.set_state(f, [[True] * 3] * 2, "no_such_state")


def test_flag_attrs_follow_cf():
    attrs = states.flag_attrs()
    assert len(attrs["flag_masks"]) == len(attrs["flag_meanings"].split())
    assert attrs["flag_masks"].dtype == states.FLAG_DTYPE


def test_availability_text_round_trip():
    mapping = {"satellite_latitude": "not_in_files", "source_file": "not_requested"}
    text = states.format_availability(mapping)
    assert states.parse_availability(text) == mapping
    assert states.parse_availability("") == {}
    with pytest.raises(ValueError):
        states.format_availability({"x": "lost"})


# ---------------------------------------------------------------- validate


def example_swath(n_along=4, n_cross=3, n_level=5):
    """A minimal valid 3-D swath Dataset."""
    p, p3 = ("along_track", "cross_track"), ("along_track", "cross_track", "level")
    temp = np.full((n_along, n_cross, n_level), 250, dtype=np.float32)
    flag = states.new_flags(temp.shape)
    temp[:1, :1, :1] = np.nan
    states.set_state(flag, np.isnan(temp), "missing")
    ds = xr.Dataset(
        {
            "latitude": (p, np.zeros((n_along, n_cross)), {"units": "degrees_north"}),
            "longitude": (p, np.zeros((n_along, n_cross)), {"units": "degrees_east"}),
            "time": (p, np.full((n_along, n_cross), np.datetime64("2007-01-01T00:00",
                                                                  "us"))),
            "altitude": (("level",), np.arange(n_level, dtype=float) * 3,
                         {"units": "km"}),
            "temperature": (p3, temp, {"units": "K", "origin": "retrieved",
                                       "ancillary_variables": "temperature_flag"}),
            "temperature_flag": (p3, flag, states.flag_attrs()),
            "granule": (("along_track",), np.ones(n_along)),
            "granule_day": (("along_track",),
                            np.full(n_along, np.datetime64("2007-01-01", "us"))),
            "segment": (("along_track",), np.ones(n_along)),
            "distance_along": (("along_track",), np.arange(n_along, dtype=float),
                               {"units": "km"}),
            "distance_cross": (("cross_track",), np.arange(n_cross, dtype=float),
                               {"units": "km"}),
        }
    )
    tracked = [n for n, v in FORMATS["swath"]["variables"].items() if v.get("tracked")]
    ds.attrs = {a: "" for a in FORMATS["swath"]["attrs"]}
    ds.attrs.update(format="swath", format_version=1, instrument="AIRS", product="3d",
                    grid="regular", characterisation="AIRS-3D-placeholder/0",
                    availability=states.format_availability(
                        {**{n: "not_requested" for n in tracked},
                         "temperature_uncertainty": "not_provided"}))
    return ds


def test_valid_example_passes():
    validate(example_swath(), "swath")
    validate(example_swath(n_along=0), "swath")  # no data is still valid


@pytest.mark.parametrize("fault, expected", [
    (lambda ds: ds.drop_vars("temperature_flag"), "missing flag variable"),
    (lambda ds: ds.assign(temperature=ds.temperature.astype(np.float64)), "float32"),
    (lambda ds: ds.assign_attrs(format_version=2), "format_version"),
    (lambda ds: ds.assign_attrs(characterisation="AIRS"), "characterisation"),
    (lambda ds: ds.drop_vars("segment"), "'segment'"),
    (lambda ds: ds.assign_attrs(availability=""), "not listed in 'availability'"),
    (lambda ds: ds.transpose("cross_track", "along_track", "level"), "dimensions"),
])
def test_faults_reported(fault, expected):
    found = problems(fault(example_swath()), "swath")
    assert any(expected in p for p in found), found


def test_unexplained_nan_reported():
    ds = example_swath()
    ds["temperature"][1, 1, 1] = np.nan
    assert any("NaN value(s) with no 'no value' state" in p
               for p in problems(ds, "swath"))


def test_latitude_out_of_range_reported():
    ds = example_swath()
    ds["latitude"][0, 0] = 95.0
    with pytest.raises(FormatError, match="outside"):
        validate(ds, "swath")


def test_availability_helper():
    ds = example_swath()
    assert availability(ds, "temperature") == "present"
    assert availability(ds, "satellite_latitude") == "not_requested"
    with pytest.raises(KeyError):
        availability(ds, "wind")


def test_unknown_format():
    with pytest.raises(ValueError):
        validate(example_swath(), "image")


def test_profile_data_may_be_float64():
    # the limb loader keeps its values exactly (float64); swath data stay float32
    from atmosgw.curate.data.convert import to_profiles

    prof = to_profiles(example_swath())
    validate(prof.assign(temperature=prof.temperature.astype(np.float64)), "profile")
    found = problems(prof.assign(temperature=prof.temperature.astype(np.float16)),
                     "profile")
    assert any("float32 or float64" in p for p in found)
