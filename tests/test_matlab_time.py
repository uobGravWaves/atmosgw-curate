"""Tests for atmosgw.curate.matlab_time."""

import datetime as dt

import numpy as np
import pytest

from atmosgw.curate.matlab_time import (
    datenum,
    datenum_to_datetime,
    datetime_to_datenum,
    make_datetime,
)


# values quoted in MathWorks documentation / MATLAB Answers
def test_documented_values():
    assert datenum(2001, 12, 19) == 731204
    assert datenum(1994, 12, 19) == 728647
    assert datenum(1963, 0, 1) == datenum(1963, 1, 1) == 716972  # month 0 -> January
    assert datenum(1970, 1, 1) == 719529


def test_documented_carryover():
    # datestr([2009 22 03 00 00 00]) -> 03-Oct-2010
    assert make_datetime(2009, 22, 3) == np.datetime64("2010-10-03")
    # datestr([2010 07 -05 00 00 00]) -> 25-Jun-2010
    assert make_datetime(2010, 7, -5) == np.datetime64("2010-06-25")
    # minutes 70 -> hours + 1, minutes 10
    assert make_datetime(2000, 1, 1, 0, 70) == np.datetime64("2000-01-01T01:10")
    # minutes -15 -> previous hour, minute 45
    assert make_datetime(2000, 1, 1, 5, -15) == np.datetime64("2000-01-01T04:45")


def test_users_example():
    assert make_datetime(2020, 1, 35) == np.datetime64("2020-02-04")
    assert datenum(2020, 1, 35) == datenum(2020, 2, 4)


def test_against_python_calendar():
    # Python's ordinal: 0001-01-01 is 1; MATLAB datenum of that date is 367
    rng = np.random.default_rng(1)
    for _ in range(2000):
        d = dt.date.fromordinal(int(rng.integers(1, dt.date(9999, 12, 31).toordinal())))
        assert datenum(d.year, d.month, d.day) == d.toordinal() + 366


def test_leap_years_and_year_zero():
    assert datenum(0, 1, 1) == 1
    assert datenum(0, 3, 1) - datenum(0, 2, 28) == 2  # year 0 is a leap year
    assert datenum(1900, 3, 1) - datenum(1900, 2, 28) == 1  # 1900 is not
    assert datenum(2000, 3, 1) - datenum(2000, 2, 28) == 2  # 2000 is


def test_arrays_broadcast():
    out = datenum(2020, 1, np.arange(1, 4))
    assert np.array_equal(out, datenum(2020, 1, 1) + np.arange(3))


def test_fractional_days():
    assert datenum(2010, 1, 1.5) == datenum(2010, 1, 1, 12)


def test_round_trip_and_nan():
    t = np.array(["2007-01-01T06:30:15.25", "NaT"], dtype="datetime64[us]")
    dn = datetime_to_datenum(t)
    assert np.isnan(dn[1])
    back = datenum_to_datetime(dn)
    assert back[0] == t[0] and np.isnat(back[1])


def test_non_integer_month_rejected():
    with pytest.raises(ValueError):
        datenum(2020, 1.5, 1)


def test_sentinel_dates_from_get_limbsounders():
    # get_limbsounders uses datenum(9999, 999, 999) as "still running"
    # month 999 -> year + 83, month 3; then day 999 of March 10082
    assert datenum(9999, 999, 999) == datenum(10082, 3, 1) + 998
    assert make_datetime(9999, 999, 999) == np.datetime64("10084-11-23")
    assert make_datetime(2004, 1, 35) == np.datetime64("2004-02-04")
