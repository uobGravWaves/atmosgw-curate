"""
MATLAB-compatible date helpers.

Convert between MATLAB serial date numbers ("datenums") and NumPy
``datetime64``, and reproduce MATLAB's ``datenum(Y, M, D, H, MI, S)``,
including its carry-over of out-of-range values (e.g. ``datenum(2020, 1, 35)``
is 4 February 2020). Used wherever ported MATLAB code or MATLAB-written files
give times as datenums; inside ``atmosgw.curate`` all times are
``datetime64[us]``.

Inputs and outputs
------------------
- ``datenum(y, m=1, d=1, h=0, mi=0, s=0)``: numbers or arrays (broadcast
  together, like MATLAB) -> float, or float64 array.
- ``make_datetime(y, m=1, d=1, h=0, mi=0, s=0)``: as ``datenum`` ->
  ``datetime64[us]``.
- ``datenum_to_datetime(dn, resolution="ms")``: datenum(s) ->
  ``datetime64[us]``, rounded to ``"ms"`` (default) or ``"us"``.
- ``datetime_to_datenum(t)``: anything NumPy can turn into ``datetime64``
  (ISO strings, Python datetimes) -> float, or float64 array.

Missing values: NaN in -> NaT out, and NaT in -> NaN out.

MATLAB rules reproduced (MathWorks documentation, "Carryover in Date Vectors
and Strings", and the ``datenum`` reference page):

- values outside the normal range carry into the next unit up
  (month 22 of 2009 -> October 2010; minute 70 -> +1 hour 10 minutes);
- days, hours, minutes and seconds may be zero, negative or fractional
  (day 0 is the last day of the previous month; day -5 of July -> 25 June);
- month values less than 1 are set to 1 (so month 0 is January, not December).

Non-integer months and years are rejected: MATLAB's behaviour for them is
undocumented.

Datenum 1 is 1 January year 0 (proleptic ISO/Gregorian calendar, which has a
year 0), so datenum 719529 is 1970-01-01T00:00.

Outputs are ``datetime64[us]``, not ``[ns]``: ``[ns]`` only covers years
1678-2261, and MATLAB code uses far-future sentinels such as
``datenum(9999, 999, 999)``.

Precision: a float64
datenum near the present resolves only about 10 microseconds, so
``datenum_to_datetime`` rounds to the nearest millisecond by default.
``make_datetime`` does not go through a float datenum and keeps microsecond
precision.

Origin
------
New code written to match documented MATLAB behaviour; not a port of MATLAB
source.

Validation
----------
``tests/test_matlab_time.py``: MathWorks' documented examples (values and
carry-over cases), 2,000 random dates against Python's own calendar, leap-year
rules including year 0, broadcasting, NaN/NaT handling, and the sentinel dates
used in ``get_limbsounders``. Not yet compared against MATLAB itself.

Status
------
Preliminary. AI-originated (Claude), reviewed by Corwin Wright.

Usage
-----
>>> from atmosgw.curate.matlab_time import datenum, make_datetime
>>> datenum(2020, 1, 35) == datenum(2020, 2, 4)
True
>>> make_datetime(2020, 1, 35)
np.datetime64('2020-02-04T00:00:00.000000')

Change log
----------
2026/11/23 [Claude] initial version

"""

import numpy as np

DATENUM_UNIX_EPOCH = 719529  # datenum(1970, 1, 1)

# days before the start of each month, for normal and leap years
_DAYS_BEFORE_MONTH = np.array(
    [
        [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
        [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335],
    ]
)


def datenum(year, month=1, day=1, hour=0, minute=0, second=0):
    """
    MATLAB-style serial date number from date parts, with MATLAB's carry-over.

    Inputs can be scalars or arrays (broadcast together, like MATLAB). Returns
    a float for scalar input, otherwise a float64 array.

    Examples
    --------
    >>> datenum(2020, 1, 35) == datenum(2020, 2, 4)
    True
    >>> datenum(2001, 12, 19)
    731204.0
    """
    first_of_month, day, hour, minute, second = _first_of_month(
        year, month, day, hour, minute, second
    )
    result = first_of_month + (day - 1) + hour / 24 + minute / 1440 + second / 86400
    return float(result) if result.ndim == 0 else result


def make_datetime(year, month=1, day=1, hour=0, minute=0, second=0):
    """
    Like datenum(), with the same carry-over rules, but returns datetime64[us].
    NaN inputs give NaT.

    Example: make_datetime(2020, 1, 35) -> 2020-02-04T00:00 (as datetime64[us])
    """
    first_of_month, day, hour, minute, second = _first_of_month(
        year, month, day, hour, minute, second
    )
    # seconds since 1970-01-01, kept in seconds (not days) to preserve precision
    seconds = (
        (first_of_month - DATENUM_UNIX_EPOCH + day - 1) * 86400
        + hour * 3600 + minute * 60 + second
    )
    return _seconds_to_datetime(seconds, resolution_us=1)


def datenum_to_datetime(dn, resolution="ms"):
    """
    Convert MATLAB datenum(s) to numpy datetime64[us]. NaN becomes NaT.

    resolution: "ms" (default) or "us". A float64 datenum only resolves about
    10 microseconds today, so finer digits are rounding noise.
    """
    resolution_us = {"ms": 1000, "us": 1}[resolution]
    dn = np.asarray(dn, dtype=float)
    return _seconds_to_datetime((dn - DATENUM_UNIX_EPOCH) * 86400, resolution_us)


def datetime_to_datenum(t):
    """
    Convert datetime64 values (or anything numpy can turn into datetime64,
    e.g. ISO strings or Python datetimes) to MATLAB datenums. NaT becomes NaN.
    """
    t = np.asarray(t, dtype="datetime64[us]")
    nat = np.isnat(t)
    days = t.astype("int64") / 86400e6 + DATENUM_UNIX_EPOCH
    out = np.where(nat, np.nan, days)
    return float(out) if out.ndim == 0 else out


# supporting functions


def _first_of_month(year, month, day, hour, minute, second):
    """Apply MATLAB's month rules; return datenum of the 1st plus remaining parts."""
    inputs = (year, month, day, hour, minute, second)
    parts = [np.asarray(x, dtype=float) for x in inputs]
    year, month, day, hour, minute, second = np.broadcast_arrays(*parts)

    if np.any(np.isfinite(month) & (month != np.floor(month))):
        raise ValueError(
            "Non-integer months are not supported (MATLAB's behaviour is undocumented)."
        )
    if np.any(np.isfinite(year) & (year != np.floor(year))):
        raise ValueError("Non-integer years are not supported.")

    # MATLAB rule: months below 1 are set to 1; months above 12 carry into years
    month = np.where(month < 1, 1, month)
    year = year + np.floor((month - 1) / 12)
    month = (month - 1) % 12 + 1

    # datenum of 1 January (formula also valid for year 0 and negative years)
    year_start = (
        365 * year + np.ceil(year / 4) - np.ceil(year / 100) + np.ceil(year / 400) + 1
    )
    leap = ((year % 4 == 0) & (year % 100 != 0)) | (year % 400 == 0)
    month_index = np.nan_to_num(month, nan=1).astype(int) - 1
    first_of_month = year_start + _DAYS_BEFORE_MONTH[leap.astype(int), month_index]
    first_of_month = np.where(np.isfinite(month), first_of_month, np.nan)
    return first_of_month, day, hour, minute, second


def _seconds_to_datetime(seconds, resolution_us):
    """Seconds since 1970 -> datetime64[us], rounded to resolution_us microseconds."""
    microseconds = np.round(np.asarray(seconds) * 1e6 / resolution_us) * resolution_us
    nat = ~np.isfinite(microseconds)
    out = np.where(nat, 0, microseconds).astype(np.int64)
    out = out.astype("datetime64[us]")
    out = np.where(nat, np.datetime64("NaT", "us"), out)
    return out[()] if out.ndim == 0 else out
