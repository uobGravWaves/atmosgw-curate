"""
Compare one sweep batch: MATLAB outputs against the Python port (see README.md).

For each call: runs ``_load_limb_arrays`` with the same arguments, compares
every field with MATLAB's (names, shapes, NaN pattern, values), appends one
JSON line to RESULTS, and deletes the MATLAB file.

Tolerances as tests/test_limb_matlab.py: bit-identical, except interpolated
pressure 1e-14 relative; pressure computed from height (ALOMAR, PW files),
heights kept, 1e-15 relative; HIRDLS latitude/longitude 1e-12 degrees; dx
1e-12 relative (all instruments).

Usage: python compare_batch.py BATCH.json OUT_FOLDER RESULTS.jsonl
"""

import inspect
import json
import os
import sys
import time
import traceback
import warnings

import numpy as np

from atmosgw.curate.obs.limb import _load_limb_arrays, load_limb_sounders
from atmosgw.curate.readers import read_mat

DEFAULTS = {name: p.default
            for name, p in inspect.signature(load_limb_sounders).parameters.items()
            if p.default is not inspect.Parameter.empty}


def main(batch_file, out_dir, results):
    with open(batch_file) as f:
        cases = json.load(f)
    with open(results, "a") as out:
        for case in cases:
            mat = os.path.join(out_dir, case["id"] + ".mat")
            out.write(json.dumps(compare_case(case, mat)) + "\n")
            out.flush()
            if os.path.exists(mat):
                os.remove(mat)


def compare_case(case, mat):
    """The result record of one call."""
    rec = {k: case[k] for k in ("id", "inst", "kind", "config", "time_range")}
    try:
        matlab = read_mat(mat)
        rec["ml_seconds"] = float(np.asarray(matlab["Seconds"]).squeeze())
        err = np.asarray(matlab["Err"])
        rec["ml_error"] = str(err.squeeze()) if err.size else ""
        start = time.time()
        rec["py_error"] = ""
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                time_range = case["time_range"]
                got, _ = _load_limb_arrays(
                    time_range if len(time_range) == 2 else time_range[0],
                    case["inst"], **{**DEFAULTS, **case["options"]})
            rec["py_warnings"] = sorted({str(w.message)[:120] for w in caught})
        except Exception as error:
            rec["py_error"] = f"{type(error).__name__}: {error}"
        rec["py_seconds"] = time.time() - start
        if rec["ml_error"] or rec["py_error"]:
            both = rec["ml_error"] and rec["py_error"]
            rec["status"] = "BOTH_ERROR" if both else "ERROR_MISMATCH"
            return rec
        want = matlab["D"]
        fields = compare_fields(got, want, tolerances(case))
        n_matlab = int(np.shape(want["Lat"])[0]) if np.size(want["Lat"]) else 0
        rec["n_profiles"] = [int(np.shape(got["Lat"])[0]), n_matlab]
        rec["fields"] = fields
        statuses = {f["status"] for f in fields.values()}
        rec["status"] = ("DIFF" if "DIFF" in statuses else
                         "within_tol" if "within_tol" in statuses else "identical")
    except Exception:
        rec["status"] = "SWEEP_ERROR"
        rec["trace"] = traceback.format_exc()[-800:]
    return rec


def tolerances(case):
    """{field: ("rtol" or "atol", value)}; other fields must be bit-identical."""
    tol = {"dx": ("rtol", 1e-12)}  # haversine trig: last-bit rounding
    if not case["options"].get("original_z"):
        tol["Pres"] = ("rtol", 1e-14)  # interpolated in log(pressure)
    elif case["inst"] == "ALOMAR" or case["kind"] == "pws":
        tol["Pres"] = ("rtol", 1e-15)  # computed from height
    if case["inst"] == "HIRDLS" and case["kind"] == "raw":
        tol.update(Lat=("atol", 1e-12), Lon=("atol", 1e-12))  # azimuth/reckon
    return tol


def compare_fields(got, want, tol):
    """Per field: status identical, within_tol or DIFF, with the details."""
    out = {}
    for name in sorted(set(got) | set(want)):
        if name not in got or name not in want:
            side = "python" if name not in got else "matlab"
            out[name] = {"status": "DIFF", "why": f"missing in {side}"}
            continue
        if name == "OriginalFiles":
            want_files = [str(np.asarray(f).squeeze())
                          for f in np.asarray(want[name], dtype=object).ravel()]
            same = list(got[name]) == want_files
            out[name] = {"status": "identical" if same else "DIFF"}
            continue
        g = np.asarray(got[name], dtype=float)
        w = np.asarray(want[name], dtype=float)
        if g.size == 0 and w.size == 0:
            out[name] = {"status": "identical"}
            continue
        if g.shape != w.shape:
            why = f"shape {g.shape} vs MATLAB {w.shape}"
            out[name] = {"status": "DIFF", "why": why}
            continue
        nan_diff = int(np.sum(np.isnan(g) != np.isnan(w)))
        both = ~np.isnan(g) & ~np.isnan(w)
        diff = np.abs(g[both] - w[both])
        max_abs = float(diff.max()) if diff.size else 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = diff / np.abs(w[both])
        max_rel = float(np.nanmax(rel)) if np.any(np.isfinite(rel)) else 0.0
        rec = {"nan_diff": nan_diff, "max_abs": max_abs, "max_rel": max_rel}
        kind, value = tol.get(name, (None, None))
        if nan_diff == 0 and max_abs == 0:
            rec["status"] = "identical"
        elif nan_diff == 0 and ((kind == "rtol" and max_rel <= value)
                                or (kind == "atol" and max_abs <= value)):
            rec["status"] = "within_tol"
        else:
            rec["status"] = "DIFF"
        out[name] = rec
    return out


if __name__ == "__main__":
    main(*sys.argv[1:4])
