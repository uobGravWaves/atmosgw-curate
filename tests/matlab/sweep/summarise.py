"""
Summarise a sweep (see README.md): counts by status and instrument, and an
explanation for every difference beyond the comparison tolerances. Known
last-bit rounding (found 2026-11-25) is accepted with the limits below;
anything else is listed as unexplained.

Usage: python summarise.py WORK_FOLDER
"""

import collections
import glob
import json
import os
import sys


def explain(rec, name, field):
    """Why a field beyond tolerance is still acceptable, or None."""
    if field.get("nan_diff", 1):
        return None
    if name == "Alt" and rec["inst"] == "MLS" and field["max_rel"] <= 1e-15:
        return "MLS heights from pressure (p2h), last bit"
    if name == "Pres" and field["max_rel"] <= 1e-14:
        return "pressure from height above 86 km (table part), last bits"
    if name == "Time" and field["max_abs"] <= 1.2e-10:
        return "time, one unit in the last place of a datenum (10 us)"
    return None


def main(work):
    paths = sorted(glob.glob(os.path.join(work, "results_batch_*.jsonl")))
    records = [json.loads(line) for path in paths for line in open(path)]
    counts = collections.Counter(r["status"] for r in records)
    print(len(records), "calls:", dict(counts))
    explained, unexplained = collections.Counter(), []
    for rec in records:
        if rec["status"] in ("identical", "within_tol"):
            continue
        if rec["status"] != "DIFF":
            unexplained.append((rec, rec["status"]))
            continue
        bad = {k: v for k, v in rec["fields"].items() if v["status"] == "DIFF"}
        why = {k: explain(rec, k, v) for k, v in bad.items()}
        if all(why.values()):
            for reason in why.values():
                explained[(rec["inst"], reason)] += 1
        else:
            unexplained.append((rec, {k: v for k, v in bad.items() if not why[k]}))
    for (inst, reason), n in sorted(explained.items()):
        print(f"  explained: {inst}: {reason}: {n}")
    totals = collections.defaultdict(lambda: [0, 0, 0.0, 0.0])
    for rec in records:
        t = totals[(rec["inst"], rec["kind"])]
        t[0] += 1
        t[1] += rec.get("n_profiles", [0])[0]
        t[2] += rec.get("ml_seconds", 0)
        t[3] += rec.get("py_seconds", 0)
    for (inst, kind), (n, prof, ml, py) in sorted(totals.items()):
        print(f"  {inst} {kind}: {n} calls, {prof} profiles, "
              f"MATLAB {ml:.0f} s, Python {py:.0f} s")
    print("UNEXPLAINED:", len(unexplained))
    for rec, detail in unexplained:
        print("   ", rec["id"], rec.get("n_profiles"), detail,
              rec.get("ml_error", "")[:200], rec.get("py_error", "")[:200])
    return len(unexplained)


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1]) else 0)
