"""
Make the cases of the MATLAB-vs-Python sweep (see README.md).

Dates are sampled (fixed seed) from the days that have files, per instrument
and year. Writes batch_NNN.json files (lists of calls) into the work folder.

Usage: python make_cases.py WORK_FOLDER
"""

import datetime as dt
import json
import os
import random
import re
import sys

import h5py
import numpy as np

from atmosgw.curate.matlab_time import datetime_to_datenum
from atmosgw.curate.obs._limb_alomar import DATA_VERSION
from atmosgw.curate.paths import dataset_dir

SEED = 20261125
BATCH_SIZE = 30
RAW_PER_YEAR = {"ACE": 4, "GNSS": 4, "HIRDLS": 12, "MIPAS": 4, "MLS": 4, "SABER": 4,
                "SOFIE": 4, "ALOMAR": 4}
PWS = ("HIRDLS", "SABER", "MLS", "GNSS")
MONTHS = ("January February March April May June July August September October "
          "November December").split()
# file-name patterns giving the date of each file (raw data)
PATTERNS = {"ACE": r"ace_v52_(\d{4})d(\d{3})\.mat$",
            "GNSS": r"merged_ro_(\d{4})_(\d{3})\.mat$",
            "HIRDLS": r"_(\d{4})d(\d{3})", "MLS": r"_(\d{4})d(\d{3})\.he5$",
            "SOFIE": r"Level2_(\d{4})(\d{3})", "MIPAS": r"_(\d{4})(\d{2})(\d{2})_nom",
            "SABER": r"_(" + "|".join(MONTHS) + r")(\d{4})_"}


def main(work):
    rng = random.Random(SEED)
    cases = []
    for inst, per_year in RAW_PER_YEAR.items():
        days = days_with_files(inst, "raw")
        for year in sorted({d.year for d in days}):
            pool = [d for d in days if d.year == year]
            picks = rng.sample(pool, min(per_year, len(pool)))
            if inst == "SABER":  # a first of the month (reads the day before)
                firsts = [d for d in pool if d.day == 1]
                if firsts:
                    picks[0] = rng.choice(firsts)
            picks = sorted(picks)
            for i, day in enumerate(picks):
                n_days = 3 if i == len(picks) - 1 else 1  # one 3-day range a year
                for config in ("interp", "origz"):
                    cases.append(make_case(inst, "raw", config, day, n_days))
    for inst in PWS:
        days = days_with_files(inst, "pws")
        for year in sorted({d.year for d in days}):
            pool = [d for d in days if d.year == year]
            for i, day in enumerate(sorted(rng.sample(pool, min(2, len(pool))))):
                cases.append(make_case(inst, "pws", "interp", day, 2 if i else 1))

    batches, current = [], []
    for case in cases:
        if current and ((case["inst"], case["kind"]) != (current[0]["inst"],
                                                         current[0]["kind"])
                        or len(current) >= BATCH_SIZE):
            batches.append(current)
            current = []
        current.append(case)
    batches.append(current)
    for i, batch in enumerate(batches):
        with open(os.path.join(work, f"batch_{i:03d}.json"), "w") as f:
            json.dump(batch, f)
    print(f"{len(cases)} calls in {len(batches)} batches")


def days_with_files(name, kind):
    """Sorted dates (datetime.date) for which the dataset has a file."""
    root = str(dataset_dir(name, kind))
    days = set()
    if name == "ALOMAR":  # yearly files, one group per measurement day
        for file in sorted(os.listdir(root)):
            if file.endswith(f"_{DATA_VERSION}_file4juli.h5"):
                with h5py.File(os.path.join(root, file)) as h:
                    days |= {dt.date.fromisoformat(k) for k in h.keys()}
        return sorted(days)
    pattern = r"_PWs_(\d{4})d(\d{3})\.mat$" if kind == "pws" else PATTERNS[name]
    for folder, _, files in os.walk(root):
        for file in files:
            found = re.search(pattern, file, re.IGNORECASE)
            if not found:
                continue
            parts = found.groups()
            if name == "SABER" and kind == "raw":  # monthly files: every day
                month = MONTHS.index(parts[0].capitalize()) + 1
                day = dt.date(int(parts[1]), month, 1)
                while day.month == month:
                    days.add(day)
                    day += dt.timedelta(days=1)
            elif name == "MIPAS" and kind == "raw":
                days.add(dt.date(int(parts[0]), int(parts[1]), int(parts[2])))
            else:
                days.add(dt.date(int(parts[0]), 1, 1)
                         + dt.timedelta(days=int(parts[1]) - 1))
    return sorted(days)


def make_case(inst, kind, config, first, n_days):
    """One call: ISO dates for Python, datenums for MATLAB, and the options."""
    dates = [first + dt.timedelta(days=k) for k in (0, n_days - 1)][: 1 + (n_days > 1)]
    options = {"file_source": True}
    if config == "interp":
        options["dx"] = True
    else:
        options["original_z"] = True
    if kind == "pws":
        options["hindley_pws"] = True
    return {"id": f"{inst}_{kind}_{config}_{first.isoformat()}_{n_days}d",
            "inst": inst, "kind": kind, "config": config,
            "time_range": [d.isoformat() for d in dates],
            "datenums": [float(datetime_to_datenum(np.datetime64(d.isoformat())))
                         for d in dates],
            "options": options}


if __name__ == "__main__":
    main(sys.argv[1])
