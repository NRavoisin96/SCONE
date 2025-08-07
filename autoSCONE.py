#!/usr/bin/env python3
"""
Automate SCONE (Fortran) runs for several source-population (‘pop’) sizes,
record timing statistics, and collate everything in a CSV.

Author : ChatGPT (OpenAI), 7 Aug 2025
"""

from __future__ import annotations

import csv
import datetime as _dt
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

# ────────────────────────────── USER-TUNEABLE PARAMETERS ───────────────────────
SCONE_DIR = Path("/home/daeyeun/MPhil_project/SCONE")          # project root
EXECUTABLE = SCONE_DIR / "Build" / "scone.out"                 # compiled code
INPUT_FILE = SCONE_DIR / "InputFiles" / "SCONE_ToyProblem"     # main deck
POP_VALUES = [100, 200, 300]                                   # populations
CSV_PREFIX = "scone_run_summary"                               # output stem
# ───────────────────────────────────────────────────────────────────────────────

# Compiled regexes (multiline / dot-all)
_INIT_RE = re.compile(
    r"Initialisation procedure time.*?"
    r"CPU\s+time:\s+([0-9.E+-]+).*?seconds.*?"
    r"Wall time:\s+([0-9:]+)",
    re.S,
)
_IN_CYCLE_RE = re.compile(
    r"In-cycle procedure time.*?"
    r"CPU\s+time:\s+([0-9.E+-]+).*?seconds.*?"
    r"Wall time:\s+([0-9:]+)",
    re.S,
)

# ------------------------------------------------------------------------------
def _set_pop_in_deck(content: List[str], new_pop: int) -> List[str]:
    """Return a *new* list of lines with the pop value replaced."""
    out = content.copy()
    for idx, line in enumerate(out):
        if line.strip().startswith("pop"):
            out[idx] = re.sub(r"pop\s+\d+;", f"pop      {new_pop};", line)
            break
    else:
        raise RuntimeError("No 'pop' line found in input deck.")
    return out


def _run_scone(pop: int) -> Dict[str, str | float | int]:
    """Run SCONE once for a given population size and return timing data."""
    # 1. Rewrite input file
    INPUT_FILE.write_text("".join(_set_pop_in_deck(original_deck, pop)))

    # 2. Launch solver & capture stdout
    proc = subprocess.run(
        [str(EXECUTABLE), str(INPUT_FILE.resolve())],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        cwd=SCONE_DIR,
        check=True,
    )
    out = proc.stdout

    # 3. Parse timing blocks
    init_m = _INIT_RE.search(out)
    cycle_ms = _IN_CYCLE_RE.findall(out)

    if not init_m or len(cycle_ms) < 2:
        raise ValueError(
            f"Unexpected SCONE output format for pop={pop}. "
            "Could not locate all requested timing blocks."
        )

    init_cpu, init_wall = init_m.groups()
    (cyc1_cpu, cyc1_wall), (cyc2_cpu, cyc2_wall) = cycle_ms[:2]

    return {
        "pop": pop,
        "init_cpu_s": float(init_cpu),
        "init_wall_hms": init_wall,
        "cycle1_cpu_s": float(cyc1_cpu),
        "cycle1_wall_hms": cyc1_wall,
        "cycle2_cpu_s": float(cyc2_cpu),
        "cycle2_wall_hms": cyc2_wall,
    }


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    if not EXECUTABLE.exists():
        sys.exit(f"Executable not found at {EXECUTABLE!s}")

    # Preserve a pristine copy of the original deck
    backup_path = INPUT_FILE.with_suffix(".bak")
    if not backup_path.exists():
        shutil.copy2(INPUT_FILE, backup_path)
    original_deck = INPUT_FILE.read_text().splitlines(keepends=True)

    results: List[Dict[str, str | float | int]] = []
    try:
        for pop in POP_VALUES:
            print(f"⇒ Running SCONE with pop = {pop} …", flush=True)
            stats = _run_scone(pop)
            results.append(stats)
            print("   ✔ completed")

    finally:  # Always restore the file, even after Ctrl-C
        INPUT_FILE.write_text("".join(original_deck))

    # Write CSV
    timestamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = SCONE_DIR / f"{CSV_PREFIX}_{timestamp}.csv"
    with csv_path.open("w", newline="") as f:
        fieldnames = list(results[0].keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)

    print(f"\nAll done. Summary saved to {csv_path}")
