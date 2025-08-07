#!/usr/bin/env python3
"""
Benchmark SCONE (Fortran Monte-Carlo solver) across a matrix of
  ▸ geometry definitions,
  ▸ acceleration schemes, and
  ▸ source-population sizes (‘pop’).

For every (geometry → accelerationMethod → pop) triple the script

  1. rewrites the appropriate input deck
  2. runs   ./Build/scone.out <deck>
  3. parses the solver’s stdout for timing blocks
  4. appends the result to an in-memory table

When all experiments finish the table is flushed to
  scone_benchmark_<YYYYMMDD_HHMMSS>.csv   in the current directory.

The input files are **always restored** to their pristine state, even
after an interruption (SIGINT / Ctrl-C).
"""

from __future__ import annotations

import csv
import datetime as _dt
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# ─────────────────────────── USER-EDITABLE CONTROL PANEL ───────────────────────
# NB  Changing these lists is the *only* thing most users ever need to do.

POP_VALUES: List[int] = [100, 200]                        # innermost
ACCELERATION_METHODS: List[str] = ["patchSingle", "octree"]         # middle
GEOMETRY_CASES: List[str] = [                                       # outermost
    # polyhedral (HEX or POLY) meshes
    "FinalFuelPinHex72",
    "FinalFuelPinHex243",
    # "FinalFuelPinHex576",
    # "FinalFuelPinHex1125",
    # "FinalFuelPinHex1944",
    # "FinalFuelPinHex3087",
    # "FinalFuelPinHex4608",
    # "FinalFuelPinPoly264",
    # "FinalFuelPinPoly436",
    # "FinalFuelPinPoly468",
    # "FinalFuelPinPoly940",
    # "FinalFuelPinPoly1560",
    # tetrahedral meshes
    #"FinalFuelPinTet137",
    # "FinalFuelPinTet298",
     "FinalFuelPinTet427",
    # "FinalFuelPinTet660",
]
# ───────────────────────── END OF USER-EDITABLE SECTION ────────────────────────


# Project layout (adapt if you move directories)
BASE_DIR = Path(__file__).resolve().parent
EXECUTABLE = BASE_DIR / "Build" / "scone.out"
INPUT_DIR = BASE_DIR / "InputFiles"
TET_DECK = INPUT_DIR / "SCONE_ToyProblemTet"
POLY_DECK = INPUT_DIR / "SCONE_ToyProblemPoly"

# Regexes for parsing SCONE stdout
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


# Helpers ───────────────────────────────────────────────────────────────────────
def _load_original_decks() -> Dict[Path, List[str]]:
    """Read the untouched decks once at start-up."""
    decks: Dict[Path, List[str]] = {}
    for p in (TET_DECK, POLY_DECK):
        decks[p] = p.read_text().splitlines(keepends=True)
        # Make an on-disk backup the first time we run, “just in case”.
        bak = p.with_suffix(".bak")
        if not bak.exists():
            shutil.copy2(p, bak)
    return decks


def _prepare_deck(
    deck_lines: List[str],
    pop: int,
    geometry: str,
    accel: str,
) -> List[str]:
    """Return *new* list of lines with the required substitutions."""
    out = deck_lines.copy()

    # Replace population line (assumed unique token 'pop')
    for i, l in enumerate(out):
        if l.strip().startswith("pop"):
            out[i] = re.sub(r"pop\s+\d+;", f"pop      {pop};", l)
            break
    else:
        raise RuntimeError("No 'pop' line found in deck.")

    # Replace geometry identifier + acceleration method in *same* line
    for i, l in enumerate(out):
        if "FinalFuelPin" in l and ("patchSingle" in l or "octree" in l):
            # geometry
            l = re.sub(
                r"FinalFuelPin(?:Tet|Hex|Poly)\d+",
                geometry,
                l,
            )
            # accel method
            l = re.sub(r"\b(patchSingle|octree)\b", accel, l)
            out[i] = l
            break
    else:
        raise RuntimeError("Could not locate geometry/acceleration line.")

    return out


def _run_scone(deck_path: Path) -> str:
    """Run SCONE once and return full stdout (stderr is suppressed)."""
    proc = subprocess.run(
        [str(EXECUTABLE), str(deck_path.resolve())],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        cwd=BASE_DIR,
        text=True,          # decode as UTF-8
        check=True,         # raise upon non-zero exit
    )
    return proc.stdout


def _extract_timings(stdout: str) -> Tuple[float, str, float, str, float, str]:
    """Return (initCPU, initWall, cyc1CPU, cyc1Wall, cyc2CPU, cyc2Wall)."""
    m_init = _INIT_RE.search(stdout)
    cycles = _IN_CYCLE_RE.findall(stdout)

    if not m_init or len(cycles) < 2:
        raise ValueError("Unexpected SCONE output – timing blocks missing.")

    init_cpu, init_wall = m_init.groups()
    (c1_cpu, c1_wall), (c2_cpu, c2_wall) = cycles[:2]

    return (
        float(init_cpu),
        init_wall,
        float(c1_cpu),
        c1_wall,
        float(c2_cpu),
        c2_wall,
    )


def _select_template(geometry: str) -> Path:
    """Heuristic: any ‘Tet…’ → tetrahedral deck, everything else → poly deck."""
    return TET_DECK if "Tet" in geometry else POLY_DECK


# Main driver ──────────────────────────────────────────────────────────────────
def main() -> None:
    if not EXECUTABLE.exists():
        sys.exit(f"Fatal: executable not found at {EXECUTABLE}")

    originals = _load_original_decks()
    results = []

    try:
        for geom in GEOMETRY_CASES:                       # outer-most
            deck_template = _select_template(geom)
            for accel in ACCELERATION_METHODS:            # intermediate
                for pop in POP_VALUES:                    # inner-most
                    print(f"→ {geom:>20s} | {accel:<11s} | pop={pop:4d}",
                          end=" … ", flush=True)

                    # 1. materialise temporary deck on disk
                    modified = _prepare_deck(
                        originals[deck_template], pop, geom, accel
                    )
                    deck_template.write_text("".join(modified))

                    # 2. run solver & harvest timings
                    stdout = _run_scone(deck_template)
                    timings = _extract_timings(stdout)

                    # 3. store statistics
                    results.append(
                        {
                            "geometry": geom,
                            "acceleration": accel,
                            "pop": pop,
                            "init_cpu_s": timings[0],
                            "init_wall_hms": timings[1],
                            "cycle1_cpu_s": timings[2],
                            "cycle1_wall_hms": timings[3],
                            "cycle2_cpu_s": timings[4],
                            "cycle2_wall_hms": timings[5],
                        }
                    )
                    print("✓")

    finally:
        # Always restore pristine decks
        for p, lines in originals.items():
            p.write_text("".join(lines))

    # ─────────────────────── report to CSV ─────────────────────────
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = BASE_DIR / f"scone_benchmark_{stamp}.csv"
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)

    print(f"\nAll experiments complete. Results → {out_path}")


if __name__ == "__main__":
    main()
