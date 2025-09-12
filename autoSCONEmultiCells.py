# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Run SCONE over GEOM_CASES and record ONLY one CSV:
    cell_counts_<timestamp>.csv with columns:
        geometry, cell_counts
where cell_counts is the exact integer list from the line immediately
after the marker:
    "Starting the procedure to calculate number of cells for each type"

Notes:
  - We still edit decks for POP / geometry / acceleration as before.
  - We restore original decks on exit, even after errors.
"""

from __future__ import annotations
import csv
import datetime as dt
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

# ───────────── USER-EDITABLE CONTROL PANEL ─────────────
POP_VALUES  = [1000]
ACCEL       = ["patchSingle"]
GEOM_CASES  = [
    # polyhedral
    "FinalFuelPinHex72",  "FinalFuelPinHex243", "FinalFuelPinHex576",
    "FinalFuelPinHex1125","FinalFuelPinHex1944","FinalFuelPinHex3087",
    "FinalFuelPinPoly264","FinalFuelPinPoly436",
    "FinalFuelPinPoly468","FinalFuelPinPoly940",
    # tetrahedral
    "FinalFuelPinTet298", "FinalFuelPinTet427", "FinalFuelPinTet660",
]
# ───────────────────────────────────────────────────────

BASE      = Path(__file__).resolve().parent
EXEC      = BASE / "Build" / "scone.out"
IN_DIR    = BASE / "InputFiles"
TET_DECK  = IN_DIR / "SCONE_ToyProblemTet"
POLY_DECK = IN_DIR / "SCONE_ToyProblemPoly"

# marker and number parsing
COUNTS_MARKER_RE = re.compile(
    r"Starting the procedure to calculate number of cells for each type",
    re.I,
)
INT_RE = re.compile(r"-?\d+")

KILLED_CODES = {9, -9, 137, -137}  # SIGKILL variants

# ───────────────────── helpers: decks ─────────────────────────
def _load_decks() -> Dict[Path, List[str]]:
    decks: Dict[Path, List[str]] = {}
    for p in (TET_DECK, POLY_DECK):
        decks[p] = p.read_text().splitlines(keepends=True)
        bak = p.with_suffix(".bak")
        if not bak.exists():
            shutil.copy2(p, bak)
    return decks

def _edit(lines: List[str], pop: int, geom: str, accel: str) -> List[str]:
    """Edit a deck: pop value, geometry case name, and acceleration keyword."""
    out = lines.copy()
    # pop
    for i, l in enumerate(out):
        if l.strip().startswith("pop"):
            out[i] = re.sub(r"pop\s+\d+;", f"pop      {pop};", l)
            break
    # geometry + acceleration (first line that mentions FinalFuelPin and an accel tag)
    for i, l in enumerate(out):
        if "FinalFuelPin" in l and ('patchSingle' in l or 'octree' in l):
            l = re.sub(r"FinalFuelPin(?:Tet|Hex|Poly)\d+", geom, l)
            l = re.sub(r"\b(patchSingle|octree)\b", accel, l)
            out[i] = l
            break
    return out

def _tmpl(geom: str) -> Path:
    return TET_DECK if "Tet" in geom else POLY_DECK

# ───────────────────── helpers: run & parse ─────────────────
def _run(deck: Path) -> str:
    cmd = [str(EXEC), str(deck.resolve()), "--omp", "1"]
    try:
        res = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=BASE,
            check=True,
        )
        return res.stdout
    except subprocess.CalledProcessError as e:
        # Even if SCONE ends with ERROR STOP, we still get stdout for parsing.
        if e.returncode in KILLED_CODES or "Killed" in (e.stdout or ""):
            return e.stdout or ""  # still return what we have
        return e.stdout or ""

def _extract_cell_counts(txt: str) -> Optional[str]:
    """
    Find the marker, then take the first subsequent non-empty line
    that contains at least one integer. Join all integers with spaces.
    """
    if not txt:
        return None
    m = COUNTS_MARKER_RE.search(txt)
    if not m:
        return None
    tail = txt[m.end():].splitlines()
    # scan a few lines to be robust to blank lines or banners
    for line in tail[:6]:
        if not line.strip():
            continue
        nums = INT_RE.findall(line)
        if nums:
            return " ".join(nums)
        # stop early if we hit another banner line of non-numeric glyphs
        if "><" in line or "<>" in line:
            break
    return None

# ───────────────────── main ────────────────────────────
def main() -> None:
    if not EXEC.exists():
        sys.exit(f"Executable missing: {EXEC}")

    originals = _load_decks()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    counts_csv = BASE / f"cell_counts_{stamp}.csv"

    with counts_csv.open("w", newline="") as fh_c:
        writer = csv.DictWriter(fh_c, fieldnames=["geometry", "cell_counts"])
        writer.writeheader()
        fh_c.flush()

        try:
            for geom in GEOM_CASES:
                for acc in ACCEL:
                    for pop in POP_VALUES:
                        deck = _tmpl(geom)
                        print(f"→ {geom} | {acc} | pop={pop} … ", end="", flush=True)

                        # edit deck and run
                        deck.write_text(''.join(_edit(originals[deck], pop, geom, acc)))
                        out = _run(deck)

                        # always attempt to parse counts (even if SCONE ERROR STOPs)
                        counts = _extract_cell_counts(out)
                        if counts:
                            writer.writerow({"geometry": geom, "cell_counts": counts})
                            fh_c.flush()
                            print("recorded ✓")
                        else:
                            print("no counts found")

        finally:
            # restore pristine input decks
            for p, lines in originals.items():
                p.write_text(''.join(lines))

    print("\nAll runs complete.")
    print(f"Cell counts CSV → {counts_csv}")

if __name__ == "__main__":
    main()
