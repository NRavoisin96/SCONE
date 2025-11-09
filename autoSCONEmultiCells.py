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
  - Robust global edits for POP / geometry / acceleration (no silent misses).
  - Restores original decks on exit, even after errors.
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
ACCEL       = ["patchMulti"]
GEOM_CASES  = [
    # polyhedral
    "FinalFuelPinHex72",  "FinalFuelPinHex243", "FinalFuelPinHex576",
    "FinalFuelPinHex1125","FinalFuelPinHex1944","FinalFuelPinHex3087",
    "FinalFuelPinHex4608","FinalFuelPinHex6561","FinalFuelPinHex9000",
    "FinalFuelPinPoly264","FinalFuelPinPoly436",
    "FinalFuelPinPoly468","FinalFuelPinPoly940",
    "FinalFuelPinPoly1560","FinalFuelPinPoly5810",
    # tetrahedral
    "FinalFuelPinTet137", "FinalFuelPinTet298",
    "FinalFuelPinTet427", "FinalFuelPinTet660", "FinalFuelPinTet1331",
    "FinalFuelPinTet1820", "FinalFuelPinTet2856",
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

# robust edit regexes
GEOM_CASE_RE = re.compile(r"\bFinalFuelPin(?:Tet|Hex|Poly)\d+\b")
ACCEL_RE     = re.compile(r"\b(patchSingle|patchMulti|octree)\b", re.I)
POP_LINE_RE  = re.compile(r"^\s*pop\s+\d+\s*;", re.I)

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
    """
    Global, robust edits:
      - pop: first 'pop <int>;' at line start
      - geometry: replace all FinalFuelPin(Tet|Hex|Poly)<digits>
      - accel: replace any of {patchSingle, patchMulti, octree} with requested accel
    """
    out: List[str] = []
    seen_pop = seen_geom = seen_accel = False

    for l in lines:
        if not seen_pop and POP_LINE_RE.search(l):
            l = POP_LINE_RE.sub(f"pop      {pop};", l)
            seen_pop = True

        if GEOM_CASE_RE.search(l):
            l = GEOM_CASE_RE.sub(geom, l)
            seen_geom = True

        if ACCEL_RE.search(l):
            l = ACCEL_RE.sub(accel, l)
            seen_accel = True

        out.append(l)

    if not seen_geom:
        raise RuntimeError(f"Could not set geometry to '{geom}': token not found in deck.")
    if not seen_pop:
        raise RuntimeError("Could not set 'pop' in deck.")
    # accel may be legitimately absent; keep as soft requirement
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
        if e.returncode in KILLED_CODES or "Killed" in (e.stdout or ""):
            return e.stdout or ""
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
    for line in tail[:20]:  # widened scan window for banners
        if not line.strip():
            continue
        nums = INT_RE.findall(line)
        if nums:
            return " ".join(nums)
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

                        edited = _edit(originals[deck], pop, geom, acc)
                        deck.write_text(''.join(edited))
                        out = _run(deck)

                        counts = _extract_cell_counts(out)
                        if counts:
                            writer.writerow({"geometry": geom, "cell_counts": counts})
                            fh_c.flush()
                            print("recorded ✓")
                        else:
                            print("no counts found")

        finally:
            for p, lines in originals.items():
                p.write_text(''.join(lines))

    print("\nAll runs complete.")
    print(f"Cell counts CSV → {counts_csv}")

if __name__ == "__main__":
    main()