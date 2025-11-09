# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Extract SCONE grid sizes (x, y, z) into one CSV.

• Runs ./Build/scone.out … --omp 1
• Continues after SIGKILL / OOM, recording “NA” values
• Writes each result row to the CSV as soon as the run finishes
• Restores pristine input decks even if interrupted.
"""

from __future__ import annotations
import csv, datetime as dt, re, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List

# ───────────── USER-EDITABLE CONTROL PANEL ─────────────
POP_VALUES  = [1000]
ACCEL       = ["patchSingle"]
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

# regex for grid size lines
_GRID_RE  = re.compile(r"Grid size in x\s*:\s*(\d+).*?Grid size in y\s*:\s*(\d+).*?Grid size in z\s*:\s*(\d+)", re.S)
_KILLS    = {9, -9, 137, -137}   # SIGKILL returncodes

FIELDS = ["geometry", "acceleration", "pop", "grid_x", "grid_y", "grid_z", "note"]

# ───────────────────── helpers ─────────────────────────
def _load_decks() -> Dict[Path, List[str]]:
    decks = {}
    for p in (TET_DECK, POLY_DECK):
        decks[p] = p.read_text().splitlines(keepends=True)
        bak = p.with_suffix(".bak")
        if not bak.exists(): shutil.copy2(p, bak)
    return decks

def _edit(lines:List[str], pop:int, geom:str, accel:str) -> List[str]:
    out = lines.copy()
    # pop
    for i,l in enumerate(out):
        if l.strip().startswith("pop"):
            out[i] = re.sub(r"pop\s+\d+;", f"pop      {pop};", l); break
    # geometry + accel
    for i,l in enumerate(out):
        if "FinalFuelPin" in l:
            l = re.sub(r"FinalFuelPin(?:Tet|Hex|Poly)\d+", geom, l)
            l = re.sub(r"\b(patchSingle|octree)\b", accel, l)
            out[i] = l; break
    return out

def _tmpl(geom:str)->Path: 
    return TET_DECK if "Tet" in geom else POLY_DECK

def _run(deck:Path)->str:
    cmd=[str(EXEC), str(deck.resolve()), "--omp","1"]
    try:
        res=subprocess.run(cmd, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True,
                           cwd=BASE, check=True)
        return res.stdout
    except subprocess.CalledProcessError as e:
        if e.returncode in _KILLS or "Killed" in e.stdout: return "OOM"
        return e.stdout or "ERROR"

def _parse(txt:str):
    m=_GRID_RE.search(txt)
    if not m: raise ValueError("grid size missing")
    return int(m.group(1)), int(m.group(2)), int(m.group(3))

# ───────────────────── main ────────────────────────────
def main()->None:
    if not EXEC.exists(): sys.exit(f"Executable missing: {EXEC}")
    originals=_load_decks()

    csv_path = BASE / f"scone_grid_{dt.datetime.now():%Y%m%d_%H%M%S}.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader(); fh.flush()

        try:
            for geom in GEOM_CASES:
                deck=_tmpl(geom)
                for acc in ACCEL:
                    for pop in POP_VALUES:
                        tag=f"{geom:>22s} | {acc:<11s} | pop={pop:<7d}"
                        print(f"→ {tag} … ", end="", flush=True)

                        deck.write_text(''.join(_edit(originals[deck], pop, geom, acc)))
                        out=_run(deck)

                        if out=="OOM" or out.startswith("ERROR"):
                            row=dict(geometry=geom, acceleration=acc, pop=pop,
                                     grid_x="NA", grid_y="NA", grid_z="NA",
                                     note="killed/failed")
                            print("⚠")
                        else:
                            try:
                                gx,gy,gz=_parse(out)
                                row=dict(geometry=geom, acceleration=acc, pop=pop,
                                         grid_x=gx, grid_y=gy, grid_z=gz, note="")
                                print("✓")
                            except Exception as ex:
                                row=dict(geometry=geom, acceleration=acc, pop=pop,
                                         grid_x="NA", grid_y="NA", grid_z="NA",
                                         note=f"parse error: {ex}")
                                print("✗ parse error")

                        writer.writerow(row); fh.flush()

        finally:
            for p,lines in originals.items(): p.write_text(''.join(lines))

    print(f"\nAll experiments complete. Results → {csv_path}")

if __name__=="__main__":
    main()
