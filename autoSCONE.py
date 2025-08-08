# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Benchmark SCONE across (geometry × acceleration × pop).

• Runs ./Build/scone.out … --omp 1
• Continues after SIGKILL / OOM, recording “NA” timings + note
• **Writes each result row to the CSV as soon as the run finishes**
  → almost no Python-side memory usage for bookkeeping
• Restores pristine input decks even if interrupted.
"""

from __future__ import annotations
import csv, datetime as dt, re, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List, Tuple

# ───────────── USER-EDITABLE CONTROL PANEL ─────────────
POP_VALUES  = [
    500,
    1000,
    1500,
    2000,
    2500,
    3000,
    3500,
    4000,
    4500,
    5000,
    5500,
    6000,
    6500,
    7000,
    8000,
    9000,
    10000,
    15000,
    20000,
    40000,
    60000,
    80000, 
    100000,
    150000
]
ACCEL       = ["patchSingle", "octree"]
GEOM_CASES  = [
    # polyhedral
    "FinalFuelPinHex72",  "FinalFuelPinHex243", "FinalFuelPinHex576",
    "FinalFuelPinHex1125","FinalFuelPinHex1944","FinalFuelPinHex3087",
    "FinalFuelPinHex4608","FinalFuelPinPoly264","FinalFuelPinPoly436",
    "FinalFuelPinPoly468","FinalFuelPinPoly940","FinalFuelPinPoly1560",
    # tetrahedral
    "FinalFuelPinTet137", "FinalFuelPinTet298",
    "FinalFuelPinTet427", "FinalFuelPinTet660",
]
# ───────────────────────────────────────────────────────

BASE      = Path(__file__).resolve().parent
EXEC      = BASE / "Build" / "scone.out"
IN_DIR    = BASE / "InputFiles"
TET_DECK  = IN_DIR / "SCONE_ToyProblemTet"
POLY_DECK = IN_DIR / "SCONE_ToyProblemPoly"

_INIT_RE   = re.compile(r"Initialisation procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)", re.S)
_CYCLE_RE  = re.compile(r"In-cycle procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)", re.S)
_KILLS     = {9, -9, 137, -137}   # SIGKILL returncodes

FIELDS = [
    "geometry", "acceleration", "pop",
    "init_cpu_s", "init_wall_hms",
    "cycle1_cpu_s", "cycle1_wall_hms",
    "cycle2_cpu_s", "cycle2_wall_hms",
    "note",
]

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
    for i,l in enumerate(out):
        if l.strip().startswith("pop"):
            out[i] = re.sub(r"pop\s+\d+;", f"pop      {pop};", l); break
    for i,l in enumerate(out):
        if "FinalFuelPin" in l and ('patchSingle' in l or 'octree' in l):
            l = re.sub(r"FinalFuelPin(?:Tet|Hex|Poly)\d+", geom, l)
            l = re.sub(r"\b(patchSingle|octree)\b", accel, l)
            out[i] = l; break
    return out

def _tmpl(geom:str)->Path: return TET_DECK if "Tet" in geom else POLY_DECK

def _run(deck:Path)->str:
    cmd=[str(EXEC), str(deck.resolve()), "--omp","1"]
    try:
        res=subprocess.run(cmd, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True,
                           cwd=BASE, check=True)
        return res.stdout
    except subprocess.CalledProcessError as e:
        if e.returncode in _KILLS or "Killed" in e.stdout: return "OOM"
        raise

def _parse(txt:str)->Tuple[float,str,float,str,float,str]:
    m0=_INIT_RE.search(txt); cyc=_CYCLE_RE.findall(txt)
    if not m0 or len(cyc)<2: raise ValueError("timing blocks missing")
    ic,iw=m0.groups(); (c1c,c1w),(c2c,c2w)=cyc[:2]
    return float(ic),iw,float(c1c),c1w,float(c2c),c2w

# ───────────────────── main ────────────────────────────
def main()->None:
    if not EXEC.exists(): sys.exit(f"Executable missing: {EXEC}")
    originals=_load_decks()

    csv_path = BASE / f"scone_benchmark_{dt.datetime.now():%Y%m%d_%H%M%S}.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        fh.flush()                                      # ensure header on disk

        try:
            for geom in GEOM_CASES:
                deck=_tmpl(geom)
                for acc in ACCEL:
                    for pop in POP_VALUES:
                        tag=f"{geom:>22s} | {acc:<11s} | pop={pop:<7d}"
                        print(f"→ {tag} … ", end="", flush=True)

                        deck.write_text(''.join(_edit(originals[deck], pop, geom, acc)))
                        out=_run(deck)

                        if out=="OOM":
                            row=dict(geometry=geom, acceleration=acc, pop=pop,
                                     init_cpu_s="NA", init_wall_hms="NA",
                                     cycle1_cpu_s="NA", cycle1_wall_hms="NA",
                                     cycle2_cpu_s="NA", cycle2_wall_hms="NA",
                                     note="terminated (out-of-memory)")
                            print("⚠ killed – continuing")
                        else:
                            ic,iw,c1c,c1w,c2c,c2w=_parse(out)
                            row=dict(geometry=geom, acceleration=acc, pop=pop,
                                     init_cpu_s=ic, init_wall_hms=iw,
                                     cycle1_cpu_s=c1c, cycle1_wall_hms=c1w,
                                     cycle2_cpu_s=c2c, cycle2_wall_hms=c2w,
                                     note="")
                            print("✓")

                        writer.writerow(row)
                        fh.flush()                      # write immediately

        finally:
            for p,lines in originals.items(): p.write_text(''.join(lines))

    print(f"\nAll experiments complete. Results → {csv_path}")

if __name__=="__main__":
    main()
