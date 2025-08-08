# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Benchmark SCONE across:
    ▸ geometry cases
    ▸ acceleration methods
    ▸ population sizes (‘pop’)

Key features
────────────
• Always calls ./Build/scone.out … --omp 1
• Continues after out-of-memory SIGKILL (exit 9, -9, 137 or -137, or “Killed”)
• Writes ‘NA’ timings plus a note in that situation
• Restores pristine input decks even after Ctrl-C
"""

from __future__ import annotations
import csv, datetime as _dt, re, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List, Tuple

# ───────────── USER-EDITABLE CONTROL PANEL ──────────────
POP_VALUES        = [
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
ACCEL_METHODS     = ["patchSingle", "octree"]
GEOMETRY_CASES    = [
    # polyhedral
    "FinalFuelPinHex72", "FinalFuelPinHex243", "FinalFuelPinHex576",
    "FinalFuelPinHex1125", "FinalFuelPinHex1944", "FinalFuelPinHex3087",
    "FinalFuelPinHex4608", "FinalFuelPinPoly264", "FinalFuelPinPoly436",
    "FinalFuelPinPoly468", "FinalFuelPinPoly940", "FinalFuelPinPoly1560",
    # tetrahedral
    "FinalFuelPinTet137", "FinalFuelPinTet298",
    "FinalFuelPinTet427", "FinalFuelPinTet660",
]
# ─────────────────────────────────────────────────────────

BASE_DIR   = Path(__file__).resolve().parent
EXECUTABLE = BASE_DIR / "Build" / "scone.out"
INPUT_DIR  = BASE_DIR / "InputFiles"
TET_DECK   = INPUT_DIR / "SCONE_ToyProblemTet"
POLY_DECK  = INPUT_DIR / "SCONE_ToyProblemPoly"

_INIT_RE    = re.compile(r"Initialisation procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)", re.S)
_IN_CYCLE_RE = re.compile(r"In-cycle procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)", re.S)
_KILL_CODES  = {9, -9, 137, -137}   # plain SIGKILL or 128+SIGKILL

# ───────────────────────── helpers ───────────────────────
def _load_original_decks() -> Dict[Path, List[str]]:
    decks = {}
    for p in (TET_DECK, POLY_DECK):
        decks[p] = p.read_text().splitlines(keepends=True)
        bak = p.with_suffix(".bak")
        if not bak.exists(): shutil.copy2(p, bak)
    return decks

def _prepare_deck(lines: List[str], pop: int, geom: str, accel: str) -> List[str]:
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

def _template_for(geom:str)->Path: return TET_DECK if 'Tet' in geom else POLY_DECK

def _run_scone(deck:Path)->str:
    cmd=[str(EXECUTABLE), str(deck.resolve()), "--omp","1"]
    try:
        proc=subprocess.run(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            cwd=BASE_DIR, check=True)
        return proc.stdout
    except subprocess.CalledProcessError as e:
        killed = (e.returncode in _KILL_CODES) or ("Killed" in e.stdout)
        if killed: return "OOM-KILLED"
        raise

def _extract(t:str)->Tuple[float,str,float,str,float,str]:
    m0=_INIT_RE.search(t); cycles=_IN_CYCLE_RE.findall(t)
    if not m0 or len(cycles)<2: raise ValueError("Timing blocks missing")
    init_cpu,init_wall=m0.groups()
    (c1_cpu,c1_wall),(c2_cpu,c2_wall)=cycles[:2]
    return float(init_cpu),init_wall,float(c1_cpu),c1_wall,float(c2_cpu),c2_wall

# ───────────────────────── main ──────────────────────────
def main()->None:
    if not EXECUTABLE.exists(): sys.exit(f"Missing {EXECUTABLE}")
    originals=_load_original_decks(); rows=[]
    try:
        for geom in GEOMETRY_CASES:
            deck=_template_for(geom)
            for accel in ACCEL_METHODS:
                for pop in POP_VALUES:
                    print(f"→ {geom:>22s} | {accel:<11s} | pop={pop:<7d} … ", end='', flush=True)
                    deck.write_text(''.join(_prepare_deck(originals[deck], pop, geom, accel)))
                    out=_run_scone(deck)
                    if out=="OOM-KILLED":
                        rows.append(dict(geometry=geom,acceleration=accel,pop=pop,
                                         init_cpu_s="NA",init_wall_hms="NA",
                                         cycle1_cpu_s="NA",cycle1_wall_hms="NA",
                                         cycle2_cpu_s="NA",cycle2_wall_hms="NA",
                                         note="terminated (out-of-memory)"))
                        print("⚠ killed – continuing")
                    else:
                        ic,iw,c1c,c1w,c2c,c2w=_extract(out)
                        rows.append(dict(geometry=geom,acceleration=accel,pop=pop,
                                         init_cpu_s=ic,init_wall_hms=iw,
                                         cycle1_cpu_s=c1c,cycle1_wall_hms=c1w,
                                         cycle2_cpu_s=c2c,cycle2_wall_hms=c2w,
                                         note=""))
                        print("✓")
    finally:
        for p,lines in originals.items(): p.write_text(''.join(lines))

    stamp=_dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    csv=BASE_DIR/f"scone_benchmark_{stamp}.csv"
    with csv.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    print(f"\nAll experiments complete. Results → {csv}")

if __name__=="__main__": main()
