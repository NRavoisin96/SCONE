# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Benchmark SCONE across (geometry × acceleration × pop), writing:
  • timings CSV (unchanged from your working version)
  • grid_metrics CSV: min angle, min edge length, n_vertices/edges/faces, avg edge length
  • vertex_valence CSV: full vertexValence array (as a single string)
  • edge_valence   CSV: full edgeValence array   (as a single string)

Notes:
  - We DO NOT add POP or ACCEL to the new CSVs (per request).
  - We keep your original behaviour: continue on OOM, flush each row immediately,
    and restore pristine decks at the end.
"""

from __future__ import annotations
import csv, datetime as dt, re, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# ───────────── USER-EDITABLE CONTROL PANEL ─────────────
POP_VALUES  = [
    1000
]
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

# --- regexes for timings (unchanged) ---
_INIT_RE   = re.compile(r"Initialisation procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)", re.S)
_CYCLE_RE  = re.compile(r"In-cycle procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)", re.S)
_KILLS     = {9, -9, 137, -137}   # SIGKILL returncodes

# --- new regexes for mesh-quality block & valence arrays ---
# Works whether numbers are fixed or scientific notation
MIN_ANGLE_RE   = re.compile(r"Minimum angle\s*:\s*([0-9.Ee+\-]+)")
MIN_EDGE_RE    = re.compile(r"Minimum edge length\s*:\s*([0-9.Ee+\-]+)")
N_VERT_RE      = re.compile(r"No\. of vertices\s*:\s*(\d+)")
N_EDGES_RE     = re.compile(r"No\. of edges\s*:\s*(\d+)")
N_FACES_RE     = re.compile(r"No\. of faces\s*:\s*(\d+)")
AVG_EDGE_RE    = re.compile(r"average edge length\s*([0-9.Ee+\-]+)")

# We capture the numeric runs that follow the labels
VERTEX_LABEL_RE = re.compile(r"^\s*vertexValence\s*$", re.M)
EDGE_LABEL_RE   = re.compile(r"^\s*edgeValence\s*$",   re.M)
# A helper to grab a line (or a couple of lines) of integers after a label:
INT_LINE_RE     = re.compile(r"(?:^|\s)(-?\d+)(?=\s|$)")

# CSV schemas
TIMING_FIELDS = [
    "geometry", "acceleration", "pop",
    "init_cpu_s", "init_wall_hms",
    "cycle1_cpu_s", "cycle1_wall_hms",
    "cycle2_cpu_s", "cycle2_wall_hms",
    "note",
]
GRID_FIELDS   = ["geometry", "min_angle", "min_edge_length", "n_vertices", "n_edges", "n_faces", "avg_edge_length"]
VVAL_FIELDS   = ["geometry", "vertexValence"]
EVAL_FIELDS   = ["geometry", "edgeValence"]

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
    """Your original editing logic; unchanged."""
    out = lines.copy()
    # pop
    for i, l in enumerate(out):
        if l.strip().startswith("pop"):
            out[i] = re.sub(r"pop\s+\d+;", f"pop      {pop};", l)
            break
    # geometry + acceleration
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
            cmd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True,
            cwd=BASE, check=True
        )
        return res.stdout
    except subprocess.CalledProcessError as e:
        # Return special token for OOM so the loop continues
        if e.returncode in _KILLS or "Killed" in (e.stdout or ""):
            return "OOM"
        # Otherwise, still return the stdout so we can parse grid metrics / valences from failures
        return e.stdout or ""

def _parse_timings(txt: str) -> Tuple[float, str, float, str, float, str]:
    m0 = _INIT_RE.search(txt); cyc = _CYCLE_RE.findall(txt)
    if not m0 or len(cyc) < 2:
        raise ValueError("timing blocks missing")
    ic, iw = m0.groups()
    (c1c, c1w), (c2c, c2w) = cyc[:2]
    return float(ic), iw, float(c1c), c1w, float(c2c), c2w

# New: parse the mesh-quality block
def _parse_grid_metrics(txt: str) -> Optional[dict]:
    a = MIN_ANGLE_RE.search(txt)
    b = MIN_EDGE_RE.search(txt)
    c = N_VERT_RE.search(txt)
    d = N_EDGES_RE.search(txt)
    e = N_FACES_RE.search(txt)
    f = AVG_EDGE_RE.search(txt)
    if not all([a, b, c, d, e, f]):
        return None
    return {
        "min_angle":        float(a.group(1)),
        "min_edge_length":  float(b.group(1)),
        "n_vertices":       int(c.group(1)),
        "n_edges":          int(d.group(1)),
        "n_faces":          int(e.group(1)),
        "avg_edge_length":  float(f.group(1)),
    }

# New: parse arrays after labels; supports values on one or multiple lines
def _extract_array_after_label(txt: str, label_re: re.Pattern) -> Optional[str]:
    m = label_re.search(txt)
    if not m:
        return None
    # Take up to two subsequent lines after the label; grab all ints
    lines = txt[m.end():].splitlines()
    vals: List[str] = []
    for k in range(min(3, len(lines))):  # safety window; usually 1 line is enough
        found = INT_LINE_RE.findall(lines[k])
        if found:
            vals.extend(found)
        # Stop if next line looks like another section header or blank line
        if k > 0 and (not found or lines[k].strip() == "" or "><" in lines[k]):
            break
    if not vals:
        return None
    return " ".join(vals)

# ───────────────────── main ────────────────────────────
def main() -> None:
    if not EXEC.exists():
        sys.exit(f"Executable missing: {EXEC}")
    originals = _load_decks()

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    timing_csv = BASE / f"scone_benchmark_{stamp}.csv"
    grid_csv   = BASE / f"grid_metrics_{stamp}.csv"
    vval_csv   = BASE / f"vertex_valence_{stamp}.csv"
    eval_csv   = BASE / f"edge_valence_{stamp}.csv"

    with timing_csv.open("w", newline="") as fh_t, \
         grid_csv.open("w", newline="")   as fh_g, \
         vval_csv.open("w", newline="")   as fh_v, \
         eval_csv.open("w", newline="")   as fh_e:

        w_t = csv.DictWriter(fh_t, fieldnames=TIMING_FIELDS)
        w_g = csv.DictWriter(fh_g, fieldnames=GRID_FIELDS)
        w_v = csv.DictWriter(fh_v, fieldnames=VVAL_FIELDS)
        w_e = csv.DictWriter(fh_e, fieldnames=EVAL_FIELDS)

        w_t.writeheader(); fh_t.flush()
        w_g.writeheader(); fh_g.flush()
        w_v.writeheader(); fh_v.flush()
        w_e.writeheader(); fh_e.flush()

        try:
            for geom in GEOM_CASES:
                deck = _tmpl(geom)
                for acc in ACCEL:
                    for pop in POP_VALUES:
                        tag = f"{geom:>22s} | {acc:<11s} | pop={pop:<7d}"
                        print(f"→ {tag} … ", end="", flush=True)

                        # Edit deck and run
                        deck.write_text(''.join(_edit(originals[deck], pop, geom, acc)))
                        out = _run(deck)

                        # Always try to parse grid metrics / valences, even on failure
                        gm = _parse_grid_metrics(out) if isinstance(out, str) and out else None
                        if gm:
                            w_g.writerow({"geometry": geom, **gm}); fh_g.flush()

                        vvals = _extract_array_after_label(out, VERTEX_LABEL_RE) if isinstance(out, str) and out else None
                        if vvals:
                            w_v.writerow({"geometry": geom, "vertexValence": vvals}); fh_v.flush()

                        evals = _extract_array_after_label(out, EDGE_LABEL_RE) if isinstance(out, str) and out else None
                        if evals:
                            w_e.writerow({"geometry": geom, "edgeValence": evals}); fh_e.flush()

                        # Handle timings CSV as before
                        if out == "OOM":
                            row = dict(
                                geometry=geom, acceleration=acc, pop=pop,
                                init_cpu_s="NA", init_wall_hms="NA",
                                cycle1_cpu_s="NA", cycle1_wall_hms="NA",
                                cycle2_cpu_s="NA", cycle2_wall_hms="NA",
                                note="terminated (out-of-memory)",
                            )
                            print("⚠ killed – continuing")
                        else:
                            try:
                                ic, iw, c1c, c1w, c2c, c2w = _parse_timings(out)
                                row = dict(
                                    geometry=geom, acceleration=acc, pop=pop,
                                    init_cpu_s=ic, init_wall_hms=iw,
                                    cycle1_cpu_s=c1c, cycle1_wall_hms=c1w,
                                    cycle2_cpu_s=c2c, cycle2_wall_hms=c2w,
                                    note="",
                                )
                                print("✓")
                            except Exception:
                                # No timings (e.g., early fatal) – still continue; note left blank
                                row = dict(
                                    geometry=geom, acceleration=acc, pop=pop,
                                    init_cpu_s="NA", init_wall_hms="NA",
                                    cycle1_cpu_s="NA", cycle1_wall_hms="NA",
                                    cycle2_cpu_s="NA", cycle2_wall_hms="NA",
                                    note="no timing section",
                                )
                                print("… (no timings)")

                        w_t.writerow(row); fh_t.flush()

        finally:
            # Restore pristine input decks even on Ctrl-C
            for p, lines in originals.items():
                p.write_text(''.join(lines))

    print("\nAll experiments complete.")
    print(f"Timings CSV       → {timing_csv}")
    print(f"Grid metrics CSV  → {grid_csv}")
    print(f"vertexValence CSV → {vval_csv}")
    print(f"edgeValence CSV   → {eval_csv}")
    

if __name__ == "__main__":
    main()
