# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Sweep maxFacesNumber × maxRefinementLevel inside octree_class.f90,
re-compile SCONE for every combination, run one benchmark, and stream the
timings directly to CSV (so Python never holds the whole matrix in RAM).

Outer loop  : maxFacesNumber
Inner loop  : maxRefinementLevel
"""

from __future__ import annotations
import csv, datetime as dt, re, shutil, subprocess, sys
from pathlib import Path
from typing import List, Tuple

# ─────────────── CONTROL PANEL ──────────────────────────────────
FACES_VALUES   = [1, 4, 7, 9, 11, 13, 15, 17, 19]    # outer loop
REFINE_VALUES  = [2, 4, 6, 8, 10, 12, 14, 16, 18]                     # inner loop

BASE_DIR   = Path(__file__).resolve().parent      # /home/dk725/SCONE/SCONE
BUILD_DIR  = BASE_DIR / "Build"
EXECUTABLE = BUILD_DIR / "scone.out"

# choose whichever deck you normally run after (re)building SCONE
DECK_FILE  = BASE_DIR / "InputFiles" / "SCONE_ToyProblemTet1"

# *** this is the file that contains the parameters on line 20 ***
TREE_SOURCE = BASE_DIR / "DataStructures" / "Trees" / "octree_class.f90"
# ────────────────────────────────────────────────────────────────

# regexes to patch the variables irrespective of spacing / case
FACES_RE   = re.compile(r"(maxFacesNumber\s*=\s*)\d+", re.I)
REFINE_RE  = re.compile(r"(maxRefinementLevel\s*=\s*)\d+", re.I)

# regexes to extract timings from SCONE output
INIT_RE  = re.compile(
    r"Initialisation procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)",
    re.S,
)
CYCLE_RE = re.compile(
    r"In-cycle procedure time.*?CPU\s+time:\s+([0-9.E+-]+).*?Wall time:\s+([0-9:]+)",
    re.S,
)
KILL_RC  = {9, -9, 137, -137}          # SIGKILL & 128+SIGKILL

CSV_FIELDS = [
    "maxFacesNumber", "maxRefLevel",
    "init_cpu_s", "init_wall_hms",
    "cycle1_cpu_s", "cycle1_wall_hms",
    "cycle2_cpu_s", "cycle2_wall_hms",
    "note",
]

# ───────────────────── helpers ────────────────────────────────
def backup_original(src: Path) -> List[str]:
    text = src.read_text().splitlines(keepends=True)
    bak = src.with_suffix(".bak")
    if not bak.exists():
        shutil.copy2(src, bak)
    return text

def patch_source(lines: List[str], faces: int, refine: int) -> List[str]:
    """Return NEW source lines with maxFacesNumber and maxRefinementLevel swapped in."""
    out = lines.copy()
    for i, l in enumerate(out):
        if "maxFacesNumber" in l:
            # use a lambda so we don’t generate “\14…” backreferences
            out[i] = FACES_RE.sub(lambda m: m.group(1) + str(faces), l)
        if "maxRefinementLevel" in l:
            out[i] = REFINE_RE.sub(lambda m: m.group(1) + str(refine), l)
    return out

def compile_scone() -> bool:
    """Return True on successful rebuild, False otherwise."""
    try:
        subprocess.run(
            ["make", "-C", str(BUILD_DIR)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"\n❌ Compilation failed:\n{e.stdout}\n")
        return False

def run_scone() -> Tuple[str, str]:
    """Return ('OK', stdout)  or ('OOM', '' )  or ('FAIL', stdout)."""
    try:
        out = subprocess.run(
            [str(EXECUTABLE), str(DECK_FILE), "--omp", "1"],
            cwd=BASE_DIR,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        ).stdout
        return "OK", out
    except subprocess.CalledProcessError as e:
        if e.returncode in KILL_RC or "Killed" in e.stdout:
            return "OOM", ""
        sys.stderr.write(f"\n❌ Runtime error:\n{e.stdout}\n")
        return "FAIL", e.stdout

def parse_timings(stdout: str) -> Tuple[float, str, float, str, float, str]:
    m0 = INIT_RE.search(stdout)
    cycles = CYCLE_RE.findall(stdout)
    if not m0 or len(cycles) < 2:
        raise ValueError("Timing blocks missing")
    ic, iw = m0.groups()
    (c1c, c1w), (c2c, c2w) = cycles[:2]
    return float(ic), iw, float(c1c), c1w, float(c2c), c2w

# ───────────────────── main driver ─────────────────────────────
def main() -> None:
    if not EXECUTABLE.exists():
        sys.exit(f"Executable not found at {EXECUTABLE}")

    pristine_src = backup_original(TREE_SOURCE)
    csv_path = BASE_DIR / f"treeSweepTet137_{dt.datetime.now():%Y%m%d_%H%M%S}.csv"

    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader(); fh.flush()

        try:
            for faces in FACES_VALUES:           # outer
                for ref in REFINE_VALUES:        # inner
                    tag = f"faces={faces:<7d} | ref={ref}"
                    print(f"→ {tag} … ", end="", flush=True)

                    # 1 ─ patch .f90
                    TREE_SOURCE.write_text("".join(patch_source(pristine_src, faces, ref)))

                    # 2 ─ build
                    if not compile_scone():
                        writer.writerow(dict(
                            maxFacesNumber=faces, maxRefLevel=ref,
                            init_cpu_s="NA", init_wall_hms="NA",
                            cycle1_cpu_s="NA", cycle1_wall_hms="NA",
                            cycle2_cpu_s="NA", cycle2_wall_hms="NA",
                            note="compilation failed",
                        )); fh.flush()
                        print("✗ compile") ; continue

                    # 3 ─ run
                    status, out = run_scone()
                    if status == "OOM":
                        writer.writerow(dict(
                            maxFacesNumber=faces, maxRefLevel=ref,
                            init_cpu_s="NA", init_wall_hms="NA",
                            cycle1_cpu_s="NA", cycle1_wall_hms="NA",
                            cycle2_cpu_s="NA", cycle2_wall_hms="NA",
                            note="terminated (out-of-memory)",
                        )); fh.flush()
                        print("⚠ OOM") ; continue
                    elif status != "OK":
                        writer.writerow(dict(
                            maxFacesNumber=faces, maxRefLevel=ref,
                            init_cpu_s="NA", init_wall_hms="NA",
                            cycle1_cpu_s="NA", cycle1_wall_hms="NA",
                            cycle2_cpu_s="NA", cycle2_wall_hms="NA",
                            note="runtime failure",
                        )); fh.flush()
                        print("✗ runtime") ; continue

                    # 4 ─ parse timings
                    try:
                        ic, iw, c1c, c1w, c2c, c2w = parse_timings(out)
                        writer.writerow(dict(
                            maxFacesNumber=faces, maxRefLevel=ref,
                            init_cpu_s=ic, init_wall_hms=iw,
                            cycle1_cpu_s=c1c, cycle1_wall_hms=c1w,
                            cycle2_cpu_s=c2c, cycle2_wall_hms=c2w,
                            note="",
                        )); fh.flush()
                        print("✓")
                    except Exception as err:
                        writer.writerow(dict(
                            maxFacesNumber=faces, maxRefLevel=ref,
                            init_cpu_s="NA", init_wall_hms="NA",
                            cycle1_cpu_s="NA", cycle1_wall_hms="NA",
                            cycle2_cpu_s="NA", cycle2_wall_hms="NA",
                            note=f"parsing error: {err}",
                        )); fh.flush()
                        print("✗ parse")

        finally:
            # restore original .f90 even on Ctrl-C
            TREE_SOURCE.write_text("".join(pristine_src))

    print(f"\nSweep finished. Results → {csv_path}")

if __name__ == "__main__":
    main()
