#!/usr/bin/env python3
"""
Octree sensitivity sweep driver.

Reproduces the octree parameter study of Section 4.3.1 (Figs. 7 and 8), which
establishes the optimal parameter pair {Dmax = 6, Nfaces,max = 1} used for the
octree baseline in all subsequent comparisons. One SCONE process is launched
per geometry; the octreeOptimisationPackage then loops over the full
depths x nMaxFaces x runs grid internally, so a single run produces the entire
sensitivity surface for that mesh.

Parameters
----------
The swept ranges and the query population are fixed here to the values reported
in the manuscript: 5e6 query points, Dmax in {4, ..., 8} and Nfaces,max in
{1, ..., 16}, with a fixed RNG seed per parameter pair. Ensemble sizes follow
the figure legends: 100 independent runs on the finer meshes (Tet1331, Tet1820,
Hex1125, Hex1944, Poly1560) and 50 on the remainder.

Requirements
------------
A SCONE build. The driver looks for the executable at ../Build/scone.out
relative to the working directory; set the SCONE_BIN environment variable to
override.

Input templates
---------------
Templates/OctreeOptimisation/<geometry>.inp are ordinary SCONE optimisation
input files in which the swept quantities have been replaced by placeholder
tokens:

  pop        @POPULATION@;
  seed       @SEED@;
  depths     @DEPTHS@;
  nMaxFaces  @NMAXFACES@;
  nRuns      @NRUNS@;
  outputFile @OUTPUT@;

Everything else (geometry and nuclear data blocks) is left as-is. A missing
token is treated as an error rather than silently ignored.

Usage
-----
  python3 octreeOptimisation.py --dry-run # render inputs, run nothing
  python3 octreeOptimisation.py           # full sweep (long; use tmux)

Outputs are written to Results/OctreeOptimisation/<geometry>/, each containing
the SCONE output file, the rendered input and the captured stdout. Completed
geometries drop a DONE marker and are skipped when the driver is re-invoked, so
an interrupted sweep can simply be restarted.
"""

import argparse
import datetime as dt
import os
import pathlib
import subprocess
import time

SCONE_EXECUTABLE = os.environ.get("SCONE_BIN", "../Build/scone.out")
TEMPLATE_DIRECTORY = pathlib.Path("Templates/OctreeOptimisation")
RESULTS_DIRECTORY = pathlib.Path("Results/OctreeOptimisation")

# Ensemble sizes per the published figure legends: 100 runs on the finer
# meshes, 50 elsewhere.
DEFAULT_NRUNS = 50
NRUNS_OVERRIDES = {"Tet1331": 100, "Tet1820": 100, "Hex1125": 100, "Hex1944": 100, "Poly1560": 100}

# Swept ranges and query population, fixed to the values reported in the
# manuscript. The seed is held fixed so that each {Dmax, Nfaces,max} pair is
# evaluated on the same query points.
POPULATION = 5000000
SEED = 1
DEPTHS = "(4 5 6 7 8)"
NMAXFACES = "(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)"

GEOMETRIES = [
  "Hex72", "Hex243", "Poly264", "Tet137", "Tet298", "Tet427",
  "Poly468", "Hex576", "Poly940", "Hex1125",
  "Tet1331", "Tet1820", "Hex1944", "Poly1560"
]

# --------------------------------------------------------------------------
# Driver -------------------------------------------------------------------
# --------------------------------------------------------------------------

def log(msg: str) -> None:
  """Write a timestamped line to stdout and to the sweep log."""
  stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
  line = f"[{stamp}] {msg}"
  print(line, flush = True)
  with open(RESULTS_DIRECTORY / "sweep.log", "a") as fh:
    fh.write(line + "\n")

def run_one(geometry: str, dry_run: bool) -> None:
  """Render and run the full sensitivity grid for one geometry."""
  template = TEMPLATE_DIRECTORY / f"{geometry}.inp"
  if not template.exists():
    log(f"Skip {geometry}: template not found.")
    return

  run_dir = RESULTS_DIRECTORY / geometry
  if (run_dir / "DONE").exists():
    log(f"Skip {geometry}: already complete.")
    return

  run_dir.mkdir(parents = True, exist_ok = True)
  nruns = NRUNS_OVERRIDES.get(geometry, DEFAULT_NRUNS)

  text = template.read_text()
  for token, value in [("@POPULATION@", str(POPULATION)),
                       ("@SEED@", str(SEED)),
                       ("@DEPTHS@", DEPTHS),
                       ("@NMAXFACES@", NMAXFACES),
                       ("@NRUNS@", str(nruns)),
                       ("@OUTPUT@", str(run_dir / "output"))]:
    # An absent token means the template predates this driver: fail loudly
    # rather than run a configuration that was not actually applied.
    if token not in text:
      raise SystemExit(f"Token {token} missing in {template}.")
    text = text.replace(token, value)

    rendered = run_dir / f"{geometry}.inp"
    rendered.write_text(text)

    if dry_run:
      log(f"DRY-RUN {geometry} (nRuns = {nruns}).")
      return

    log(f"START {geometry} (nRuns = {nruns}).")
    t0 = time.monotonic()
    with open(run_dir / "stdout.log", "w") as out:
      proc = subprocess.run([SCONE_EXECUTABLE, str(rendered)], stdout = out, stderr = subprocess.STDOUT)
    elapsed = (time.monotonic() - t0) / 3600

    # A failed run is left without a DONE marker so that re-invoking the driver
    # retries it; the partial output and stdout log are kept for inspection.
    if proc.returncode != 0:
      log(f"FAIL {geometry} (return code = {proc.returncode}, {elapsed:.2f} h).")
      return

    (run_dir / "DONE").write_text(dt.datetime.now().isoformat())
    log(f"DONE {geometry} in {elapsed:.2f} h.")

def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--dry-run", action = "store_true", help = "render inputs and print the plan; run nothing")
  args = parser.parse_args()

  RESULTS_DIRECTORY.mkdir(exist_ok = True)
  pending = [g for g in GEOMETRIES
              if not (RESULTS_DIRECTORY / g / "DONE").exists()]
  log(f"Octree sensitivity sweep: {len(pending)} of {len(GEOMETRIES)} geometries pending.")

  for g in GEOMETRIES:
    run_one(g, args.dry_run)

  log("Optimisation sweep complete.")

if __name__ == "__main__":
  main()