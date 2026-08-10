#!/usr/bin/env python3
"""
Benchmark sweep driver for the Patch-Search host-element determination study.

Reproduces the host-determination, ablation and initialisation measurements
reported in Sections 4.3.2-4.3.5 of the manuscript. Each (configuration,
geometry) pair is one rendered SCONE input file and one SCONE run; runs are
executed sequentially so that timings within a sweep are comparable, and the
sweep is resumable (completed runs are marked with a DONE file and skipped on
re-invocation).

Configurations
--------------
Each entry of CONFIGS selects an acceleration structure and the state of the
single-face intersection optimisation, and maps onto the manuscript as follows:

  patch_full            AMLG Patch-Search, all optimisations enabled.
                        Configuration (a) of Section 4.3.5; the main
                        benchmark of Section 4.3.2.
  patch_noshort         As above with the single-face optimisation disabled,
                        both as a refinement-termination criterion and as a
                        query-time branch. Configuration (b) of Section 4.3.5.
  standard_ascg         AMLG hierarchy retained, pseudo-angular sector
                        searches replaced by conventional point-in-element
                        tests. Configuration (c) of Section 4.3.5. (The identifier 
                        is a legacy harness name kept for consistency with the 
                        archived result directories; it does NOT denote the uniform
                        ASCG method of the original formulation.)
  patch_optimised_init  Initialisation-cost study of Section 4.3.4, using the
                        cached separating-axis quantities and the centroid-
                        based inclusion test.
  patch_naive_init      The same structure built with the naive procedures
                        (separating-axis quantities recomputed per test,
                        vertex-based inclusion test, no restriction to
                        unprocessed cells). Produces an identical
                        acceleration structure; only the build cost differs.

The two initialisation configurations use a reduced query population, since
only initialisation time is measured from them.

Requirements
------------
A SCONE build. The driver looks for the executable at ../Build/scone.out
relative to the working directory; set the SCONE_BIN environment variable to
override. Note that the mapping-frequency statistics of Section 4.2 require a
separate build (-DPATCH_SEARCH_STATS=ON -DOPENMP=OFF) and are collected by
patchSearchStatisticsParser.py, not by this driver.

Input templates
---------------
templates/<geometry>.inp are ordinary SCONE input files in which the
configuration-dependent entries have been replaced by placeholder tokens:

  pop                @POPULATION@;
  patchType          @PATCHTYPE@;
  naiveInit          @NAIVE@;
  singleFaceShortcut @SHORTCUT@;
  depths             @DEPTHS@;
  seeds              @SEEDS@;
  outputFile         @OUTPUT@;

Everything else (octree settings, geometry and nuclear data blocks) is left
as-is. A missing token is treated as an error rather than silently ignored.

Usage
-----
  python3 patchSearchSweeps.py --dry-run    # render inputs, run nothing
  python3 patchSearchSweeps.py              # full sweep (long; use tmux)

Outputs are written to results/<config>/<geometry>/, each containing the SCONE
output file, the rendered input and the captured stdout. results/sweep.log
records progress and results/sweep_manifest.csv the wall time per run.
"""

import argparse
import csv
import datetime as dt
import os
import pathlib
import subprocess
import sys
import time

SCONE_EXECUTABLE = os.environ.get("SCONE_BIN", "../Build/scone.out")
TEMPLATES_DIRECTORY = pathlib.Path("templates")
RESULTS_DIRECTORY = pathlib.Path("results")

GEOMETRIES = [
  "Hex72", "Hex243", "Poly264", "Tet137", "Tet298", "Tet427",
  "Poly468", "Hex576", "Poly940", "Hex1125",
  "Tet1331", "Tet1820", "Hex1944", "Poly1560"
]

# Maximum AMLG depths swept per geometry. Tet1331 omits Dmax = 1: the uniform
# single-layer grid exhausts memory on this geometry (Section 4.3.3).
DEPTHS_OVERRIDE = {"Tet1331": "(2 3 4 5)"}
DEFAULT_DEPTHS = "(1 2 3 4 5)"

# Twenty RNG seeds per configuration, shared across configurations so that
# comparisons between them are paired.
SEEDS = "(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)"

# (name, query population, acceleration structure, naive initialisation flag,
#  single-face optimisation flag, seeds). See the module docstring for how each
# configuration maps onto the manuscript.
CONFIGS = [
  ("patch_full",           "1000000000", "patchSearchAcceleration",  "0", "1", SEEDS),
  ("patch_noshort",        "1000000000", "patchSearchAcceleration",  "0", "0", SEEDS),
  ("standard_ascg",        "1000000000", "standardASCGAcceleration", "0", "1", SEEDS),
  ("patch_optimised_init", "1000000",    "patchSearchAcceleration",  "0", "1", SEEDS),
  ("patch_naive_init",     "1000000",    "patchSearchAcceleration",  "1", "1", SEEDS)
]

# Rough per-run estimate (seconds), used only for the ETA printout.
EST_SECONDS_PER_RUN = 90 * 60

# --------------------------------------------------------------------------
# Driver -------------------------------------------------------------------
# --------------------------------------------------------------------------

def log(msg: str) -> None:
  """Write a timestamped line to stdout and to results/sweep.log."""
  stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
  line = f"[{stamp}] {msg}"
  print(line, flush = True)
  with open(RESULTS_DIRECTORY / "sweep.log", "a") as fh:
    fh.write(line + "\n")

def render_template(template: pathlib.Path, dest: pathlib.Path, population: str, patch_type: str, naive_initialisation: str,
                    shortcut: str, depths: str, seeds: str, output_path: pathlib.Path) -> None:
  """Substitute the configuration tokens of a template into a SCONE input."""
  text = template.read_text()
  for token, value in [
    ("@POPULATION@", population),
    ("@PATCHTYPE@", patch_type),
    ("@NAIVE@", naive_initialisation),
    ("@SHORTCUT@", shortcut),
    ("@DEPTHS@", depths),
    ("@SEEDS@", seeds),
    ("@OUTPUT@", str(output_path)),
  ]:
    # An absent token means the template predates this driver: fail loudly
    # rather than run a configuration that was not actually applied.
    if token not in text:
      sys.exit(f"Error: token {token} missing from template {template}. Add the placeholder (see module docstring) and re-run.")
    text = text.replace(token, value)
  dest.write_text(text)

def run_one(config_name: str, population: str, patch_type: str, naive_initialisation: str, shortcut: str,
            seeds: str, geometry: str, dry_run: bool) -> float | None:
  """Run one (configuration, geometry) pair. Returns wall time, or None if skipped."""
  template = TEMPLATES_DIRECTORY / f"{geometry}.inp"
  if not template.exists():
    log(f"Skip {config_name} / {geometry}: template {template} not found.")
    return None

  run_dir = RESULTS_DIRECTORY / config_name / geometry
  done_marker = run_dir / "DONE"
  if done_marker.exists():
    log(f"Skip {config_name} / {geometry}: already complete.")
    return None

  run_dir.mkdir(parents = True, exist_ok = True)
  rendered = run_dir / f"{geometry}_{config_name}.inp"
  scone_output = run_dir / "output"
  depths = DEPTHS_OVERRIDE.get(geometry, DEFAULT_DEPTHS)
  render_template(template, rendered, population, patch_type, naive_initialisation, shortcut, depths, seeds, scone_output)

  if dry_run:
    log(f"DRY-RUN {config_name}/{geometry}: would execute {SCONE_EXECUTABLE} {rendered}.")
    return None

  log(f"START {config_name}/{geometry}.")
  t0 = time.monotonic()
  stdout_log = run_dir / "stdout.log"
  with open(stdout_log, "w") as out:
    proc = subprocess.run([SCONE_EXECUTABLE, str(rendered)], stdout = out, stderr = subprocess.STDOUT)
  elapsed = time.monotonic() - t0

  # A failed run is left without a DONE marker so that re-invoking the driver
  # retries it; the partial output and stdout log are kept for inspection.
  if proc.returncode != 0:
    log(f"FAIL {config_name} / {geometry} (return code = {proc.returncode}, {elapsed / 3600:.2f} h) -- see {stdout_log}.")
    return elapsed

  done_marker.write_text(dt.datetime.now().isoformat())
  log(f"DONE  {config_name} / {geometry} in {elapsed / 3600:.2f} h.")
  return elapsed

def main() -> None:
  parser = argparse.ArgumentParser()
  parser.add_argument("--dry-run", action = "store_true", help = "render inputs and print the plan; run nothing")
  args = parser.parse_args()

  RESULTS_DIRECTORY.mkdir(exist_ok = True)
  jobs = [(c, g) for c in CONFIGS for g in GEOMETRIES]
  pending = [
    (c, g) for (c, g) in jobs
    if not (RESULTS_DIRECTORY / c[0] / g / "DONE").exists()
  ]
  eta_h = len(pending) * EST_SECONDS_PER_RUN / 3600
  log(f"Sweep: {len(jobs)} jobs total, {len(pending)} pending, crude ETA {eta_h:.1f} h.")

  # The manifest is appended to and flushed per run, so that an interrupted
  # sweep still leaves a usable record of what completed.
  manifest_path = RESULTS_DIRECTORY / "sweep_manifest.csv"
  new_manifest = not manifest_path.exists()
  with open(manifest_path, "a", newline = "") as mf:
    writer = csv.writer(mf)
    if new_manifest:
      writer.writerow(["config", "geometry", "wall_seconds", "finished_at"])
    for (name, pop, ptype, naive, shortcut, seeds), geom in jobs:
      elapsed = run_one(name, pop, ptype, naive, shortcut, seeds, geom, args.dry_run)
      if elapsed is not None and not args.dry_run:
        writer.writerow([name, geom, f"{elapsed:.1f}", dt.datetime.now().isoformat()])
        mf.flush()

  log("Sweep finished.")

if __name__ == "__main__":
  main()