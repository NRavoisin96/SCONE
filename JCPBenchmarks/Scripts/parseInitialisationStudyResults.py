#!/usr/bin/env python3
"""
Parse the initialisation-cost study of Section 4.3.4 (naive vs optimised).

Reads the SCONE outputs produced by runPatchSearchSweep.py for the
patch_naive_init and patch_optimised_init configurations. Output files are
discovered recursively by content; the mode is inferred from whether 'naive'
appears in the path.

Inputs: Results/InitialisationStudy/patch_naive_init/<geometry>/...
        Results/InitialisationStudy/patch_optimised_init/<geometry>/...

Outputs: init_study_long.csv -- one row per (mode, geometry, Dmax, run)
         init_factors.csv    -- per (geometry, Dmax): optimised and naive
                                initialisation means, the naive/optimised
                                factor with propagated error, and the
                                storage-identity verdict.
                                Feeds Fig. 12 and Tables M.17-Q.21.
         console report      -- storage-identity check plus per-depth factor
                                ranges.

The storage-identity check is the consistency check of Section 4.3.4: both
configurations must produce byte-identical acceleration structures, so any
mismatch in reported storage means the two are not comparable and the measured
factors are meaningless.

Usage: python3 parseInitialisationStudyResults.py [Results/InitialisationStudy]
"""

import re
import sys
import math
import pathlib
import numpy as np
import pandas as pd

RESULTS = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "Results/InitialisationStudy")

ARRAY_RE = re.compile(r"rawPatchInitTimes_D(\d+)_Res\s*=\s*\[([^\]]*)\]", re.S)
STORE_RE = re.compile(r"patchStorageSize_D(\d+)\s*=\s*([-+0-9.EeDd]+)")

def to_floats(raw: str) -> list[float]:
  return [float(v.replace("D", "E").replace("d", "e")) for v in re.split(r"[,\s]+", raw.strip()) if v]


def discover(root: pathlib.Path):
  hits = []
  for f in sorted(root.rglob("*")):
    if not f.is_file():
      continue
    try:
      text = f.read_text(errors="ignore")
    except OSError:
      continue
    if "rawPatchInitTimes" in text:
      mode = ("naive" if any("naive" in p.lower() for p in f.parts)
    else "optimised")
      hits.append((f, mode, f.parent.name))
  return hits

def main() -> None:
  if not RESULTS.exists():
    sys.exit(f"Error: '{RESULTS}' does not exist. Pass the results directory as the first argument.")
  files = discover(RESULTS)
  if not files:
    sys.exit(f"Error: no initialisation study outputs found under '{RESULTS.resolve()}'.")

  long_rows = []
  storage = {}
  for f, mode, geom in files:
    text = f.read_text()
    for d, s in STORE_RE.findall(text):
      storage[(mode, geom, int(d))] = float(s.replace("D", "E"))
    for d, raw in ARRAY_RE.findall(text):
      for r, t in enumerate(to_floats(raw), start = 1):
        long_rows.append({"mode": mode, "geometry": geom, "Dmax": int(d), "run": r, "init_s": t})

  long_df = pd.DataFrame(long_rows)
  long_df.to_csv("init_study_long.csv", index = False)

  stats = (long_df.groupby(["mode", "geometry", "Dmax"], as_index = False)
           .agg(init_mean = ("init_s", "mean"), init_std = ("init_s", lambda x: x.std(ddof = 1)), n_runs = ("init_s", "size")))

  rows = []
  print(f"{'geometry':>9} {'D':>2} | {'t_opt (s)':>10} {'t_naive (s)':>12} | {'factor':>7} | storage")
  print("-" * 62)
  for (geom, d), sub in stats.groupby(["geometry", "Dmax"]):
    try:
      opt = sub[sub["mode"] == "optimised"].iloc[0]
      nai = sub[sub["mode"] == "naive"].iloc[0]
    except IndexError:
      print(f"{geom:>9} {d:>2} | missing one mode -- skipped")
      continue
    factor = nai.init_mean / opt.init_mean
    ferr = factor * math.sqrt((nai.init_std / nai.init_mean) ** 2 + (opt.init_std / opt.init_mean) ** 2)
    s_opt = storage.get(("optimised", geom, d))
    s_nai = storage.get(("naive", geom, d))
    s_ok = (s_opt is not None and s_nai is not None and s_opt == s_nai)
    rows.append({"geometry": geom, "Dmax": d,
                 "init_opt_s": round(opt.init_mean, 5),
                 "init_opt_std": round(opt.init_std, 5),
                 "init_naive_s": round(nai.init_mean, 5),
                 "init_naive_std": round(nai.init_std, 5),
                 "factor": round(factor, 3),
                 "factor_err": round(ferr, 3),
                 "storage_MB": s_opt,
                 "storage_identical": s_ok})
    print(f"{geom:>9} {d:>2} | {opt.init_mean:10.4f} {nai.init_mean:12.4f} | {factor:6.2f}x | {'OK' if s_ok else 'MISMATCH'}")

  out = pd.DataFrame(rows)
  out.to_csv("init_factors.csv", index = False)

  bad = out[~out.storage_identical]
  print(f"\n[S] Storage identity: {'ALL OK' if bad.empty else f'{len(bad)} MISMATCHES -- LIST:'}")
  if not bad.empty:
    print(bad[["geometry", "Dmax"]].to_string(index = False))

  for d in sorted(out.Dmax.unique()):
    sub = out[out.Dmax == d]
    print(f"Initialisation factor range at Dmax = {d}: factors {sub.factor.min():.2f}-{sub.factor.max():.2f}x ({len(sub)} geometries).")

if __name__ == "__main__":
  main()