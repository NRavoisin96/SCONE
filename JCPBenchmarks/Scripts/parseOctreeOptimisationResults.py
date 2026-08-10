"""
Parse the octree sensitivity sweep of Section 4.3.1 (Figs. 7 and 8).

Reads the SCONE outputs produced by runOctreeOptimisation.py, one per geometry, and
determines the optimal {Dmax, Nfaces,max} pair for each mesh. This is what
establishes the {6, 1} parameter pair used for the octree baseline in all
subsequent comparisons.

Inputs: Results/OctreeOptimisation/<geometry>/output.m  (asciiMATLAB, with
        raw blocks rawHostTimes_D{d}_N{n}, rawInitTimes_D{d}_N{n} and
        storageSize_D{d}_N{n})

Outputs: octree_sensitivity_long.csv  -- one row per (geometry, Dmax, Nfaces, run)
         octree_sensitivity_stats.csv -- mean/std per (geometry, Dmax, Nfaces); feeds Figs. 7 and 8
         console report               -- per-geometry optimal pair, the
                                         published {6, 1} timing, the margin
                                         between them, and whether {6, 1}
                                         remains optimal or lies within the
                                         combined standard deviations.

Usage: python3 parseOctreeOptimisationResults.py [Results/OctreeOptimisation]
"""

import re
import sys
import pathlib
import numpy as np
import pandas as pd

RESULTS = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "Results/OctreeOptimisation")
PUBLISHED_OPTIMUM = (6, 1) # Section 4.3.1: optimal pair, invariant across all topologies.

BLOCK_RE = re.compile(r"raw(Host|Init)Times_D(\d+)_N(\d+)\w*\s*=\s*\[([^\]]*)\]", re.S)
STORE_RE = re.compile(r"storageSize_D(\d+)_N(\d+)\s*=\s*([-+0-9.Ee]+)")

def parse_geometry(path: pathlib.Path, geometry: str) -> pd.DataFrame:
  text = path.read_text()
  frames = {}
  for kind, d, n, raw in BLOCK_RE.findall(text):
      vals = [float(v) for v in re.split(r"[,\s]+", raw.strip()) if v]
      key = (int(d), int(n))
      frames.setdefault(key, {})[kind.lower()] = vals
  storage = {(int(d), int(n)): float(s)
              for d, n, s in STORE_RE.findall(text)}
  rows = []
  for (d, n), data in frames.items():
    host = data.get("host", [])
    init = data.get("init", [])
    # host and init arrays are emitted per parameter pair and are the same length; zip() pairs run k of each.
    for k, (h, i) in enumerate(zip(host, init), start=1):
      rows.append({"geometry": geometry, "Dmax": d, "Nfaces": n, "run": k, "host_s": h, "init_s": i,
                   "storage_MB": storage.get((d, n), np.nan)})
  return pd.DataFrame(rows)

def main() -> None:
  frames = []
  for gdir in sorted(RESULTS.iterdir()):
    f = gdir / "output.m"
    if gdir.is_dir() and f.exists():
      frames.append(parse_geometry(f, gdir.name))
  if not frames:
    raise SystemExit(f"No outputs found under {RESULTS}.")
  long_df = pd.concat(frames, ignore_index = True)
  long_df.to_csv("octree_sensitivity_long.csv", index = False)

  stats = (long_df
            .groupby(["geometry", "Dmax", "Nfaces"], as_index = False)
            .agg(host_mean = ("host_s", "mean"),
                host_std = ("host_s", lambda x: x.std(ddof = 1)),
                init_mean = ("init_s", "mean"),
                init_std = ("init_s", lambda x: x.std(ddof = 1)),
                storage_MB = ("storage_MB", "first"),
                n_runs = ("run", "size")))
  stats.to_csv("octree_sensitivity_stats.csv", index = False)

  print(f"{'geometry':>9} | optimum (D,N) | t_opt (s) | t(6,1) (s) | margin | {PUBLISHED_OPTIMUM} still optimal?")
  print("-" * 78)
  for g, sub in stats.groupby("geometry", sort = False):
    best = sub.loc[sub.host_mean.idxmin()]
    pub = sub[(sub.Dmax == PUBLISHED_OPTIMUM[0]) (sub.Nfaces == PUBLISHED_OPTIMUM[1])]
    if pub.empty:
      print(f"{g:>9} | published config not in grid -- check input")
      continue
    pub = pub.iloc[0]
    margin = pub.host_mean - best.host_mean
    within = margin <= best.host_std + pub.host_std
    verdict = ("YES" if (best.Dmax, best.Nfaces) == PUBLISHED_OPTIMUM
                else ("within 1 sigma" if within else "NO -- MOVED"))
    print(f"{g:>9} | ({int(best.Dmax)},{int(best.Nfaces)})"
          f"{'':<6} | {best.host_mean:8.3f} | {pub.host_mean:9.3f} | "
          f"{margin:6.3f} | {verdict}")

if __name__ == "__main__":
  main()