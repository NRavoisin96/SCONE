#!/usr/bin/env python3
"""
Parse the mapping-frequency statistics of Section 4.2 (counters, terminal-cell
volume tallies and raw times).

Reads the SCONE outputs produced by runPatchSearchSweeps.py for the patch_full
configuration, built with -DPATCH_SEARCH_STATS=ON -DOPENMP=OFF (the counters
are not thread-safe; see the build notes in the README).

Inputs: Results/MappingStatistics/patch_full/<geometry>/output.m  (asciiMATLAB)

Outputs: mapping_stats_long.csv  -- one row per (geometry, Dmax, run)
         mapping_stats_table.csv -- per (geometry, Dmax): dynamic trigger
                                    fractions (mean +/- std over runs),
                                    static terminal-cell volume fractions,
                                    their difference in percentage points,
                                    and storage. Feeds Table 4 and
                                    Appendices B-F.
         console report          -- consistency checks:
           [A] Termination invariant  (outside + psi + single_face + angular == queries - displacements, 
                                       integer-exact, every run)
           [B] Cross-depth invariance (unique queries per seed identical across all depths of a geometry)
           [C] Volume closure         (category volumes sum to the reported total, to float precision)

Displaced queries (phi) re-descend the structure and terminate via another
mechanism, which is why they are subtracted from the query count rather than
counted as a terminal category -- [A] is the machine-checkable form of that statement.

Usage: python3 parseMappingStatisticsResults.py [Results/MappingStatistics/patch_full]
"""

import re
import sys
import pathlib
import numpy as np
import pandas as pd

RESULTS = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "Results/MappingStatistics/patch_full")

NUM = r"[-+0-9.EeDd]+"
SCALAR_RE = re.compile(rf"^\s*(\w+)\s*=\s*({NUM})\s*;", re.M)
ARRAY_RE = re.compile(r"(\w+)_Res\s*=\s*\[([^\]]*)\]", re.S)

COUNTERS = {
  "rawPatchNQueries": "queries",
  "rawPatchNOutside": "outside",
  "rawPatchNDirectElement": "psi",
  "rawPatchNSingleFace": "single_face",
  "rawPatchNAngularSearch": "angular",
  "rawPatchNVertexDisplacement": "displacement"
}
TIMES = {"rawPatchInitTimes": "init_s", "rawPatchHostTimes": "host_s"}
VOLUMES = {
  "patchElementMappingVolume": "vol_psi",
  "patchSingleFaceVolume": "vol_single_face",
  "patchEdgeMappingVolume": "vol_edge",
  "patchVertexMappingVolume": "vol_vertex",
  "patchOutsideVolume": "vol_outside",
  "patchTotalVolume": "vol_total"
}

def to_floats(raw: str) -> list[float]:
  return [float(v.replace("D", "E").replace("d", "e")) for v in re.split(r"[,\s]+", raw.strip()) if v]

def parse_geometry(path: pathlib.Path, geometry: str):
  text = path.read_text()
  scalars = {k: float(v.replace("D", "E")) for k, v in SCALAR_RE.findall(text)}
  arrays: dict[str, list[float]] = {k: to_floats(raw) for k, raw in ARRAY_RE.findall(text)}

  depth_tags = sorted({int(m.group(1)) for m in re.finditer(r"rawPatchNQueries_D(\d+)", text)})
  long_rows, table_rows, checks = [], [], []

  unique_by_depth = {}
  for d in depth_tags:
    counts = {name: np.array(arrays[f"{key}_D{d}"], dtype = np.int64) for key, name in COUNTERS.items()}
    times = {name: np.array(arrays[f"{key}_D{d}"]) for key, name in TIMES.items()}
    n_runs = len(counts["queries"])
    unique = counts["queries"] - counts["displacement"]
    unique_by_depth[d] = unique

    # [A] Termination invariant, integer-exact per run
    term = (counts["outside"] + counts["psi"] + counts["single_face"] + counts["angular"])
    inv_ok = bool(np.array_equal(term, unique))

    # [C] Volume closure
    vols = {name: scalars.get(f"{key}_D{d}", np.nan) for key, name in VOLUMES.items()}
    closure = (vols["vol_psi"] + vols["vol_single_face"] + vols["vol_edge"] + vols["vol_vertex"] +
               vols["vol_outside"]) - vols["vol_total"]
    closure_ok = bool(abs(closure) < 1e-5 * max(vols["vol_total"], 1.0))

    checks.append((d, inv_ok, closure_ok, closure))

    for r in range(n_runs):
      long_rows.append({
        "geometry": geometry, "Dmax": d, "run": r + 1, "unique_queries": unique[r],
        **{k: counts[k][r] for k in counts}, **{k: times[k][r] for k in times},
      })

    # Dynamic fractions (% of unique queries), static fractions (% of in-mesh volume)
    dyn = {k: 100.0 * counts[k] / unique for k in ("psi", "single_face", "angular", "displacement")}
    in_mesh = vols["vol_total"] - vols["vol_outside"]

    # Static fractions are keyed by terminal-cell mapping: edge-mapped volume
    # corresponds to the angular-sector branch, vertex-mapped to displacement.
    stat = {
      "psi": 100.0 * vols["vol_psi"] / in_mesh,
      "single_face": 100.0 * vols["vol_single_face"] / in_mesh,
      "angular": 100.0 * vols["vol_edge"] / in_mesh,
      "displacement": 100.0 * vols["vol_vertex"] / in_mesh,
    }
    row = {"geometry": geometry, "Dmax": d, "storage_MB": scalars.get(f"patchStorageSize_D{d}", np.nan),
           "invariant_ok": inv_ok, "closure_ok": closure_ok}
    for k in ("psi", "single_face", "angular", "displacement"):
      row[f"dyn_{k}_pct"] = dyn[k].mean()
      row[f"dyn_{k}_std"] = dyn[k].std(ddof = 1)
      row[f"stat_{k}_pct"] = stat[k]
      row[f"delta_{k}_pp"] = dyn[k].mean() - stat[k]
    table_rows.append(row)

  # [B] Cross-depth invariance of unique queries per seed
  ds = list(unique_by_depth)
  cross_ok = all(np.array_equal(unique_by_depth[ds[0]], unique_by_depth[d]) for d in ds[1:])
  return long_rows, table_rows, checks, cross_ok

def main() -> None:
  all_long, all_table = [], []
  print(f"{'geometry':>9} | cross-depth | per-depth [A]invariant [C]closure")
  print("-" * 72)
  for gdir in sorted(RESULTS.iterdir()):
    f = gdir / "output.m"
    if not (gdir.is_dir() and f.exists()):
      continue
    long_rows, table_rows, checks, cross_ok = parse_geometry(f, gdir.name)
    all_long += long_rows
    all_table += table_rows
    detail = "  ".join(
      f"D{d}:{'ok' if a else 'FAIL'}/{'ok' if c else f'FAIL({res:+.1e})'}"
      for d, a, c, res in checks)
    print(f"{gdir.name:>9} | {'ok' if cross_ok else 'FAIL':>11} | {detail}")

  pd.DataFrame(all_long).to_csv("mapping_stats_long.csv", ifndex = False)
  table = pd.DataFrame(all_table)
  table.to_csv("mapping_stats_table.csv", index = False)
  print(f"\nWrote mapping_stats_long.csv ({len(all_long)} rows) and "
        f"mapping_stats_table.csv ({len(all_table)} rows).")

  # Console preview at Dmax = 3, the general-purpose default recommended in Section 5; the CSV carries every depth.
  z = table[table.Dmax == 3]
  if not z.empty:
    print("\nTable 5 preview (Dmax = 3): dynamic % (static %) per category.")
    for _, r in z.iterrows():
      print(f" {r.geometry:>9}: psi {r.dyn_psi_pct:6.2f} "
            f"({r.stat_psi_pct:6.2f}) | sf {r.dyn_single_face_pct:6.2f} "
            f"({r.stat_single_face_pct:6.2f}) | "
            f"ang {r.dyn_angular_pct:6.4f} ({r.stat_angular_pct:6.4f}) | "
            f"disp {r.dyn_displacement_pct:6.4f} "
            f"({r.stat_displacement_pct:6.4f})")

if __name__ == "__main__":
  main()