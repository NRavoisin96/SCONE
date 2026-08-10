#!/usr/bin/env python3
"""
Parse the unified measurement campaign of Section 4.3 (asciiMATLAB outputs).

Reads the SCONE outputs produced by runPatchSearchSweep.py for the three
host-determination configurations, and derives both the headline savings and
the ablation decomposition. Configurations are assigned from the provenance
echoes in each output file (patchType and singleFaceShortcut) rather than from
directory names, so a misfiled run is detected rather than silently mislabelled:

  patch_full    configuration (a), all depths, Tables H.12-L.16
  patch_noshort configuration (b), Dmax = 3,   Table R.22
  standard_ascg configuration (c), Dmax = 3,   Table R.22
                  (legacy name; this is the AMLG hierarchy with conventional
                   point-in-element tests, not the uniform ASCG method of the
                   original formulation)

Each output also carries the octree baseline runs for that geometry; these are
pooled across configurations to give the Table G.11 baseline.

Runs whose host time exceeds the array median by more than DROP_TOL are
excluded, one-sided, as background-activity contamination (Section 4.3);
the paired initialisation entry is dropped with it to preserve seed pairing.
Every exclusion is echoed in the console report.

Inputs: one or more results roots, e.g. Results/MainCampaign and
        Results/PerformanceDecompositionStudy; files are discovered recursively by content.

Outputs: campaign_stats.csv -- per (config, geometry, Dmax): init/host/total
                               mean and std, storage, retained run count.
         octree_stats.csv   -- per (config, geometry): octree init/host stats.
         decomposition.csv  -- Dmax = 3 ablation: per-mechanism savings with
                               propagated errors, memory and initialisation
                               inflation factors. Feeds Fig. 13 and Table R.22.
         console report     -- [GATE]   completeness, provenance echoes, exclusions and outliers.
                               [DECOMP] per-geometry decomposition and the head-to-head verdict counts.
                               [RANGES] optimal-depth savings versus the pooled octree baseline, and the
                                        depth-3 versus single-layer comparison.

Usage: python3 parseMainCampaignResults.py Results/MainCampaign Results/PerformanceDecompositionStudy
"""

import csv
import math
import re
import sys
import pathlib
import statistics as st

DEFAULT_ROOTS = [pathlib.Path("Results/MainCampaign"), pathlib.Path("Results/PerformanceDecompositionStudy")]
ROOTS = ([pathlib.Path(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else DEFAULT_ROOTS)
GEOM_ORDER = ["Tet137", "Tet298", "Tet427", "Tet1331", "Tet1820",
              "Hex72", "Hex243", "Hex576", "Hex1125", "Hex1944",
              "Poly264", "Poly468", "Poly940", "Poly1560"]
EXPECTED_FLAG = {"patch_full": 1, "patch_noshort": 0}
OUTLIER_TOL = 0.05 # Report-only: runs deviating >5% from the median.
DROP_TOL = 0.15    # exclusion criterion of Section 4.3.

ARRAY_RE = re.compile(r"^\s*(\w+)\s*=\s*\[([^\]]*)\]\s*;", re.M | re.S)
SCALAR_RE = re.compile(r"^\s*(\w+)\s*=\s*([-+0-9.EeDd]+)\s*;", re.M)
PTYPE_RE = re.compile(r"patchType\s*=\s*\{'(\w+)'\}")

def to_floats(raw: str) -> list[float]:
  return [float(v.replace("D", "E").replace("d", "e")) for v in re.split(r"[,\s]+", raw.strip()) if v]


def parse_file(path: pathlib.Path) -> dict:
  text = path.read_text()
  d: dict = {}
  for key, raw in ARRAY_RE.findall(text):
    d[key] = to_floats(raw)
  for key, raw in SCALAR_RE.findall(text):
    if key not in d:
      d[key] = float(raw.replace("D", "E").replace("d", "e"))
  m = PTYPE_RE.search(text)
  if m:
    d["patchType"] = m.group(1)
  return d

def stats(vals: list[float]) -> tuple[float, float]:
  return (st.mean(vals), st.stdev(vals) if len(vals) > 1 else 0.0)


def outliers(vals: list[float]) -> list[int]:
  med = st.median(vals)
  return [i + 1 for i, v in enumerate(vals) if med and abs(v - med) / med > OUTLIER_TOL]

def drop_contaminated(host: list[float], init: list[float]):
  """One-sided exclusion: drop runs > DROP_TOL above the array median
  (host-activity contamination only ever slows runs). Drops the paired
  init entry too, preserving seed pairing."""
  med = st.median(host)
  keep = [i for i, v in enumerate(host) if not (med and v > med * (1 + DROP_TOL))]
  dropped = [i + 1 for i in range(len(host)) if i not in keep]
  kept_init = [init[i] for i in keep] if init else init
  return [host[i] for i in keep], kept_init, dropped

def ratio_with_err(a, sa, b, sb):
  """r = (a - b)/a * 100 with first-order error propagation."""
  if a == 0:
    return float("nan"), float("nan")
  r = (a - b) / a
  sr = math.sqrt((b / (a * a) * sa) ** 2 + (sb / a) ** 2)
  return 100 * r, 100 * sr


def main() -> None:
  CONFIG_ALIASES = {"patch_full": "patch_full", "patch": "patch_full",
                    "full": "patch_full",
                    "patch_noshort": "patch_noshort",
                    "no_short": "patch_noshort",
                    "noshort": "patch_noshort",
                    "no_shortcut": "patch_noshort",
                    "standard_ascg": "standard_ascg",
                    "std_ascg": "standard_ascg",
                    "ascg": "standard_ascg",
                    "standardascg": "standard_ascg"}
  files = []
  raw_names = set()
  for root in ROOTS:
    for f in sorted(root.rglob("*")):
      if not f.is_file():
        continue
      text = f.read_text(errors = "ignore")
      if "rawPatchInitTimes" not in text:
        continue
      ptype_m = PTYPE_RE.search(text)
      flag_m = re.search(r"singleFaceShortcut\s*=\s*(\d)", text)
      ptype = ptype_m.group(1) if ptype_m else ""
      if "standard" in ptype.lower():
        config = "standard_ascg"
      elif flag_m and flag_m.group(1) == "0":
        config = "patch_noshort"
      else:
        config = "patch_full"
      raw = f.parent.parent.name
      raw_names.add(raw)
      path_cfg = CONFIG_ALIASES.get(raw.lower())
      if path_cfg and path_cfg != config:
        print(f"WARNING: {f} sits under '{raw}' but its provenance echo says {config}. Trusting the echo.")
      files.append((f, config, f.parent.name))
  if not files:
    sys.exit("ERROR: no campaign outputs found under: "
             + ", ".join(str(r.resolve()) for r in ROOTS))
  print(f"Discovered config directories: {sorted(raw_names)}.")
  print("Config assignment: from provenance echoes (patchType + singleFaceShortcut).")

  gate_msgs = []
  rows, oct_rows = [], []
  seen = {}
  for f, config, geom in files:
    d = parse_file(f)
    # -- provenance echo checks ------------------------------------
    flag = d.get("singleFaceShortcut")
    want = EXPECTED_FLAG.get(config)
    if want is not None and flag is not None and int(flag) != want:
      gate_msgs.append(f"PROVENANCE MISMATCH: {config}/{geom} echoes singleFaceShortcut = {int(flag)}.")
    ptype = d.get("patchType", "")
    if config == "standard_ascg" and "standard" not in ptype.lower():
      gate_msgs.append(f"PROVENANCE MISMATCH: {config}/{geom} echoes patchType = '{ptype}'.")
    # -- octree runs ----------------------------------------------
    oh = d.get("rawOctreeHostTimes_Res", [])
    oi = d.get("rawOctreeInitTimes_Res", [])
    if oh:
      oh, oi, dropped = drop_contaminated(oh, oi)
      for r in dropped:
        gate_msgs.append(f"EXCLUDED: {config}/{geom} octree host run {r} (>{DROP_TOL:.0%} above median).")
        hm, hs = stats(oh)
        im, is_ = stats(oi)
        for r in outliers(oh):
          gate_msgs.append(f"OUTLIER: {config} / {geom} octree host run {r} ({oh[r-1]:.3f} s, median "
                           f"{st.median(oh):.3f}).")
        oct_rows.append({"config": config, "geometry": geom, "n_runs": len(oh), "oct_init_mean": im, 
                         "oct_init_std": is_, "oct_host_mean": hm, "oct_host_std": hs, "oct_storage_MB": d.get("octreeStorageSize")})
    # -- patch runs per depth --------------------------------------
    for key, vals in d.items():
      m = re.fullmatch(r"rawPatchHostTimes_D(\d+)_Res", key)
      if not m:
        continue
      depth = int(m.group(1))
      init = d.get(f"rawPatchInitTimes_D{depth}_Res", [])
      vals, init, dropped = drop_contaminated(vals, init)
      for r in dropped:
        gate_msgs.append(f"EXCLUDED: {config}/{geom} D{depth} host run {r} (>{DROP_TOL:.0%} above median)")

      tot = [a + b for a, b in zip(init, vals)] if init else []
      hm, hs = stats(vals)
      im, is_ = stats(init) if init else (float("nan"),) * 2
      tm, ts = stats(tot) if tot else (float("nan"),) * 2
      for r in outliers(vals):
        gate_msgs.append(f"OUTLIER: {config}/{geom} D{depth} host run {r} ({vals[r-1]:.3f} s, median {st.median(vals):.3f})")
      seen.setdefault(config, set()).add((geom, depth))
      rows.append({"config": config, "geometry": geom, "Dmax": depth, "n_runs": len(vals), "init_mean": im, "init_std": is_,
                   "host_mean": hm, "host_std": hs, "total_mean": tm, "total_std": ts, "storage_MB": d.get(f"patchStorageSize_D{depth}")})

  # -- completeness ---------------------------------------------------
  expect = {"patch_full": {(g, dd) for g in GEOM_ORDER for dd in range(1, 6)} - {("Tet1331", 1)},
            "patch_noshort": {(g, 3) for g in GEOM_ORDER}, "standard_ascg": {(g, 3) for g in GEOM_ORDER}}
  print("=== [GATE] ===")
  total_expected = sum(len(v) for v in expect.values())
  total_seen = sum(len(seen.get(c, set()) & expect[c]) for c in expect)
  for c, want_set in expect.items():
    missing = sorted(want_set - seen.get(c, set()))
    if missing:
      print(f"MISSING in {c}: {missing}.")
  print(f"Configs present: {total_seen} / {total_expected}.")
  if gate_msgs:
    for msg in gate_msgs:
      print(msg)
  else:
    print(f"provenance echoes: ALL OK; no outliers above {OUTLIER_TOL:.0%}.")

  with open("campaign_stats.csv", "w", newline = "") as fh:
    w = csv.DictWriter(fh, fieldnames = rows[0].keys())
    w.writeheader()
    w.writerows(rows)
  with open("octree_stats.csv", "w", newline = "") as fh:
    w = csv.DictWriter(fh, fieldnames = oct_rows[0].keys())
    w.writeheader()
    w.writerows(oct_rows)
  print(f"\nWrote campaign_stats.csv ({len(rows)} rows), octree_stats.csv ({len(oct_rows)} rows).")

  # -- pooled octree per geometry (all configs' runs) -----------------
  pooled = {}
  for g in GEOM_ORDER:
    sub = [r for r in oct_rows if r["geometry"] == g]
    if sub:
      pooled[g] = (st.mean([r["oct_host_mean"] for r in sub]), st.mean([r["oct_host_std"] for r in sub]))

  by = {(r["config"], r["geometry"], r["Dmax"]): r for r in rows}

  # -- decomposition at D3 -------------------------------------------
  print("\n=== [DECOMP] (Dmax = 3) ===")
  print(f"{'geometry':>9} | {'search %':>10} | {'shortcut %':>11} | {'mem x':>6} | {'init x':>6}")
  dec_rows = []
  n_std_beats_oct, n_oct_beats_noshort = 0, 0
  for g in GEOM_ORDER:
    try:
      full = by[("patch_full", g, 3)]
      nosc = by[("patch_noshort", g, 3)]
      std_ = by[("standard_ascg", g, 3)]
    except KeyError:
      continue
    sc_pct, sc_err = ratio_with_err(nosc["host_mean"], nosc["host_std"], full["host_mean"], full["host_std"])
    se_pct, se_err = ratio_with_err(std_["host_mean"], std_["host_std"], full["host_mean"], full["host_std"])
    memx = (nosc["storage_MB"] / full["storage_MB"] if full["storage_MB"] else float("nan"))
    initx = (nosc["init_mean"] / full["init_mean"] if full["init_mean"] else float("nan"))
    oh = pooled.get(g, (float("nan"),))[0]
    if std_["host_mean"] < oh:
      n_std_beats_oct += 1
    if oh < nosc["host_mean"]:
      n_oct_beats_noshort += 1
    dec_rows.append({"geometry": g,
                     "t_standard": round(std_["host_mean"], 4),
                     "std_standard": round(std_["host_std"], 4),
                     "t_noshort": round(nosc["host_mean"], 4),
                     "std_noshort": round(nosc["host_std"], 4),
                     "t_full": round(full["host_mean"], 4),
                     "std_full": round(full["host_std"], 4),
                     "t_octree_pooled": round(oh, 4),
                     "shortcut_saving_pct": round(sc_pct, 2),
                     "shortcut_saving_err": round(sc_err, 2),
                     "search_saving_pct": round(se_pct, 2),
                     "search_saving_err": round(se_err, 2),
                     "mem_factor_noshort": round(memx, 2),
                     "init_factor_noshort": round(initx, 2)})
    print(f"{g:>9} | {se_pct:6.2f}±{se_err:4.2f} | {sc_pct:6.2f}±{sc_err:4.2f} | {memx:6.2f} | {initx:6.2f}")
  if not dec_rows:
    print("No ablation configs matched -- check the discovered config directory names printed above.")
    return
  with open("decomposition.csv", "w", newline = "") as fh:
    w = csv.DictWriter(fh, fieldnames = dec_rows[0].keys())
    w.writeheader()
    w.writerows(dec_rows)
  scs = [r["shortcut_saving_pct"] for r in dec_rows]
  ses = [r["search_saving_pct"] for r in dec_rows]
  mems = [r["mem_factor_noshort"] for r in dec_rows]
  inits = [r["init_factor_noshort"] for r in dec_rows]
  print(f"\nShortcut contribution range: {min(scs):.1f}-{max(scs):.1f}%.")
  print(f"Search contribution range: {min(ses):.1f}-{max(ses):.1f}%.")
  print(f"Memory factor range: {min(mems):.1f}-{max(mems):.1f}x.")
  print(f"Initialisation factor range: {min(inits):.1f}-{max(inits):.1f}x.")
  print(f"standard_ascg beats octree on {n_std_beats_oct} of {len(dec_rows)} geometries.")
  print(f"octree beats patch_noshort on {n_oct_beats_noshort} of {len(dec_rows)} geometries.")

  # -- headline ranges ------------------------------------------------
  print("\n=== [RANGES] (patch_full vs pooled same-session octree) ===")
  savings = {}
  d3_vs_d1 = {}
  for g in GEOM_ORDER:
    depths = [by[("patch_full", g, dd)] for dd in range(1, 6) if ("patch_full", g, dd) in by]
    if not depths or g not in pooled:
      continue
    best = min(depths, key = lambda r: r["host_mean"])
    oh, ohs = pooled[g]
    sv, sv_err = ratio_with_err(oh, ohs, best["host_mean"], best["host_std"])
    savings[g] = sv
    print(f"{g:>9}: optimal D{best['Dmax']} saving {sv:5.1f} ± {sv_err:3.1f}%.")
    d1 = by.get(("patch_full", g, 1))
    d3 = by.get(("patch_full", g, 3))
    if d1 and d3:
      dv, _ = ratio_with_err(d1["host_mean"], d1["host_std"], d3["host_mean"], d3["host_std"])
      d3_vs_d1[g] = dv
  print(f"\nABSTRACT RANGE: {min(savings.values()):.0f}-{max(savings.values()):.0f}% (over {len(savings)} geometries).")
  wins = {g: v for g, v in d3_vs_d1.items() if v > 0}
  print(f"D3 vs D1 host: wins on {len(wins)} of {len(d3_vs_d1)} (best +{max(d3_vs_d1.values()):.1f}%, worst {min(d3_vs_d1.values()):.1f}%).")

if __name__ == "__main__":
  main()