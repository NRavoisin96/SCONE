# JCP benchmarks — Patch-Search host-element determination

Scripts, inputs and measurement data supporting:

> D. Kim, N. Ravoisin, G. T. Parks, *Efficient Host-Element Determination in
> Convex Polyhedral Meshes Using an Adaptive Patch-Search Algorithm*,
> Journal of Computational Physics.

Everything in this directory is specific to the paper; the rest of the
repository is SCONE itself. Section, table and figure numbers below refer to
the manuscript.

---

## 1. Build configurations

Three (mutually exclusive) builds are used.

| Purpose | CMake invocation |
| --- | --- |
| Timing and memory measurements (Sections 4.3.1–4.3.5) | `cmake .. ` |
| Mapping-frequency statistics (Section 4.2) | `cmake .. -DPATCH_SEARCH_STATS=ON -DOPENMP=OFF` |
| Host-element validation against brute force (Section 4.1) | `cmake .. -DVALIDATE=ON` |

`PATCH_SEARCH_STATS` instruments the query path with mapping-trigger counters. This build must never
be used for timing runs, as the counters are incremented in hot code paths. Furthermore, these counters
are **not** thread-safe; as such, attempting to use this build with OpenMP enabled throws an error.

`VALIDATE` compares every host-element index returned by the acceleration
structure against a brute-force search over all elements (Eq. 9), and throws an error on the first
disagreement between the two methods. This reproduces the validation of Section 4.1.

All measurements reported in the paper were obtained on a single CPU core (Table 1).

---

## 2. Layout

```
JCPBenchmarks/
  Scripts/       sweep drivers and parsers
  Templates/     SCONE input templates with placeholder tokens
  Meshes/        the mesh geometries used in the study
  Results/       raw SCONE outputs and the processed CSVs
```

**Scripts resolve their input and output paths relative to the working**
**directory and must be run from `JCPBenchmarks/`, not from `Scripts/`.**

The SCONE executable is expected at `../Build/scone.out`. Set the `SCONE_BIN`
environment variable to point elsewhere.

---

## 3. Reproducing the measurements

### Host-determination, ablation and initialisation studies

```
python3 Scripts/runPatchSearchSweeps.py --dry-run  # render inputs, run nothing
python3 Scripts/runPatchSearchSweeps.py            # full sweep
```

One SCONE run per (configuration, geometry). Sweeps are long—run under
`tmux` or `screen`. Completed runs hold a `DONE` marker and are skipped on driver re-invocation, 
so that an interrupted sweep can easily be resumed.

The five configurations below map onto the manuscript as follows:

| Configuration | Description | Manuscript |
| --- | --- | --- |
| `patch_full` | AMLG Patch-Search, all optimisations enabled | Configuration (a), Section 4.3.5; main benchmark of Section 4.3.2 |
| `patch_noshort` | Single-face optimisation disabled, both as a refinement-termination criterion and as a query-time branch | Configuration (b), Section 4.3.5 |
| `standard_ascg` | AMLG hierarchy retained; pseudo-angular sector searches replaced by conventional point-in-element tests | Configuration (c), Section 4.3.5 |
| `patch_optimised_init` | Initialisation built with cached separating-axis quantities and the centroid-based inclusion test | Section 4.3.4 |
| `patch_naive_init` | Initialisation built with the naive procedures | Section 4.3.4 |

`standard_ascg` is a legacy identifier retained for consistency with
the archived result directories, and denotes configuration (c)—the adaptive
multi-layered hierarchy without pseudo-angular sector searches, **not**
the uniform auxiliary structured Cartesian grid method of the original Patch-Search
formulation.

The two initialisation configurations use a reduced query population
(one million rather than one billion), since only initialisation time is
measured from them. `patch_naive_init` produces a byte-identical acceleration
structure to `patch_optimised_init`; only the build cost differs, which is checked
by the parser.

`Tet1331` is swept from `Dmax = 2` only as the uniform, single-layer grid exhausts
memory on this geometry (Section 4.3.3).

### Octree sensitivity study

```
python3 Scripts/runOctreeOptimisation.py
```

One SCONE run per geometry; the octree optimisation package loops over the full
depth × maximum-face-count × run grid internally, which establishes the
`{Dmax = 6, Nfaces,max = 1}` parameter pair used for the octree baseline
throughout the paper (Section 4.3.1).

### Mapping-frequency statistics

Requires the `PATCH_SEARCH_STATS` build. Run `runPatchSearchSweeps.py` restricted
to `patch_full`, then parse with `parseMappingStatisticsResults.py`.

---

## 4. Parsing

Each parser reads raw SCONE outputs and writes the CSVs used by the paper's tables and figures.

| Script | Reads | Writes |
| --- | --- | --- |
| `parseMainCampaignResults.py` | `Results/MainCampaign`, `Results/PerformanceDecompositionStudy` | `campaign_stats.csv`, `octree_stats.csv`, `decomposition.csv` |
| `parseInitialisationStudyResults.py` | `Results/InitialisationStudy` | `init_study_long.csv`, `init_factors.csv` |
| `parseMappingStatisticsResults.py` | `Results/MappingStatistics/patch_full` | `mapping_stats_long.csv`, `mapping_stats_table.csv` |
| `parseOctreeOptimisationResults.py` | `Results/OctreeOptimisation` | `octree_sensitivity_long.csv`, `octree_sensitivity_stats.csv` |

`parseMainCampaignResults.py` takes both roots as arguments:

```
python3 Scripts/parseMainCampaignResults.py Results/MainCampaign Results/PerformanceDecompositionStudy
```

Configurations are identified from the values of the `patchType` and `singleFaceShortcut` 
fields written into each SCONE output rather than from directory names, so that a misfiled 
run is reported rather than silently mislabelled.

Each parser also prints a consistency report: completeness against the expected
run matrix, excluded runs, as well as study-specific checks (termination
closure and volume closure for the mapping statistics; identical storage sizes between
the naive and optimised initialisation configurations).

---

## 5. Data files

`Results/` holds the raw SCONE outputs alongside the processed CSVs.

| File | Feeds |
| --- | --- |
| `campaign_stats.csv` | Tables H.12–L.16; Figures 9, 10, 11 |
| `octree_stats.csv` | Table G.11 |
| `decomposition.csv` | Figure 13; Table R.22 |
| `init_factors.csv` | Figure 12; Tables M.17–Q.21 |
| `init_study_long.csv` | per-run initialisation times |
| `mapping_stats_table.csv` | Table 4; Appendices B–F |
| `mapping_stats_long.csv` | per-run trigger counts |
| `octree_sensitivity_stats.csv` | Figures 7, 8 |
| `octree_sensitivity_long.csv` | per-run octree sensitivity data |

Figures were produced from these CSVs with separate plotting scripts, which are
not included here.

Per-run timings for the main campaign are listed in the `rawPatchHostTimes_D*_Res`, 
`rawPatchInitTimes_D*_Res` and `rawOctreeHostTimes_Res` arrays inside each SCONE 
output file under `Results/MainCampaign` and `Results/PerformanceDecompositionStudy`. 
`campaign_stats.csv` contains the aggregated statistics derived from them.

### Measurement protocol

Timing results come from a single measurement campaign, with 20 independent
runs per configuration and a shared set of RNG seeds so that comparisons
between configurations are paired. Runs whose host-determination time exceeded
the per-configuration median by more than 15% were excluded as background
activity contamination; the parsers apply this criterion and report every
exclusion. Reported statistics are therefore computed over at least 19 runs per
configuration.

---

## 6. Input templates

`Templates/` holds ordinary SCONE input files in which configuration-dependent 
entries have been replaced by placeholder tokens (`@POPULATION@`, `@PATCHTYPE@`, 
`@NAIVE@`, `@SHORTCUT@`, `@DEPTHS@`, `@SEEDS@`, `@OUTPUT@`, and `@NMAXFACES@` / 
`@NRUNS@` for the octree sweep). The drivers substitute these and write the 
rendered input alongside each run's output. A missing token is treated as an error.

Everything else in the templates—geometry and nuclear data blocks, octree
settings—is left as written.