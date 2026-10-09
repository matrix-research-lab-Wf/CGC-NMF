# Manuscript-to-result audit

> Historical audit of the older August manuscript. Superseded for the current
> submission by `FINAL_MANUSCRIPT_RESULT_AUDIT.md`.

## Audit target and scope

Audited manuscript source:

`D:\codex\CGC_GOCNMF_local_manuscript_archive_20261008\CGC_GOCNMF.tex`

This source is dated 4 August 2026 and uses ERDNMF as the fourth baseline. It is
not the later ESWA/SNMFWLP manuscript reconstructed in the September session.
No physical copy of that later manuscript source or PDF is currently present in
the searched package or `D:\codex` tree. Consequently, the recovered SNMFWLP
three-seed results cannot yet be checked against their intended manuscript table.

Status terms:

- **PASS**: manuscript values agree with an included raw/summary artifact at the
  displayed precision.
- **PARTIAL**: only part of the table has a traceable artifact.
- **UNVERIFIED**: the result artifact required for numerical checking is absent.
- **PROTOCOL MISMATCH**: the available result uses a different experimental
  contract and must not be substituted.

## Table-by-table findings

| Manuscript item | Status | Evidence and finding |
|---|---|---|
| Table 1, representative-method comparison | UNVERIFIED | This is a literature-comparison table. Its citations and qualitative attributes require a separate source audit; no numerical experiment file can verify it. |
| Dataset statistics and graph parameters | PASS | All six datasets' sample counts, feature counts, class counts, labeled counts, graph-neighbor values, and GOC/CGC coefficients match `three_direct_methods_20seed_raw.csv` and `three_direct_methods_20seed_parameters.csv`. |
| Four-method clustering table | PARTIAL | Every displayed ACC/NMI mean and sample SD for GNMFLD, GOCNMF, and CGC-GOCNMF matches `three_direct_methods_20seed_summary.csv` after percentage conversion and two-decimal rounding. The six-dataset ERDNMF raw/summary files are absent, so all ERDNMF cells remain unverified. |
| Holm-adjusted p-value table | PARTIAL | Trial-level data exist for comparisons with GOCNMF and GNMFLD, but no generated six-dataset significance file is included. Comparisons with ERDNMF cannot be recomputed because the six-dataset ERDNMF trial-level file is missing. |
| External-baseline parameter table | PARTIAL | All GNMFLD parameters and the two explicitly transferred settings match the included parameter CSV. The stated six-dataset ERDNMF settings/provenance do not have a corresponding included selection/audit artifact. |
| NEU-CLS performance table | PASS | GNMFLD, GOCNMF, and CGC-GOCNMF entries match `NEU_application_20seed_summary.csv`; ERDNMF entries match `ERDNMF_NEU_formal_20seed_summary.csv`. Values agree after percentage conversion and two-decimal rounding. |
| NEU-CLS calibration table | PASS | Mean theta, sample SD, and fallback counts recomputed from `NEU_application_20seed_raw.csv` agree with 0.296/0.287/7, 0.646/0.183/0, and 0.790/0.063/0. |
| NEU-CLS confusion-matrix figure | PARTIAL | Aggregated count and normalized CSV files exist for all four methods. The manuscript's `fig-000.png`--`fig-003.png` files are absent from the reproducibility package, so rendered-panel identity and common color scaling cannot be checked. |
| LBP representation table | UNVERIFIED | The package contains the LBP feature matrix and an LBP experiment script, but no LBP raw, summary, paired-statistics, parameter, decision, or confusion output files. None of the displayed LBP values is currently reproducible from an archived result artifact. |
| Six-dataset ablation table | PASS | Theta, fallbacks, LOCAL/HARD/SHRINK changes, and loss counts match `CGC_GOCNMF_ablation6_20seed_summary.csv` at the displayed precision. |
| Calibration-sensitivity table | PASS | All displayed F=5 kappa rows and kappa=1 fold rows match `CGC_GOCNMF_sensitivity3_20seed_summary.csv` at three-decimal precision. |
| COIL100 graph-stress table | PASS | Edge correlations, theta, out-of-fold differences, ACC/NMI changes, and loss counts match `GCGOCNMF_graph_stress_20seed_summary.csv`. |
| Graph-stress figures | UNVERIFIED | The numerical summary exists, but the manuscript image files `fig-004.png`--`fig-006.png` are absent from the package. |
| Prescribed-protocol timing table | PARTIAL | GNMFLD, GOCNMF, CGC-GOCNMF, calibration time, overhead, and calibration share agree with `three_methods_end_to_end_timing_20seed_summary.csv`. The six ERDNMF timing means/SDs have no included timing raw/summary artifact. |

## Recovered later experiments

The recovered common-stopping outputs use three seeds and convergence-based
stopping. They are not replacements for the manuscript's 20-trial,
fixed-iteration table. The manuscript explicitly states 200 iterations for
GNMFLD/ERDNMF and 50 iterations for GOCNMF/CGC-GOCNMF, whereas the recovered
common-stopping experiment permits up to 1000 iterations. Mixing the two would
invalidate the table caption, significance tests, and protocol description.

The recovered SNMFWLP outputs have the following audit result:

- all 18 runs returned finite outputs and reached the residual stopping rule;
- objective nonincrease passed for PIE 3/3;
- objective nonincrease failed for YaleB, COIL20, Optdigits, MNIST, and COIL100
  (0/3 for each dataset).

Therefore, the later manuscript must not claim monotonic objective behavior for
these SNMFWLP runs unless the objective implementation or update equations are
corrected and the experiments are rerun.

## Blocking gaps before public submission

1. Recover the exact later manuscript PDF and LaTeX source that contains the
   SNMFWLP/common-stopping table.
2. Restore or rerun the six-dataset ERDNMF raw, summary, parameter-selection,
   significance, and timing artifacts if the ERDNMF manuscript is retained.
3. Rerun and archive the complete NEU-CLS LBP experiment outputs.
4. Restore or regenerate all seven manuscript figure files from named scripts and
   archived numerical inputs.
5. Recompute and archive the six-dataset Holm-adjusted significance table from the
   exact trial-level files used by the manuscript.
6. Do not merge fixed-iteration 20-trial values with common-stopping three-seed
   values in one performance table.

## Current decision

The locally available August manuscript is numerically well supported for the
three principal methods, NEU-CLS raw-pixel experiment, ablation, sensitivity,
and graph-stress results. It is not yet a complete PAA reproducibility package
because the main ERDNMF evidence, LBP outputs, several figure files, and exact
later manuscript version are missing.
