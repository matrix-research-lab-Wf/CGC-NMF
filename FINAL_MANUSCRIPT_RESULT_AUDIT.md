# Final manuscript-to-result audit

## Resolution update (2026-10-09)

The public-release manuscript now corrects the SNMFWLP stopping statement. A
base-MATLAB R2019a script and generated CSV archive the exact CGC-GNMFLD
Wilcoxon--Holm results, and a fixed-seed 10,000-resample CSV archives the paired
bootstrap intervals. The manuscript table was synchronized to that fixed-seed
output and the PDF was rebuilt and visually inspected page by page. The detailed
findings below describe the pre-correction audit and are retained as an audit
trail.

## Audited files

- Source: `manuscript/Evidence-guided conservative graph calibration.tex`
- Rendered article: `manuscript/Evidence-guided conservative graph calibration.pdf`
- PDF length: 16 pages

The PDF contains the common-stopping table and the same SNMFWLP stopping-text
issue found in the LaTeX source. Thus, the discrepancy is present in the rendered
submission artifact, not only in an unused source file.

## Overall result

Most displayed numerical tables are traceable and agree with located raw or
summary files at the displayed precision. One factual statement is wrong, one
statistical-output provenance mismatch remains, and the release package still
needs explicit provenance notes for experiments run under different MATLAB
versions.

## Itemized audit

| Manuscript item | Verdict | Evidence |
|---|---|---|
| Method-positioning table | NOT AUDITED | Requires a paper-by-paper citation and feature audit rather than experimental CSV checking. |
| Dataset characteristics | PASS | Six datasets' n, m, c, labeled counts, p, and alpha match the main 20-seed raw/parameter files. |
| Four-method fixed-budget comparison | PASS | GNMFLD/GOCNMF/CGC values match `three_direct_methods_20seed_summary.csv`; ERDNMF values match `results/erdnmf_main/ERDNMF_final_six_datasets_20seed_summary.csv`. |
| Common-stopping table | PASS | All 48 displayed ACC/NMI mean and sample-SD values match the regenerated six-dataset three-seed summaries after percentage conversion and two-decimal rounding. |
| Common-stopping narrative | FAIL | The manuscript says SNMFWLP reached the threshold in only a subset of runs. The regenerated raw files report `converged=1` for 18/18 runs. Replace this sentence. |
| Common-stop objective diagnostics | PASS | The statement that SNMFWLP's monitored objective was not monotone on five datasets matches the audit: PIE 3/3 passes; each other dataset 0/3. |
| Fixed-budget Holm table | PASS WITH REGENERATION ADVISED | Main trial-level data and ERDNMF paired data are present and support the reported direction/significance. A single generated all-baseline Holm CSV should still be archived. |
| Paired-effect table | PARTIAL | Mean changes and W/T/L agree with the paired CSV. The displayed 10,000-resample confidence intervals have no separately archived bootstrap-output file in the package. |
| External-baseline parameters | PASS | GNMFLD parameters match the main parameter CSV; ERDNMF coefficients match the recovered final summary and timing files. |
| CGC-GNMFLD transfer table | NUMBERS PASS; P-VALUE ARTIFACT CONFLICT | Means, SDs, theta, changes, and W/T/L match `CGC_GNMFLD_transfer_raw_FIXED.csv`. The manuscript reports exact Wilcoxon-Holm `2.29e-5`, which is mathematically consistent with 20/20 positive paired differences across 12 tests. The archived transfer summary instead prints `0.00106289`; regenerate and archive the exact-test statistics file so code and paper agree. |
| NEU-CLS raw-pixel performance | PASS | Four-method means/SDs match the NEU and ERDNMF 20-seed summaries. |
| NEU-CLS calibration | PASS | Theta means/SDs and fallback counts recompute exactly from raw CSV. |
| NEU-CLS confusion figure | PASS AT FILE LEVEL | Seven manuscript PNGs were found and copied. Count and normalized confusion CSVs are present. A pixel-level regeneration comparison has not yet been performed. |
| NEU-CLS LBP table | PASS | The recovered 20-seed LBP summary, paired tests, decision file, and raw file match every displayed mean/SD, theta, fallback, and significance marker. |
| Six-dataset ablation | PASS | All displayed theta, fallback, delta, and loss-count entries match the 20-seed ablation summary. |
| Sensitivity table | PASS | Displayed F/kappa rows match the sensitivity summary at three-decimal precision. |
| Graph-stress table | PASS | Correlation, theta, delta, metric changes, and loss counts match the stress summary. |
| Stress figures | PASS AT FILE LEVEL | `fig-004.png`--`fig-006.png` were found and copied; numerical source data are present. Pixel-level regeneration has not yet been performed. |
| Timing table | PASS NUMERICALLY | Three graph methods match the timing summary; ERDNMF matches the recovered dedicated timing summary. The ERDNMF timing audit records MATLAB R2009a, whereas the recovered common-stopping work uses R2019a. This environment difference must be disclosed. |

## Mandatory manuscript correction

Current sentence:

> SNMFWLP reached the stopping threshold in only a subset of runs, and its
> monitored unconstrained objective was not monotone on five datasets.

Evidence-consistent replacement:

> SNMFWLP reached the common stopping threshold in all 18 formal runs, while its
> monitored unconstrained objective was not monotone on five of the six datasets;
> these diagnostics are reported rather than interpreted as a stationary-point
> convergence guarantee.

## Files newly recovered into the package

- Exact manuscript PDF and LaTeX source
- `fig-000.png`--`fig-006.png`
- Six-dataset ERDNMF 20-seed raw and summary files
- CGC-GNMFLD transfer raw and summary files
- NEU-CLS LBP 20-seed raw, summary, paired-test, decision, and confusion files
- ERDNMF 20-seed timing raw, summary, decision, checkpoint, and table column
- MATLAB scripts for the recovered LBP and CGC-GNMFLD experiments

## Remaining release actions

1. Correct the false SNMFWLP stopping sentence in both source and PDF.
2. Generate one exact-Wilcoxon/Holm CSV for the CGC-GNMFLD transfer table and
   remove or supersede the conflicting `0.00106289` audit output.
3. Archive the 10,000-resample paired-effect bootstrap intervals with seeds and
   resampling code.
4. State explicitly that the ERDNMF timing measurements were produced under
   MATLAB R2009a; do not describe every experiment as R2019a-only.
5. Perform clean-clone execution and figure-regeneration checks before publishing
   the repository URL.
