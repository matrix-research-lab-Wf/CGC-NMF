# Upload inventory

Prepared for: <https://github.com/matrix-research-lab-Wf/CGC-NMF>

| Category | Location | Included content | Release state |
|---|---|---|---|
| MATLAB source | `code/`, `run_smoke_test.m` | Core routines, comparisons, ablations, stress tests, timing, NEU-CLS, recovered protocols | Pending third-party license audit |
| Experimental evidence | `results/` | Trial-level CSV, summaries, decisions, logs, checkpoints, and LaTeX table fragments | Included; known audit exceptions documented |
| Manuscript | `manuscript/` | Current LaTeX source and 16-page rendered PDF | Included; one mandatory factual correction remains |
| Figures | `figures/` | Seven figures used by the manuscript | Included; pixel-level regeneration comparison pending |
| Datasets | `data/` | Six benchmark matrices and two NEU-CLS matrices | Local staging only; redistribution rights unresolved |
| Reproducibility records | root Markdown files | Provenance, recovery notes, result audit, release checklist, license scope | Included |
| Integrity record | `MANIFEST_SHA256.csv` | Relative path, byte count, and SHA-256 for each packaged file except the manifest itself | Included |
| Release utilities | `tools/` | Recovery helpers, release builder, and pre-upload verifier | Included |

The local staging directory contains the complete collected package. This inventory
does not itself establish ownership, dataset redistribution permission, numerical
correctness beyond the cited audit, or successful execution from a clean clone.
