# Upload inventory

Prepared for: <https://github.com/matrix-research-lab-Wf/CGC-NMF>

| Category | Location | Included content | Release state |
|---|---|---|---|
| MATLAB source | `code/`, `run_smoke_test.m` | Core routines, comparisons, ablations, stress tests, timing, official data preparation, NEU-CLS, recovered protocols | Included; provenance boundaries documented |
| Experimental evidence | `results/` | Trial-level CSV, summaries, decisions, logs, checkpoints, and LaTeX table fragments | Included; known audit exceptions documented |
| Figures | `figures/` | Seven figures used by the manuscript | Included; pixel-level regeneration comparison pending |
| Datasets | `data/` | Redistributable Optdigits matrix plus official acquisition/preparation instructions | Optdigits included; restricted matrices excluded from Git |
| Reproducibility records | root Markdown files | Provenance, recovery notes, result audit, release checklist, license scope | Included |
| Integrity record | `MANIFEST_SHA256.csv` | Relative path, byte count, and SHA-256 for each packaged file except the manifest itself | Included |
| Release utilities | `tools/` | Recovery helpers, release builder, pre-upload verifier, and data-preparation release gate | Included |

Restricted author-side matrices are held outside the Git checkout. This inventory
does not itself establish dataset redistribution permission or numerical correctness
beyond the cited audit and explicit regeneration checks.
