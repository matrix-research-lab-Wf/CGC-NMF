# CGC-NMF reproducibility package

[![MATLAB R2019a](https://img.shields.io/badge/MATLAB-R2019a-orange.svg)](https://www.mathworks.com/products/matlab.html)

Reproducibility materials for **Evidence-guided conservative graph calibration
for direct semi-supervised clustering** by Deshu Sun and Feng Wang.

This directory is a cleaned local package assembled from the authors' experiment
files for the manuscript on evidence-guided conservative graph calibration.

## Important status

This repository contains the public-release package. Dataset licensing was
checked before release: the attributed CC BY 4.0 Optdigits matrix is included,
whereas matrices without verified redistribution permission are excluded and
documented in `DATA_PROVENANCE.md`. The formerly
missing unified common-stopping and SNMFWLP scripts were reconstructed from the
original Codex session patch history and rerun in MATLAB R2019a. Six-dataset,
three-seed outputs are available under `results/recovered_six_dataset_3seed`.
See `RECOVERY_NOTES.md` for provenance, validation evidence, and limitations.

## Environment

- Microsoft Windows
- MATLAB R2019a
- No GPU is assumed
- The main six-dataset scripts use basic MATLAB functionality

Toolbox requirements for the NEU-CLS feature-preparation scripts must still be
checked on a clean MATLAB R2019a installation.

## Directory structure

- `code/core`: evaluation utilities and baseline routines
- `code/main_comparison`: paired GNMFLD/GOCNMF/CGC-GOCNMF experiment
- `code/erdnmf`: ERDNMF and terminal-validation experiments
- `code/ablation_sensitivity`: ablation and sensitivity experiments
- `code/graph_stress`: controlled graph-corruption experiments
- `code/timing`: prescribed-protocol timing experiments
- `code/neu_cls`: NEU-CLS experiments and figure scripts
- `code/data_preparation`: official download, deterministic conversion, and validation programs
- `code/legacy_runs`: older dataset-specific experiment scripts retained for audit
- `code/recovered_latest`: recovered common-stopping and SNMFWLP scripts
- `data/main_six`: redistributable Optdigits matrix; other expected paths are documented
- `data/neu_cls`: NEU-CLS preparation scripts consume author-downloaded images locally
- `results`: archived raw outputs, summaries, checkpoints, and LaTeX tables
- `results/recovered_six_dataset_3seed`: regenerated six-dataset, three-seed outputs
- `results/erdnmf_main`: recovered six-dataset ERDNMF evidence
- `results/cgc_gnmfld_transfer`: recovered second-parent transfer evidence
- `results/neu_cls_lbp`: recovered 20-seed LBP robustness evidence
- `results/erdnmf_timing`: recovered ERDNMF timing evidence
- `figures`: the seven figures used by the current manuscript

The current table-by-table evidence audit is `FINAL_MANUSCRIPT_RESULT_AUDIT.md`.
It must be read before citing or redistributing the archived results.

## Public smoke test

1. Start MATLAB R2019a.
2. Change the current folder to this repository root.
3. Run:

```matlab
run_public_smoke_test
```

The script runs the common-stopping comparison on the included Optdigits matrix.
New outputs are written to `results/reproduced_public_smoke`. The complete
six-dataset smoke test, `run_smoke_test`, requires the locally acquired matrices
listed in `data/README.md`.

To acquire and prepare the other datasets without unofficial mirrors, see
`data/README.md` and run `prepare_all_official_data_R2019a`. Exact
pixel-and-label reproduction is verified for MNIST_Han and COIL20. The
repository explicitly reports that legacy preprocessing identity is not established for
COIL100, PIE, and YaleB rather than presenting structural compatibility as an
exact reproduction.

## Full archived protocols

The main scripts accept `runs = 1`, `3`, or `20`. For example:

```matlab
rootDir = pwd;
dataDir = fullfile(rootDir,'data','main_six');
outputDir = fullfile(rootDir,'results','reproduced_20seed');
addpath(genpath(fullfile(rootDir,'code')));
run_three_direct_methods_R2009a(20,dataDir,outputDir);
run_ERDNMF_six_datasets_R2009a(20,dataDir,outputDir);
run_CGC_GOCNMF_ablation6_R2009a(20,dataDir,outputDir);
run_CGC_GOCNMF_sensitivity3_R2009a(20,dataDir,outputDir);
```

Random seeds and fixed parameters are recorded inside the experiment scripts.
The archived CSV files in `results` provide the previously generated trial-level
and summary outputs.

## Data and redistribution

`DATA_PROVENANCE.md` records authoritative sources, redistribution decisions,
matrix shapes, and SHA-256 values. Do not commit locally obtained PIE, YaleB,
COIL, MNIST, or NEU-CLS matrices; their expected names are blocked by `.gitignore`.

## Release gates

See `PRE_UPLOAD_CHECKLIST.md`. The repository's MIT license applies only to
original author-owned code; dataset and third-party boundaries are described in
`LICENSE_SCOPE.md`, `DATA_PROVENANCE.md`, and `THIRD_PARTY_NOTICES.md`.

## Recovered common-stopping experiments

The recovered scripts use formal seeds 20260617--20260619 and direct row-wise
argmax assignment. All six three-method decision files report
`THREE_SEED_INTEGRITY_PASS=1`. All SNMFWLP runs reached the stated residual-based
stopping rule and returned finite outputs. However, the recorded SNMFWLP objective
was nonincreasing for PIE only; this limitation must remain visible when reporting
the baseline. COIL100 requires substantial memory, so the reconstruction objective
is evaluated in sample blocks and the SNMFWLP run is launched in a separate fresh
MATLAB process.

## Code language

The packaged MATLAB comments were scanned for non-ASCII text. The Chinese comments
found in `Accuracy.m` were translated in this packaged copy; the original source tree
was not modified.
