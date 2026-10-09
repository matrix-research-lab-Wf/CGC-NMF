NEU-CLS engineering application -- Stage 1
===========================================

Programs
--------
1. prepare_NEU_CLS_R2009a.m
2. run_NEU_application_R2009a.m
3. summarize_NEU_application_R2009a.m

Step 1: prepare NEU-CLS
-----------------------
Extract the image-only NEU-CLS dataset, then run:

prepare_NEU_CLS_R2009a( ...
    'D:\datasets\NEU-CLS', ...
    fullfile(pwd,'NEU_CLS_32x32.mat'));

Required output:
NEU_DATA_PREPARATION_PASS=1

The loader supports standard filename prefixes Cr, In, Pa, PS, RS, and Sc,
or class subfolders. It checks for exactly 1,800 images and 300 per class.
Images are converted to 32x32 grayscale raw intensities. No deep features or
Image Processing Toolbox resize function is used.

Step 2: smoke test
------------------
run_NEU_application_R2009a(1,pwd,pwd);

The decision file must contain:
NEU_ONE_SEED_SMOKE_PASS=1

Step 3: integrity gate
----------------------
Delete the 1-seed checkpoint only if you want a clean directory, then run:

run_NEU_application_R2009a(3,pwd,pwd);

The decision file must contain:
NEU_THREE_SEED_INTEGRITY_PASS=1

Step 4: final experiment
------------------------
run_NEU_application_R2009a(20,pwd,pwd);

The decision file must contain:
NEU_TWENTY_SEED_APPLICATION_COMPLETE=1

Locked Stage-1 protocol
-----------------------
- Label fractions: 5%, 10%, 20% per class
- Seeds: 20260731 onward
- Primary evaluation: unlabeled samples only
- Metrics: ACC, NMI, ARI
- Direct row-wise argmax; no K-means
- GOCNMF and CGC-GOCNMF: p=5, alpha=10, 50 iterations
- GNMFLD: p=5, alpha=1e4, beta=10, 200 iterations
- Five-fold out-of-fold graph calibration
- Outputs: theta, Delta, standard error, fallback count, confusion matrices,
  objective audit, raw and summarized CSV files

Editorial limitation
--------------------
These are pre-specified pilot parameters, not claimed NEU-optimal. No
unlabeled ground-truth label is used for tuning. Do not state performance
superiority until the 20-seed run has completed and passed all audit gates.
