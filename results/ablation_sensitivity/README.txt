CGC-GOCNMF revision experiments
MATLAB R2009a
================================

Purpose
-------
This package adds the two experiments requested for the revised manuscript:

1. Six-dataset ablation:
   BASE, LOCAL, HARD, and SHRINK on PIE, YaleB, COIL20, COIL100,
   Optdigits, and MNIST.

2. Shrinkage-rule sensitivity:
   theta_kappa = max(0,(Delta-kappa*s)/Delta) for Delta>0,
   with kappa in {0,0.5,1,1.5} and F in {3,5,10}.
   The performance sensitivity is run on three pre-specified regimes:
   PIE (fallback), YaleB (borderline), and COIL100 (clear benefit).

Frozen protocol
---------------
- 10% labeled samples per class, at least two.
- Unlabeled-sample evaluation only.
- Direct row-wise argmax; no K-means.
- 50 factorization iterations.
- Same labeled split and initialization in every paired comparison.
- GOCNMF graph parameters from the manuscript.
- MATLAB R2009a-compatible syntax.
- Checkpoint after every completed dataset-seed run.

Required files
--------------
Place the three .m files in the same folder. The dataset files may be in
that folder or its subfolders. Accepted names include:

CMU_PIE_fac.mat / CMU_PIE.mat / PIE.mat
YaleB.mat
COIL20_Obj.mat
COIL100_Obj.mat / COIL100_Obj(1).mat
Optdigits_Han.mat
MNIST_Han.mat

Run order
---------
A. Three-seed integrity tests

clear functions;
rehash path;
run_CGC_GOCNMF_ablation6_R2009a(3,pwd,pwd);
run_CGC_GOCNMF_sensitivity3_R2009a(3,pwd,pwd);

Confirm in the decision files:
THREE_SEED_ABLATION_INTEGRITY_PASS=1
THREE_SEED_SENSITIVITY_INTEGRITY_PASS=1

B. Formal 20-seed experiments

run_CGC_GOCNMF_ablation6_R2009a(20,pwd,pwd);
run_CGC_GOCNMF_sensitivity3_R2009a(20,pwd,pwd);

Confirm:
TWENTY_SEED_ABLATION_COMPLETE=1
TWENTY_SEED_SENSITIVITY_COMPLETE=1

Checkpointing
-------------
Do not delete the checkpoint MAT files after interruption. Re-running the
same command resumes completed dataset-seed runs and regenerates CSV files.

Upload after completion
-----------------------
CGC_GOCNMF_ablation6_20seed_raw.csv
CGC_GOCNMF_ablation6_20seed_summary.csv
CGC_GOCNMF_ablation6_20seed_decision.txt
CGC_GOCNMF_sensitivity3_20seed_raw.csv
CGC_GOCNMF_sensitivity3_20seed_summary.csv
CGC_GOCNMF_sensitivity3_20seed_decision.txt

The programs also generate LaTeX tables automatically.
