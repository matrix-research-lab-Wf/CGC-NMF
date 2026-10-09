Minimal three-dataset ablation for conservative graph-calibrated GOCNMF
MATLAB R2009a

Required data files in the same directory:
  COIL100_Obj(1).mat  (also accepts COIL100_Obj.mat or COIL100.mat)
  PIE.mat
  YaleB.mat           (also accepts YaleB(1).mat)

Compared graph variants:
  BASE   = original binary p-nearest-neighbor graph
  LOCAL  = self-tuning locally weighted graph
  HARD   = LOCAL when Delta_cv>0, otherwise BASE
  SHRINK = (1-theta)BASE + theta LOCAL

HARD reuses the BASE or LOCAL result and does not require another NMF run.

Step 1: three-seed integrity/mechanism check
  clear functions;
  rehash path;
  which run_GCGOCNMF_min_ablation3_R2009a -all
  run_GCGOCNMF_min_ablation3_R2009a(3, pwd, pwd);

Required output:
  CV_FINITE_PASS=1
  SHRINK_SAFETY_PASS=1
  ABLATION_MECHANISM_PASS=1

Do not run 20 seeds until all three values equal 1.

Step 2: twenty seeds
  run_GCGOCNMF_min_ablation3_R2009a(20, pwd, pwd);

Generated files:
  GCGOCNMF_min_ablation3_3seed_raw.csv
  GCGOCNMF_min_ablation3_3seed_summary.csv
  GCGOCNMF_min_ablation3_3seed_decision.txt
  GCGOCNMF_min_ablation3_3seed_results.mat

The formulas, graph construction, 10% label rate, five-fold correction,
p, alpha, initialization seeds, and 50 iterations are frozen.
