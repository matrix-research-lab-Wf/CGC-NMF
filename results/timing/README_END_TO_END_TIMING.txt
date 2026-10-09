Dedicated end-to-end timing for three direct methods
MATLAB R2009a

Purpose
-------
Replace the incomplete run_seconds values from the accuracy experiment with
a formal, independent timing experiment.

Included in every reported total
--------------------------------
1. Method-specific graph construction.
2. Complete five-fold out-of-fold Brier calibration for CGC-GOCNMF.
3. Matrix factorization with the locked iteration counts.

Excluded
--------
1. Disk loading.
2. Common column normalization.
3. Post-training ACC/NMI evaluation.

Fairness controls
-----------------
- Every method independently constructs its graph.
- No graph is shared between methods.
- CGC-GOCNMF never reuses GOCNMF factors or outputs.
- CGC-GOCNMF runs its own factorization even when theta=0.
- The same labeled subset and initial factors are supplied to all methods.
- Method execution order is randomized reproducibly for every dataset-seed.
- One full untimed warm-up is run for every method and dataset.
- Console output is suppressed inside timed graph construction.
- All methods run on the same MATLAB session and hardware.

Run sequence
------------
1. Smoke test:
   clear functions;
   rehash path;
   which run_three_methods_end_to_end_timing_R2009a -all
   run_three_methods_end_to_end_timing_R2009a(1,pwd,pwd);

   Required:
   ONE_SEED_TIMING_SMOKE_PASS=1

2. Three-seed integrity:
   run_three_methods_end_to_end_timing_R2009a(3,pwd,pwd);

   Required:
   THREE_SEED_TIMING_INTEGRITY_PASS=1

3. Final timing:
   run_three_methods_end_to_end_timing_R2009a(20,pwd,pwd);

   Required:
   TWENTY_SEED_END_TO_END_TIMING_COMPLETE=1

Outputs
-------
three_methods_end_to_end_timing_*seed_raw.csv
three_methods_end_to_end_timing_*seed_summary.csv
three_methods_end_to_end_timing_*seed_decision.txt
three_methods_end_to_end_timing_*seed_table.tex
three_methods_end_to_end_timing_*seed_checkpoint.mat

The generated LaTeX table is ready to replace the blank timing table in:
CGC_GOCNMF_three_method_validation_and_timing_no_threeparttable.tex

Do not use the old run_seconds column in the manuscript.
