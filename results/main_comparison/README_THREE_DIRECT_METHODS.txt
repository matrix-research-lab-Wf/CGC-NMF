Unified comparison of GNMFLD, GOCNMF, and CGC-GOCNMF
MATLAB R2009a

1. Purpose
----------
This package implements the first-stage unified comparison of three methods
whose representation matrix directly supplies cluster assignments:

    GNMFLD
    GOCNMF
    CGC-GOCNMF

No K-means and no external classifier are used.

2. Fixed protocol
-----------------
Datasets:
    PIE, YaleB, COIL20, COIL100, Optdigits, MNIST

Label fraction:
    10% selected independently within each class

Seeds:
    20260617 onward

Primary evaluation:
    ACC and NMI on unlabeled samples only

Secondary contextual evaluation:
    ACC and NMI on all samples

Iterations:
    GNMFLD: 200
    GOCNMF: 50
    CGC-GOCNMF: 50

Initialization:
    the same U0 and V0 are supplied to all three methods for each paired run

Graph:
    GNMFLD uses a binary 5-nearest-neighbor graph
    GOCNMF uses its locked dataset-specific binary graph
    CGC-GOCNMF calibrates the corresponding local-scale graph against the
    GOCNMF binary graph using the frozen five-fold out-of-fold Brier rule

3. GNMFLD model
---------------
The implementation minimizes

    ||X-U V'||_F^2
    + alpha ||P_Omega(V)-P_Omega(Y)||_F^2
    + beta Tr(V' L V)

using the published multiplicative updates. Prediction is row-wise argmax(V).

4. GNMFLD parameter status
--------------------------
Directly supported by later published baseline reproductions:
    PIE       alpha=1e5, beta=10
    YaleB     alpha=1e5, beta=100
    COIL20    alpha=1e4, beta=10
    MNIST     alpha=1e5, beta=10

Explicit same-domain transfers:
    COIL100   alpha=1e4, beta=10, transferred from COIL20
    Optdigits alpha=1e5, beta=10, transferred from MNIST

The transferred values are not claimed to be original-paper settings.
The generated parameter CSV and decision file preserve this disclosure.

5. Run order
------------
Place all six datasets in the same directory as the program, or in its
immediate subdirectories.

Smoke test:
    clear functions;
    rehash path;
    which run_three_direct_methods_R2009a -all
    run_three_direct_methods_R2009a(1,pwd,pwd);

Required:
    ONE_SEED_SMOKE_PASS=1

Three-seed integrity gate:
    run_three_direct_methods_R2009a(3,pwd,pwd);

Required:
    THREE_SEED_INTEGRITY_PASS=1

Final paired experiment:
    run_three_direct_methods_R2009a(20,pwd,pwd);

Required:
    TWENTY_SEED_COMPARISON_COMPLETE=1

6. Resume behavior
------------------
The program writes a checkpoint MAT file after every dataset-seed pair.
Rerunning the same command resumes incomplete work and does not repeat completed
pairs. Delete the corresponding checkpoint only when intentionally restarting
the whole experiment.

7. Outputs
----------
three_direct_methods_*seed_raw.csv
three_direct_methods_*seed_summary.csv
three_direct_methods_*seed_paired.csv
three_direct_methods_*seed_decision.txt
three_direct_methods_*seed_parameters.csv
three_direct_methods_*seed_checkpoint.mat

8. Interpretation discipline
----------------------------
The 1-seed and 3-seed runs are implementation gates, not performance evidence.
Only the 20-seed paired output may be used for the final comparison table.

A weak GNMFLD result must not be hidden. First inspect:
    labeled_ACC
    objective_nonincrease
    parameter_source
    data dimensions and class count

Do not tune GNMFLD using unlabeled ground-truth labels.
