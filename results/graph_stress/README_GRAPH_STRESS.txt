Controlled graph-reliability stress test for CGC-GOCNMF
MATLAB R2009a

Purpose
-------
Test whether the calibration coefficient responds to controlled deterioration
of the candidate graph, rather than merely fitting one favorable dataset.

Dataset
-------
COIL100 only. This dataset is used because the clean local graph has already
shown a stable positive contribution.

Accepted data files
-------------------
COIL100_Obj(1).mat
COIL100_Obj.mat
COIL100.mat

Controlled corruption
---------------------
Let W1 be the clean locally scaled graph. Its undirected nonzero weights are
randomly permuted over the same edge support to form Wperm. For each run,

    Wrho = (1-rho) W1 + rho Wperm,

where rho = 0, 0.2, 0.4, 0.6, 0.8, 1.

The perturbation preserves:
- the p-nearest-neighbor edge support;
- symmetry and nonnegativity;
- the multiset and mean of edge weights at rho=1;
- the original NMF model, labels, alpha, p, and initialization protocol.

Only the correspondence between edge location and local similarity is
progressively destroyed.

Step 1: three-seed integrity run
-------------------------------
clear functions;
rehash path;
which run_GCGOCNMF_graph_stress_R2009a -all
run_GCGOCNMF_graph_stress_R2009a(3, pwd, pwd);

Required before the 20-seed run:
STRESS_FINITE_PASS=1
CORRUPTION_CONSTRUCTION_PASS=1
THREE_SEED_INTEGRITY_PASS=1

TREND_RESPONSE_PASS and RISK_CONTROL_PASS are displayed at three seeds but
are not used to stop the experiment, because three trials are too few for a
stable mechanism decision.

Step 2: twenty-seed stress test
-------------------------------
run_GCGOCNMF_graph_stress_R2009a(20, pwd, pwd);

The predeclared mechanism checks are:
1. theta at rho=1 is lower than theta at rho=0;
2. Delta at rho=1 is lower than Delta at rho=0;
3. candidate ACC and NMI gains at rho=1 are below their rho=0 values;
4. at rho=1, SHRINK has no more losses than the corrupted candidate;
5. at rho=1, SHRINK is no farther from BASE than the corrupted candidate.

The final gate is:
TWENTY_SEED_STRESS_PASS=1

Output
------
GCGOCNMF_graph_stress_*seed_raw.csv
GCGOCNMF_graph_stress_*seed_summary.csv
GCGOCNMF_graph_stress_*seed_decision.txt
GCGOCNMF_graph_stress_*seed_results.mat
GCGOCNMF_graph_stress_*seed_theta.png
GCGOCNMF_graph_stress_*seed_cv.png
GCGOCNMF_graph_stress_*seed_dACC.png
GCGOCNMF_graph_stress_*seed_dNMI.png

Interpretation discipline
-------------------------
This experiment tests whether the calibration reacts to a controlled loss of
edge-weight reliability. It does not establish a universal performance
guarantee. A failed trend gate must be reported, not repaired by changing rho,
the corruption rule, or the pass criteria after seeing the results.
