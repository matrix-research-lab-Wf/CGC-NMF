ERDNMF baseline for revising Table 6 (NEU-CLS)

Purpose
-------
Add ERDNMF to the existing 5%, 10%, and 20% NEU-CLS comparison without
using unlabeled ground-truth labels for parameter selection.

Required data
-------------
NEU_CLS_32x32.mat containing:
- fea: 1800 x 1024 or 1024 x 1800 nonnegative features
- gnd: 1800 labels in six classes

Stage 1: labeled-only alpha selection
-------------------------------------
Run one smoke test:

clear functions;
rehash path;
run_ERDNMF_NEU_labeledCV_select_alpha_R2009a(1,pwd,pwd);

Then the full pilot selection:

run_ERDNMF_NEU_labeledCV_select_alpha_R2009a(20,pwd,pwd);

Upload/check:
- ERDNMF_NEU_labeledCV_alpha_selection_20seed_raw.csv
- ERDNMF_NEU_labeledCV_alpha_selection_20seed_summary.csv
- ERDNMF_NEU_labeledCV_alpha_selection_20seed_selected.csv
- ERDNMF_NEU_labeledCV_alpha_selection_20seed_decision.txt

The formal experiment should be run only when FINAL_READY=1.

Stage 2: formal Table-6 experiment
----------------------------------
Run:

run_ERDNMF_NEU_formal_R2009a(1,pwd,pwd, ...
    'ERDNMF_NEU_labeledCV_alpha_selection_20seed_selected.csv');

After the smoke test passes:

run_ERDNMF_NEU_formal_R2009a(20,pwd,pwd, ...
    'ERDNMF_NEU_labeledCV_alpha_selection_20seed_selected.csv');

Upload:
- ERDNMF_NEU_formal_20seed_raw.csv
- ERDNMF_NEU_formal_20seed_summary.csv
- ERDNMF_NEU_formal_20seed_parameters.csv
- ERDNMF_NEU_formal_20seed_decision.txt

The formal program also writes ERDNMF confusion matrices for Figure 1.
After these files are returned, Table 6, Holm-adjusted tests, the surrounding
text, and the four-method Figure 1 can be finalized.
