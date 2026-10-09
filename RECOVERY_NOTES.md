# Recovery and validation notes

## What was missing

No physical copies of `run_common_stop_core_R2019a.m`,
`run_SNMFWLP_common_stop_R2019a.m`, or their latest three-seed outputs remained in
the searched local source directories.

## How the scripts were recovered

The original Codex session log retained the base-file copies and every subsequent
patch applied during the September 2026 experiment work. The scripts in
`code/recovered_latest` were reconstructed by replaying that recorded patch history.
They are therefore recovered historical implementations, not newly invented
substitutes. Helper drivers added in October 2026 are clearly named as recovery
drivers.

## MATLAB R2019a validation

- MATLAB executable: `D:\Matlab R2019a\bin\matlab.exe`
- Formal seeds: 20260617, 20260618, and 20260619
- Evaluation: unlabeled samples only
- Assignment: direct row-wise argmax, without K-means
- Common stopping rule: residual tolerance and patience recorded by each raw CSV
- Datasets: PIE, YaleB, COIL20, Optdigits, MNIST, and COIL100

All six common-stopping decision files report finite results, direct-argmax use,
objective nonincrease, fixed labeled rows, three-method integrity, and
`THREE_SEED_INTEGRITY_PASS=1`.

All eighteen SNMFWLP runs converged under the recorded residual rule and returned
finite outputs. The objective-nonincreasing audit passed for all three PIE runs but
failed for the other fifteen runs. This is a disclosed algorithm/audit result, not
a missing-file error, and no monotonicity claim should be made for those runs.

## Numerical and resource observations

MATLAB emitted ill-conditioning or singular-matrix warnings in some harmonic label
propagation solves, particularly for PIE and COIL100. These warnings are retained in
the logs. COIL100 also exceeded available memory when the full reconstruction
residual matrix was formed repeatedly. The recovered code now accumulates the same
Frobenius reconstruction error in blocks of 512 samples. This changes summation
order only and avoids the full temporary residual matrix. COIL100 SNMFWLP was run in
a separate fresh `-nojvm` MATLAB process to avoid memory fragmentation.

## Output locations

- `results/recovered_six_dataset_3seed/common_stop_<dataset>`
- `results/recovered_six_dataset_3seed/snmfwlp_<dataset>`
- `results/recovered_common_stop_3seed_summary.csv`
- `results/recovered_snmfwlp_3seed_summary.csv`
- `results/recovered_3seed_audit.txt`

These regenerated results must still be reconciled line by line with the final
manuscript tables before submission.
