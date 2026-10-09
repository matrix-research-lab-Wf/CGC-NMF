# Data setup

Only `main_six/Optdigits_Han.mat` is distributed in this repository because its
authoritative UCI record provides a CC BY 4.0 license. See
`../DATA_PROVENANCE.md` for attribution, modification details, expected shapes,
checksums, and the reasons other matrices are excluded.

To reproduce all archived experiments, obtain the remaining datasets from their
authoritative sources and place locally prepared matrices at these paths:

```text
data/main_six/PIE.mat
data/main_six/YaleB.mat
data/main_six/COIL20_Obj.mat
data/main_six/COIL100_Obj.mat
data/main_six/MNIST_Han.mat
data/neu_cls/NEU_CLS_32x32.mat
data/neu_cls/NEU_CLS_LBP59_4x4.mat
```

Every main matrix must contain `fea` (samples by nonnegative features) and `gnd`
(samples by one class label). Exact author-side dimensions and SHA-256 values are
listed in `../DATA_PROVENANCE.md`. The NEU-CLS preparation routines under
`code/neu_cls` generate the application matrices after the official images have
been downloaded locally.

Do not commit locally acquired restricted datasets. The repository `.gitignore`
explicitly blocks their expected filenames.
