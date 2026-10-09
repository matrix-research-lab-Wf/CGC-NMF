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

## Official download and conversion programs

MATLAB R2019a programs are provided under `code/data_preparation`:

```matlab
addpath(genpath(fullfile(pwd,'code')));
prepare_all_official_data_R2019a( ...
    fullfile(pwd,'official_raw'), fullfile(pwd,'data'));
```

The combined program automatically downloads COIL-20, COIL-100, and MNIST from
the provider or an authorized MNIST mirror. It also attempts the official
NEU-CLS Google Drive download. CMU PIE and Cropped Extended Yale B require the
provider's access/permission procedure; after lawful extraction to
`official_raw/PIE` and `official_raw/YaleB`, the same program converts them.

`MNIST_Han.mat` and `COIL20_Obj.mat` produced by these programs were compared
element by element with the author-side matrices in MATLAB R2019a and matched
exactly. `mnist_han_indices.csv` is the recovered 6996-row public index list; it
contains no image pixels. The legacy COIL-100, PIE, and YaleB MAT files were
created by older pipelines whose complete grayscale/cropping metadata were not
preserved. The supplied programs give deterministic 32-by-32 matrices with the
required sample/class structure, but byte identity with those three archived
MAT files is not claimed. Use `verify_prepared_datasets_R2019a` to report both
structural validation and, when reference matrices are locally available,
element-wise equality.

The NEU-CLS official-image wrapper completed on all 1800 images in MATLAB
R2019a. Its 32-by-32 features agreed with the archived matrix to a maximum
absolute difference of `7.22e-15` (floating-point roundoff). The derived LBP
matrix passed dimension, nonnegativity, and class-balance checks, but did not
match the archived LBP matrix exactly; this difference remains visible rather
than being treated as an exact reproduction.

Do not commit locally acquired restricted datasets. The repository `.gitignore`
explicitly blocks their expected filenames.
