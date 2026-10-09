# PAA public-release checklist

The following items reflect the Pattern Analysis and Applications reproducibility
requirements checked on 8 October 2026.

- [x] Code and locally available data are organized in one directory.
- [x] An English README gives a reproducibility entry point.
- [x] MATLAB source comments contain no detected non-ASCII text.
- [x] Main data files are deduplicated.
- [x] Archived raw and summary result files are included.
- [x] SHA-256 checksums are supplied in `MANIFEST_SHA256.csv`.
- [x] Run recovered-script execution checks in MATLAB R2019a from this directory.
- [x] Restore common-stopping and SNMFWLP files from the original patch history.
- [x] Regenerate six-dataset, three-seed common-stopping and SNMFWLP outputs.
- [ ] Reconcile every manuscript table and figure with a named output file.
- [ ] Verify dataset provenance and redistribution rights.
- [ ] Document third-party code licenses and acknowledgments.
- [ ] Choose a license for author-owned code.
- [ ] Push the uncompressed directory contents to a public repository.
- [ ] Add the public repository URL to the manuscript Declaration section.
- [ ] Verify the repository from a clean clone on another computer.

The ZIP archive prepared locally is only a transfer copy. PAA states that the
repository files must not be provided only as an archive.
