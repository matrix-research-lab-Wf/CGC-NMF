# Pre-upload checklist

Target repository: <https://github.com/matrix-research-lab-Wf/CGC-NMF>

## Hard blockers

- [x] Resolve the official source, license, and redistribution permission for all
  eight MAT files in `DATA_PROVENANCE.md`.
- [x] Remove any dataset that cannot legally be redistributed and replace it with
  an official download link, deterministic preparation script, preprocessing
  description, and expected checksum.
- [x] Correct the false SNMFWLP stopping sentence identified in
  `FINAL_MANUSCRIPT_RESULT_AUDIT.md`, then rebuild and visually inspect the PDF.
- [x] Regenerate and archive the exact Wilcoxon--Holm output for the CGC-GNMFLD
  transfer experiment so that code, output, and manuscript report the same test.
- [x] Archive the code, seed, and generated output for the reported 10,000-resample
  paired bootstrap confidence intervals.
- [x] Replace the license-unclear Hungarian implementation and document external
  method provenance in `THIRD_PARTY_NOTICES.md`.

## Final verification

- [x] Run `tools/verify_before_upload.ps1`; it returned exit code 0.
- [x] Run `run_public_smoke_test` using MATLAB R2019a; the final clean-clone
  repetition is recorded in the release commit workflow.
- [x] Confirm all seven manuscript figures are present and inspect their rendered
  placement. Full pixel-level regeneration remains outside this release because
  several source datasets cannot be redistributed.
- [x] Recompile the exact manuscript source twice and inspect every PDF page.
- [x] Confirm that no credentials, access tokens, private keys, or personal data
  are present.
- [x] Regenerate `MANIFEST_SHA256.csv` after the final file set is frozen.
- [x] Add the public repository URL to the manuscript's data and code availability
  statement. The immutable commit identifier will be recorded in the Git history.
- [x] Add official download/conversion programs and the recovered MNIST_Han row
  indices; verify exact MATLAB R2019a equality for MNIST_Han and COIL20.

## Current decision

The license-checked public package may be pushed after the automated verifier and
clean-clone public smoke test pass. Restricted matrices remain outside Git history.
