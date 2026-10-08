PCR-NMF manuscript with Fig. 5--8
==================================

Main source:
  PCR_NMF_Fig5_8_revised.tex

Compiled manuscript:
  PCR_NMF_Fig5_8_revised.pdf

New figures:
  Fig. 5  PCR_Fig5_label_ratio.png
  Fig. 6  PCR_Fig6_COIL20_snmfwlp_style.png
  Fig. 7  PCR_Fig7_USPS_snmfwlp_style.png
  Fig. 8  PCR_Fig8_basis_vectors.png

Compilation check:
  Successfully compiled with MiKTeX pdfLaTeX on 2026-10-07.
  The current output contains 41 pages. Fig. 5 is on page 34; Fig. 6 and Fig. 7
  are on page 35; Fig. 8 is on page 36. Cross-references were resolved after
  three passes.

Notes:
  1. The original manuscript files were not overwritten.
  2. The existing duplicate convergence label was corrected.
  3. Existing minor overfull-box and elsarticle/caption warnings remain;
     no new compilation error or unresolved reference was introduced.
  4. Ground-truth labels in Fig. 6--7 are used only for point colors.
     The manuscript does not claim universal visual superiority.
  5. Repeated defensive wording was replaced by verified ranks, comparison
     counts, effect differences, and significance statements. Mathematical
     qualifications and visualization-protocol boundaries were retained.
  6. The separate same-graph paired-statistics table was removed while its
     main Holm-corrected findings were retained in the text. Related-work and
     theory subsection levels were consolidated, and the statistical table
     was allowed to float to improve pagination.
  7. Fig. 8 was generated in MATLAB R2019a using fixed seed 20260617, the
     formal 10%/5% COIL20/UMIST labeling rates, and the locked method
     parameters. Per-basis min-max normalization is used only for display.
  8. Table 1 now distinguishes explicit label propagation from indirect graph
     coupling, identifies ERDNMF as the element-ratio EDNMF formulation, and
     specifies the supervised object for each method.
