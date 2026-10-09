FIX FOR MATLAB R2009a DLMWRITE ERROR

Problem:
  Invalid attribute tag: ,

Cause:
  MATLAB R2009a does not accept the mixed dlmwrite syntax
  dlmwrite(path,M,',','precision','%.15g').

Fix:
  The two dlmwrite calls were replaced by a version-independent CSV writer
  based on fopen/fprintf.

No formal experiment needs to be recomputed.
Keep the existing checkpoint file in the output directory, replace the old
run_ERDNMF_NEU_formal_R2009a.m with this corrected file, then rerun the same
20-run command. The program will resume from the completed checkpoint and
only regenerate the output files.
