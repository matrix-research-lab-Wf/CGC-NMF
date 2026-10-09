function run_smoke_test()
%RUN_SMOKE_TEST Run the one-seed six-dataset reproducibility checks.
% This entry point is compatible with MATLAB R2019a and writes new outputs
% below results/reproduced_smoke without modifying the archived results.

rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir,'data','main_six');
outputDir = fullfile(rootDir,'results','reproduced_smoke');

addpath(fullfile(rootDir,'code','core'));
addpath(fullfile(rootDir,'code','main_comparison'));
addpath(fullfile(rootDir,'code','erdnmf'));
addpath(fullfile(rootDir,'code','ablation_sensitivity'));

if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

fprintf('Running the one-seed three-method comparison...\n');
run_three_direct_methods_R2009a(1,dataDir,outputDir);

fprintf('Running the one-seed ERDNMF comparison...\n');
run_ERDNMF_six_datasets_R2009a(1,dataDir,outputDir);

fprintf('Running the one-seed ablation check...\n');
run_CGC_GOCNMF_ablation6_R2009a(1,dataDir,outputDir);

fprintf('Smoke tests completed. Outputs: %s\n',outputDir);
end
