function run_public_smoke_test()
%RUN_PUBLIC_SMOKE_TEST One-seed test using the redistributable Optdigits data.
% Compatible with MATLAB R2019a and base MATLAB.

rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir,'data','main_six');
outputDir = fullfile(rootDir,'results','reproduced_public_smoke');

addpath(fullfile(rootDir,'code','core'));
addpath(fullfile(rootDir,'code','recovered_latest'));
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

run_common_stop_core_R2019a(1,dataDir,outputDir,'Optdigits');
fprintf('Public smoke test completed. Outputs: %s\n',outputDir);
end
