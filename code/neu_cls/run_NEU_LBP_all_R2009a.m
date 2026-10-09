function run_NEU_LBP_all_R2009a(imageRoot,outputDir)
%RUN_NEU_LBP_ALL_R2009A
% Prepare NEU-CLS LBP features, run the three-seed integrity test, and
% then run the formal 20-seed GOCNMF/CGC-GOCNMF comparison.
%
% Usage:
%   run_NEU_LBP_all_R2009a;
%
% MATLAB R2009a compatible.

if nargin<1 || isempty(imageRoot)
    imageRoot=uigetdir(pwd,'Select the NEU-CLS image folder');
    if isequal(imageRoot,0)
        error('Image-folder selection was cancelled.');
    end
end

if nargin<2 || isempty(outputDir)
    outputDir=uigetdir(pwd,'Select the output folder');
    if isequal(outputDir,0)
        error('Output-folder selection was cancelled.');
    end
end

if exist(outputDir,'dir')~=7
    mkdir(outputDir);
end

programDir=fileparts(mfilename('fullpath'));
addpath(programDir);
rehash path;

featureFile=fullfile(outputDir,'NEU_CLS_LBP59_4x4.mat');
if exist(featureFile,'file')~=2
    prepare_NEU_CLS_LBP_R2009a(imageRoot,featureFile);
else
    fprintf('Using existing LBP feature file: %s\n',featureFile);
end

fprintf('\nSTEP 1/2: three-seed integrity test.\n');
run_NEU_LBP_GOC_CGC_R2009a(3,outputDir,outputDir);

fprintf('\nSTEP 2/2: formal 20-seed experiment.\n');
run_NEU_LBP_GOC_CGC_R2009a(20,outputDir,outputDir);

fprintf('\nNEU_LBP_ALL_STAGES_COMPLETE=1\n');
end
