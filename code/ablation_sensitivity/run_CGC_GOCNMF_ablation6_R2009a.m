function run_CGC_GOCNMF_ablation6_R2009a(runs,dataDir,outputDir)
%RUN_CGC_GOCNMF_ABLATION6_R2009A
% Six-dataset BASE/LOCAL/HARD/SHRINK ablation.
% MATLAB R2009a compatible.
%
% Usage:
%   run_CGC_GOCNMF_ablation6_R2009a(3,pwd,pwd);
%   run_CGC_GOCNMF_ablation6_R2009a(20,pwd,pwd);

if nargin < 1 || isempty(runs), runs = 3; end
if nargin < 2 || isempty(dataDir), dataDir = pwd; end
if nargin < 3 || isempty(outputDir), outputDir = pwd; end

CGC_GOCNMF_revision_core_R2009a('ablation6',runs,dataDir,outputDir);
end
