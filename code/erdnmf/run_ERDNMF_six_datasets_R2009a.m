function run_ERDNMF_six_datasets_R2009a(runs, dataDir, outputDir)
%RUN_ERDNMF_SIX_DATASETS_R2009A
% Unified ERDNMF baseline for the CGC-GOCNMF six-dataset protocol.
%
% MATLAB R2009a compatible. No K-means or external classifier is used.
% The class assignment is the largest coordinate in each row of V.
% Primary metrics are computed on UNLABELED samples only.
%
% The multiplicative updates follow the public ERDNMF.m implementation.
% Initialization and labeled-index handling are adapted only to match the
% paired CGC-GOCNMF protocol exactly.
%
% Usage:
%   run_ERDNMF_six_datasets_R2009a(1,  pwd, pwd);  % smoke test
%   run_ERDNMF_six_datasets_R2009a(3,  pwd, pwd);  % integrity gate
%   run_ERDNMF_six_datasets_R2009a(20, pwd, pwd);  % final experiment
%
% The program checkpoints after every dataset-seed pair and can resume.

if nargin < 1 || isempty(runs)
    runs = 1;
end
if nargin < 2 || isempty(dataDir)
    dataDir = pwd;
end
if nargin < 3 || isempty(outputDir)
    outputDir = pwd;
end
if ~(runs == 1 || runs == 3 || runs == 20)
    error('runs must be 1, 3, or 20.');
end
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

EPSILON = 1e-12;
LABEL_FRACTION = 0.10;
MAX_ITERATIONS = 200;
SEED_START = 20260617;

% Alpha values directly listed in the public EDNMF main program:
% PIE=400, COIL20=800, COIL100=1000, Optdigits=600.
% YaleB and MNIST are absent from that benchmark and use fixed same-domain
% transfers from Yale and Optdigits, respectively.
D(1).name = 'PIE';
D(1).files = {'CMU_PIE_fac.mat','CMU_PIE.mat','PIE.mat','PIE_fac.mat'};
D(1).expectedSamples = 2856;
D(1).expectedClasses = 68;
D(1).alpha = 400;
D(1).alphaSource = 'public-code-PIE';

D(2).name = 'YaleB';
D(2).files = {'YaleB.mat','YaleB(1).mat','YaleB_32x32.mat'};
D(2).expectedSamples = 2414;
D(2).expectedClasses = 38;
D(2).alpha = 10;
D(2).alphaSource = 'same-domain-transfer-from-Yale';

D(3).name = 'COIL20';
D(3).files = {'COIL20_Obj.mat','COIL20.mat','COIL20_Obj(1).mat'};
D(3).expectedSamples = 1440;
D(3).expectedClasses = 20;
D(3).alpha = 800;
D(3).alphaSource = 'public-code-COIL20';

D(4).name = 'COIL100';
D(4).files = {'COIL100_Obj.mat','COIL100.mat','COIL100_Obj(1).mat'};
D(4).expectedSamples = 7200;
D(4).expectedClasses = 100;
D(4).alpha = 1000;
D(4).alphaSource = 'public-code-COIL100';

D(5).name = 'Optdigits';
D(5).files = {'Optdigits_Han.mat','Optdigits.mat','optdigits.mat'};
D(5).expectedSamples = 5620;
D(5).expectedClasses = 10;
D(5).alpha = 600;
D(5).alphaSource = 'public-code-Optdigits';

D(6).name = 'MNIST';
D(6).files = {'MNIST_Han.mat','MNIST.mat','mnist.mat'};
D(6).expectedSamples = 6996;
D(6).expectedClasses = 10;
D(6).alpha = 600;
D(6).alphaSource = 'same-domain-transfer-from-Optdigits';

prefix = sprintf('ERDNMF_six_datasets_%dseed',runs);
rawPath = fullfile(outputDir,[prefix '_raw.csv']);
summaryPath = fullfile(outputDir,[prefix '_summary.csv']);
decisionPath = fullfile(outputDir,[prefix '_decision.txt']);
parameterPath = fullfile(outputDir,[prefix '_parameters.csv']);
progressPath = fullfile(outputDir,[prefix '_checkpoint.mat']);

rawTemplate = struct( ...
    'dataset','', 'seed',0, 'method','ERDNMF', 'data_file','', ...
    'samples',0, 'features',0, 'classes',0, 'labeled_count',0, ...
    'alpha',0, 'parameter_source','', 'iterations',MAX_ITERATIONS, ...
    'unlabeled_ACC',0, 'unlabeled_NMI',0, ...
    'all_ACC',0, 'all_NMI',0, 'labeled_ACC',0, ...
    'objective_initial',0, 'objective_final',0, ...
    'objective_nonincrease',0, 'run_seconds',0);

numberOfDatasets = length(D);
totalRows = numberOfDatasets*runs;
rawRows = repmat(rawTemplate,1,totalRows);
completed = false(numberOfDatasets,runs);

if exist(progressPath,'file') == 2
    P = load(progressPath);
    if isfield(P,'rawRows') && isfield(P,'completed') && ...
            length(P.rawRows)==totalRows && ...
            all(size(P.completed)==[numberOfDatasets runs])
        rawRows = P.rawRows;
        completed = P.completed;
        fprintf('Resuming checkpoint: %s\n',progressPath);
    else
        error('Existing checkpoint is incompatible with requested run count.');
    end
end

erd_write_parameter_csv(parameterPath,D,MAX_ITERATIONS);

fprintf('\nERDNMF SIX-DATASET EXPERIMENT -- 2026-08-01\n');
fprintf('Runs per dataset: %d\n',runs);
fprintf('Seeds: %d to %d\n',SEED_START,SEED_START+runs-1);
fprintf('Primary evaluation: unlabeled samples only\n');
fprintf('Direct row-wise argmax; no K-means\n');
fprintf('Iterations: %d\n\n',MAX_ITERATIONS);

for d = 1:numberOfDatasets
    if all(completed(d,:))
        fprintf('%s already complete; skipping.\n',D(d).name);
        continue;
    end

    dataPath = erd_find_file_recursive(dataDir,D(d).files);
    fprintf('\nLoading %s: %s\n',D(d).name,dataPath);
    [X,y] = erd_load_dataset(dataPath,D(d).expectedSamples,D(d).expectedClasses);
    [m,n] = size(X);
    c = length(unique(y));

    for r = 1:runs
        if completed(d,r)
            fprintf('%s seed %d already complete; skipping.\n', ...
                D(d).name,SEED_START+r-1);
            continue;
        end

        seed = SEED_START+r-1;
        fprintf('%s seed %d (%d/%d)\n',D(d).name,seed,r,runs);

        L = erd_labeled_indices(y,seed,LABEL_FRACTION);
        [U0,V0] = erd_initial_factors(m,n,c,seed,EPSILON);

        methodClock = tic;
        [V,obj0,obj1,nonincrease] = erd_train_official( ...
            X,y,L,D(d).alpha,U0,V0,MAX_ITERATIONS,EPSILON);
        runSeconds = toc(methodClock);
        metrics = erd_evaluate(V,y,L,EPSILON);

        idx = (d-1)*runs+r;
        row = rawTemplate;
        row.dataset = D(d).name;
        row.seed = seed;
        row.data_file = erd_filename(dataPath);
        row.samples = n;
        row.features = m;
        row.classes = c;
        row.labeled_count = length(L);
        row.alpha = D(d).alpha;
        row.parameter_source = D(d).alphaSource;
        row.unlabeled_ACC = metrics.unlabeledACC;
        row.unlabeled_NMI = metrics.unlabeledNMI;
        row.all_ACC = metrics.allACC;
        row.all_NMI = metrics.allNMI;
        row.labeled_ACC = metrics.labeledACC;
        row.objective_initial = obj0;
        row.objective_final = obj1;
        row.objective_nonincrease = nonincrease;
        row.run_seconds = runSeconds;
        rawRows(idx) = row;

        completed(d,r) = true;
        save(progressPath,'rawRows','completed','D','runs');
        erd_write_raw_csv(rawPath,rawRows,completed,runs);

        fprintf(['  ACC_U=%.6f NMI_U=%.6f ACC_L=%.6f ' ...
            'obj_ratio=%.6g time=%.2fs\n'], ...
            metrics.unlabeledACC,metrics.unlabeledNMI, ...
            metrics.labeledACC,obj1/max(obj0,EPSILON),runSeconds);
    end
end

summaryRows = erd_summarize(rawRows,D,runs);
erd_write_raw_csv(rawPath,rawRows,completed,runs);
erd_write_summary_csv(summaryPath,summaryRows);
erd_write_decision(decisionPath,runs,rawRows,summaryRows,D);

if all(completed(:))
    if exist(progressPath,'file') == 2
        delete(progressPath);
    end
end

fprintf('\nERDNMF_COMPLETE=%d\n',all(completed(:)));
fprintf('RAW=%s\n',rawPath);
fprintf('SUMMARY=%s\n',summaryPath);
fprintf('DECISION=%s\n',decisionPath);
end


function path = erd_find_file_recursive(folder,candidates)
path = '';
for i=1:length(candidates)
    direct = fullfile(folder,candidates{i});
    if exist(direct,'file')==2
        path = direct;
        return;
    end
end
listing = dir(folder);
for i=1:length(listing)
    if listing(i).isdir && ~strcmp(listing(i).name,'.') && ...
            ~strcmp(listing(i).name,'..')
        child = fullfile(folder,listing(i).name);
        try
            path = erd_find_file_recursive(child,candidates);
        catch
            path = '';
        end
        if ~isempty(path)
            return;
        end
    end
end
error('None of the required files was found under %s.',folder);
end


function name = erd_filename(path)
[dummy,name,ext] = fileparts(path); %#ok<ASGLU>
name = [name ext];
end


function [X,y] = erd_load_dataset(path,expectedSamples,expectedClasses)
S = load(path);
featureNames = {'fea','X','data','features'};
labelNames = {'gnd','labels','label','y','Y','truth'};
fea = [];
gnd = [];
for i=1:length(featureNames)
    if isfield(S,featureNames{i})
        fea = S.(featureNames{i});
        break;
    end
end
for i=1:length(labelNames)
    if isfield(S,labelNames{i})
        gnd = S.(labelNames{i});
        break;
    end
end
if isempty(fea), error('%s: no feature field found.',path); end
if isempty(gnd), error('%s: no label field found.',path); end

F = double(fea);
y0 = double(gnd(:));
if ndims(F)~=2
    error('%s: feature array must be two-dimensional.',path);
end
if size(F,1)==length(y0)
    X = F';
elseif size(F,2)==length(y0)
    X = F;
else
    error('%s: feature shape does not match labels.',path);
end
if min(X(:)) < -1e-12
    error('%s: negative features found; ERDNMF requires X>=0.',path);
end
[values,dummy,y] = unique(y0); %#ok<ASGLU>
y = double(y(:));
if size(X,2)~=expectedSamples || length(values)~=expectedClasses
    error('%s: expected %d samples/%d classes, obtained %d/%d.', ...
        path,expectedSamples,expectedClasses,size(X,2),length(values));
end
norms = sqrt(sum(X.^2,1));
norms(norms==0)=1;
X = bsxfun(@rdivide,X,norms);
if any(~isfinite(X(:)))
    error('%s: non-finite normalized feature found.',path);
end
end


function erd_set_seed(seed)
rand('twister',double(seed)); %#ok<RAND>
end


function L = erd_labeled_indices(y,seed,fraction)
erd_set_seed(seed);
classes = unique(y);
L = zeros(0,1);
for k=1:length(classes)
    ids = find(y==classes(k));
    order = randperm(length(ids));
    count = max(2,floor(fraction*length(ids)));
    count = min(count,length(ids));
    L = [L; ids(order(1:count))]; %#ok<AGROW>
end
L = double(L(:));
end


function [U0,V0] = erd_initial_factors(m,n,c,seed,epsilon)
erd_set_seed(seed+500001);
U0 = max(rand(m,c),epsilon);
V0 = max(rand(n,c),epsilon);
end


function [V,obj0,obj1,nonincrease] = erd_train_official( ...
    X,y,L,alpha,U,V,maxIterations,epsilon)
% Follows the public ERDNMF.m multiplicative update.
[n,c] = size(V);
Y = zeros(n,c);
for t=1:length(L)
    Y(L(t),y(L(t))) = 1;
end
H = zeros(n,c);
H(L,:) = 1-Y(L,:);

obj0 = erd_objective(X,U,V,Y,H,alpha,epsilon);
previous = obj0;
nonincrease = 1;
for iter=1:maxIterations
    numeratorU = X*V;
    denominatorU = U*(V'*V);
    U = U.*(numeratorU./max(denominatorU,epsilon));
    U = max(U,epsilon);

    sumCorrect = sum(sum(Y.*V));
    sumWrong = sum(sum(H.*V));
    sumCorrect = max(sumCorrect,epsilon);

    numeratorV = X'*U + alpha*(sumWrong/(sumCorrect^2))*Y;
    denominatorV = V*(U'*U) + alpha*(1/sumCorrect)*H;
    V = V.*(numeratorV./max(denominatorV,epsilon));
    V = max(V,epsilon);

    current = erd_objective(X,U,V,Y,H,alpha,epsilon);
    if current > previous + 1e-10*max(1,abs(previous))
        nonincrease = 0;
    end
    previous = current;
end
obj1 = previous;
end


function value = erd_objective(X,U,V,Y,H,alpha,epsilon)
residual = X-U*V';
ratio = sum(sum(H.*V))/max(sum(sum(Y.*V)),epsilon);
value = sum(sum(residual.^2)) + 2*alpha*ratio;
end


function M = erd_evaluate(V,y,L,epsilon)
[dummy,prediction] = max(V,[],2); %#ok<ASGLU>
n = length(y);
unlabeledMask = true(n,1);
unlabeledMask(L)=false;
U = find(unlabeledMask);
M.unlabeledACC = mean(prediction(U)==y(U));
M.unlabeledNMI = erd_nmi(y(U),prediction(U),epsilon);
M.allACC = mean(prediction==y);
M.allNMI = erd_nmi(y,prediction,epsilon);
M.labeledACC = mean(prediction(L)==y(L));
end


function value = erd_nmi(trueLabel,predictedLabel,epsilon)
trueLabel = trueLabel(:);
predictedLabel = predictedLabel(:);
[trueValues,dummy1,trueIndex] = unique(trueLabel); %#ok<ASGLU>
[predValues,dummy2,predIndex] = unique(predictedLabel); %#ok<ASGLU>
nt = length(trueValues);
np = length(predValues);
N = length(trueLabel);
contingency = accumarray([trueIndex predIndex],1,[nt np]);
pij = contingency/N;
pi = sum(pij,2);
pj = sum(pij,1);
mi = 0;
for i=1:nt
    for j=1:np
        if pij(i,j)>0
            mi = mi+pij(i,j)*log(pij(i,j)/max(pi(i)*pj(j),epsilon));
        end
    end
end
ht = -sum(pi(pi>0).*log(pi(pi>0)));
hp = -sum(pj(pj>0).*log(pj(pj>0)));
denominator = 0.5*(ht+hp);
if denominator<=epsilon
    value = 1;
else
    value = mi/denominator;
end
value = max(0,min(1,value));
end


function erd_write_parameter_csv(path,D,maxIterations)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'dataset,method,alpha,parameter_source,iterations,label_fraction,seeds,direct_assignment,evaluation\n');
for i=1:length(D)
    fprintf(fid,'%s,ERDNMF,%.15g,%s,%d,0.10,20260617-20260636,row-wise-argmax,unlabeled-only\n', ...
        D(i).name,D(i).alpha,D(i).alphaSource,maxIterations);
end
fclose(fid);
end


function erd_write_raw_csv(path,rows,completed,runs)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,seed,method,data_file,samples,features,classes,' ...
    'labeled_count,alpha,parameter_source,iterations,unlabeled_ACC,' ...
    'unlabeled_NMI,all_ACC,all_NMI,labeled_ACC,objective_initial,' ...
    'objective_final,objective_nonincrease,run_seconds\n']);
for d=1:size(completed,1)
    for r=1:runs
        if completed(d,r)
            idx = (d-1)*runs+r;
            S = rows(idx);
            fprintf(fid,['%s,%d,%s,%s,%d,%d,%d,%d,%.15g,%s,%d,' ...
                '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%.6f\n'], ...
                S.dataset,S.seed,S.method,S.data_file,S.samples,S.features, ...
                S.classes,S.labeled_count,S.alpha,S.parameter_source, ...
                S.iterations,S.unlabeled_ACC,S.unlabeled_NMI,S.all_ACC, ...
                S.all_NMI,S.labeled_ACC,S.objective_initial, ...
                S.objective_final,S.objective_nonincrease,S.run_seconds);
        end
    end
end
fclose(fid);
end


function rows = erd_summarize(rawRows,D,runs)
template = struct('dataset','','method','ERDNMF','runs',runs, ...
    'unlabeled_ACC_mean',0,'unlabeled_ACC_sd',0, ...
    'unlabeled_NMI_mean',0,'unlabeled_NMI_sd',0, ...
    'all_ACC_mean',0,'all_ACC_sd',0, ...
    'all_NMI_mean',0,'all_NMI_sd',0, ...
    'labeled_ACC_mean',0,'objective_nonincrease_pass',0, ...
    'run_seconds_mean',0);
rows = repmat(template,1,length(D));
for d=1:length(D)
    ids = (d-1)*runs+(1:runs);
    A = [rawRows(ids).unlabeled_ACC];
    N = [rawRows(ids).unlabeled_NMI];
    AA = [rawRows(ids).all_ACC];
    NN = [rawRows(ids).all_NMI];
    L = [rawRows(ids).labeled_ACC];
    O = [rawRows(ids).objective_nonincrease];
    T = [rawRows(ids).run_seconds];
    rows(d).dataset = D(d).name;
    rows(d).unlabeled_ACC_mean = mean(A);
    rows(d).unlabeled_ACC_sd = erd_safe_std(A);
    rows(d).unlabeled_NMI_mean = mean(N);
    rows(d).unlabeled_NMI_sd = erd_safe_std(N);
    rows(d).all_ACC_mean = mean(AA);
    rows(d).all_ACC_sd = erd_safe_std(AA);
    rows(d).all_NMI_mean = mean(NN);
    rows(d).all_NMI_sd = erd_safe_std(NN);
    rows(d).labeled_ACC_mean = mean(L);
    rows(d).objective_nonincrease_pass = sum(O==1);
    rows(d).run_seconds_mean = mean(T);
end
end


function value = erd_safe_std(x)
if length(x)>1
    value = std(x,0);
else
    value = 0;
end
end


function erd_write_summary_csv(path,rows)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,method,runs,unlabeled_ACC_mean,unlabeled_ACC_sd,' ...
    'unlabeled_NMI_mean,unlabeled_NMI_sd,all_ACC_mean,all_ACC_sd,' ...
    'all_NMI_mean,all_NMI_sd,labeled_ACC_mean,' ...
    'objective_nonincrease_pass,run_seconds_mean\n']);
for i=1:length(rows)
    S = rows(i);
    fprintf(fid,['%s,%s,%d,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%d,%.6f\n'], ...
        S.dataset,S.method,S.runs,S.unlabeled_ACC_mean, ...
        S.unlabeled_ACC_sd,S.unlabeled_NMI_mean,S.unlabeled_NMI_sd, ...
        S.all_ACC_mean,S.all_ACC_sd,S.all_NMI_mean,S.all_NMI_sd, ...
        S.labeled_ACC_mean,S.objective_nonincrease_pass, ...
        S.run_seconds_mean);
end
fclose(fid);
end


function erd_write_decision(path,runs,rawRows,summaryRows,D)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'ERDNMF SIX-DATASET AUDIT\n');
fprintf(fid,'runs=%d\n',runs);
fprintf(fid,'seeds=20260617-%d\n',20260617+runs-1);
fprintf(fid,'label_fraction=0.10 with at least two labeled samples per class\n');
fprintf(fid,'evaluation=unlabeled samples only\n');
fprintf(fid,'assignment=row-wise argmax; no K-means\n');
fprintf(fid,'iterations=200\n');
fprintf(fid,'implementation=public ERDNMF multiplicative update\n\n');
for d=1:length(D)
    S = summaryRows(d);
    fprintf(fid,['%s alpha=%g source=%s ACC=%.6f+-%.6f ' ...
        'NMI=%.6f+-%.6f labeled_ACC=%.6f objective_nonincrease=%d/%d ' ...
        'time=%.3fs\n'], ...
        D(d).name,D(d).alpha,D(d).alphaSource, ...
        S.unlabeled_ACC_mean,S.unlabeled_ACC_sd, ...
        S.unlabeled_NMI_mean,S.unlabeled_NMI_sd, ...
        S.labeled_ACC_mean,S.objective_nonincrease_pass,runs, ...
        S.run_seconds_mean);
end
finitePass = 1;
for i=1:length(rawRows)
    fields = [rawRows(i).unlabeled_ACC rawRows(i).unlabeled_NMI ...
        rawRows(i).all_ACC rawRows(i).all_NMI rawRows(i).labeled_ACC ...
        rawRows(i).objective_initial rawRows(i).objective_final];
    if any(~isfinite(fields))
        finitePass = 0;
    end
end
fprintf(fid,'\nFINITE_OUTPUT_PASS=%d\n',finitePass);
fprintf(fid,'FINAL_READY=%d\n',runs==20 && finitePass==1);
fclose(fid);
end
