function run_SNMFWLP_common_stop_R2019a(runs, dataDir, outputDir, datasetName,alphaOverride,betaOverride,validationMode,seedStart)
%RUN_SNMFWLP_COMMON_STOP_R2019A
% Integration of the authors' public SNMFWLP update equations into the
% CGC-GOCNMF paired protocol. The public repository supplies only the core
% SNMFWLP.m function; graph construction, splitting, stopping, evaluation,
% and reporting are therefore explicitly provided by this wrapper.
%
% MATLAB R2009a compatible. No K-means or external classifier is used.
% The class assignment is the largest coordinate in each row of V.
% Primary metrics are computed on UNLABELED samples only.
%
% The multiplicative updates follow the public SNMFWLP.m implementation.
% Initialization, graph construction, labeled-index handling, and stopping
% are supplied by this transparent integration wrapper.
%
% Usage:
%   run_SNMFWLP_common_stop_R2019a(1,pwd,pwd);  % smoke test
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
if nargin < 4
    datasetName = '';
end
if nargin < 5, alphaOverride = []; end
if nargin < 6, betaOverride = []; end
if nargin < 7 || isempty(validationMode), validationMode = false; end
if nargin < 8 || isempty(seedStart), seedStart = 20260617; end
if ~(runs == 1 || runs == 3 || runs == 20)
    error('runs must be 1, 3, or 20.');
end
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

EPSILON = 1e-12;
LABEL_FRACTION = 0.10;
MAX_ITERATIONS = 1000;
MIN_ITERATIONS = 20;
STOP_TOLERANCE = 1e-4;
STOP_PATIENCE = 5;
GRAPH_P = 5;
BLOCK_SIZE = 256;
SEED_START = seedStart;

% Alpha values directly listed in the public EDNMF main program:
% PIE=400, COIL20=800, COIL100=1000, Optdigits=600.
% YaleB and MNIST are absent from that benchmark and use fixed same-domain
% transfers from Yale and Optdigits, respectively.
D(1).name = 'PIE';
D(1).files = {'CMU_PIE_fac.mat','CMU_PIE.mat','PIE.mat','PIE_fac.mat'};
D(1).expectedSamples = 2856;
D(1).expectedClasses = 68;
D(1).alpha = 1; D(1).beta = 1;
D(1).alphaSource = 'smoke-default-before-labeled-only-pilot';

D(2).name = 'YaleB';
D(2).files = {'YaleB.mat','YaleB(1).mat','YaleB_32x32.mat'};
D(2).expectedSamples = 2414;
D(2).expectedClasses = 38;
D(2).alpha = 1; D(2).beta = 1;
D(2).alphaSource = 'smoke-default-before-labeled-only-pilot';

D(3).name = 'COIL20';
D(3).files = {'COIL20_Obj.mat','COIL20.mat','COIL20_Obj(1).mat'};
D(3).expectedSamples = 1440;
D(3).expectedClasses = 20;
D(3).alpha = 1; D(3).beta = 1;
D(3).alphaSource = 'smoke-default-before-labeled-only-pilot';

D(4).name = 'COIL100';
D(4).files = {'COIL100_Obj.mat','COIL100.mat','COIL100_Obj(1).mat'};
D(4).expectedSamples = 7200;
D(4).expectedClasses = 100;
D(4).alpha = 1; D(4).beta = 1;
D(4).alphaSource = 'smoke-default-before-labeled-only-pilot';

D(5).name = 'Optdigits';
D(5).files = {'Optdigits_Han.mat','Optdigits.mat','optdigits.mat'};
D(5).expectedSamples = 5620;
D(5).expectedClasses = 10;
D(5).alpha = 1; D(5).beta = 1;
D(5).alphaSource = 'smoke-default-before-labeled-only-pilot';

D(6).name = 'MNIST';
D(6).files = {'MNIST_Han.mat','MNIST.mat','mnist.mat'};
D(6).expectedSamples = 6996;
D(6).expectedClasses = 10;
D(6).alpha = 1; D(6).beta = 1;
D(6).alphaSource = 'smoke-default-before-labeled-only-pilot';

if ~isempty(datasetName)
    keep = false(1,length(D));
    for q=1:length(D), keep(q)=strcmpi(D(q).name,datasetName); end
    if ~any(keep), error('Unknown datasetName: %s',datasetName); end
    D = D(keep);
end
if ~isempty(alphaOverride)
    for q=1:length(D), D(q).alpha=alphaOverride; end
end
if ~isempty(betaOverride)
    for q=1:length(D), D(q).beta=betaOverride; end
end
if ~isempty(alphaOverride) || ~isempty(betaOverride)
    for q=1:length(D)
        D(q).alphaSource = 'disjoint-labeled-only-pilot';
    end
end

prefix = sprintf('SNMFWLP_common_stop_%dseed',runs);
rawPath = fullfile(outputDir,[prefix '_raw.csv']);
summaryPath = fullfile(outputDir,[prefix '_summary.csv']);
decisionPath = fullfile(outputDir,[prefix '_decision.txt']);
parameterPath = fullfile(outputDir,[prefix '_parameters.csv']);
progressPath = fullfile(outputDir,[prefix '_checkpoint.mat']);

rawTemplate = struct( ...
    'dataset','', 'seed',0, 'method','SNMFWLP', 'data_file','', ...
    'samples',0, 'features',0, 'classes',0, 'labeled_count',0, ...
    'alpha',0, 'beta',0, 'parameter_source','', 'iterations',0, ...
    'unlabeled_ACC',0, 'unlabeled_NMI',0, ...
    'all_ACC',0, 'all_NMI',0, 'labeled_ACC',0, 'validation_ACC',NaN, ...
    'objective_initial',0, 'objective_final',0, ...
    'objective_nonincrease',0, 'stop_residual',0, 'converged',0, ...
    'graph_seconds',0, 'run_seconds',0);

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

fprintf('\nSNMFWLP COMMON-STOP EXPERIMENT\n');
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

    graphClock = tic;
    [neighborIndex,neighborDistance] = erd_knn(X,GRAPH_P,BLOCK_SIZE);
    LS = erd_local_graph(neighborIndex,neighborDistance,n,EPSILON);
    graphSeconds = toc(graphClock);

    for r = 1:runs
        if completed(d,r)
            fprintf('%s seed %d already complete; skipping.\n', ...
                D(d).name,SEED_START+r-1);
            continue;
        end

        seed = SEED_START+r-1;
        fprintf('%s seed %d (%d/%d)\n',D(d).name,seed,r,runs);

        L = erd_labeled_indices(y,seed,LABEL_FRACTION);
        validationIndices = zeros(0,1);
        if validationMode
            [L,validationIndices] = erd_split_labeled_validation(y,L,seed+700001);
        end
        [U0,V0] = erd_initial_factors(m,n,c,seed,EPSILON);

        methodClock = tic;
        [V,obj0,obj1,nonincrease,itUsed,stopResidual,converged] = ...
            snm_train_public_updates(X,y,L,LS,D(d).alpha,D(d).beta, ...
            U0,V0,MAX_ITERATIONS,MIN_ITERATIONS,STOP_TOLERANCE, ...
            STOP_PATIENCE,EPSILON);
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
        row.beta = D(d).beta;
        row.parameter_source = D(d).alphaSource;
        row.iterations = itUsed;
        row.unlabeled_ACC = metrics.unlabeledACC;
        row.unlabeled_NMI = metrics.unlabeledNMI;
        row.all_ACC = metrics.allACC;
        row.all_NMI = metrics.allNMI;
        row.labeled_ACC = metrics.labeledACC;
        if validationMode
            [dummy,prediction] = max(V,[],2); %#ok<ASGLU>
            row.validation_ACC = mean(prediction(validationIndices)==y(validationIndices));
        end
        row.objective_initial = obj0;
        row.objective_final = obj1;
        row.objective_nonincrease = nonincrease;
        row.stop_residual = stopResidual;
        row.converged = converged;
        row.graph_seconds = graphSeconds;
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
erd_write_decision(decisionPath,runs,rawRows,summaryRows,D,MAX_ITERATIONS);

if all(completed(:))
    if exist(progressPath,'file') == 2
        delete(progressPath);
    end
end

fprintf('\nSNMFWLP_COMPLETE=%d\n',all(completed(:)));
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
    error('%s: negative features found; SNMFWLP requires X>=0.',path);
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


function [trainIndices,validationIndices] = erd_split_labeled_validation(y,L,seed)
% Stratified 50/50 split of the already selected labeled subset.  Only the
% training half enters SNMFWLP; the held-out half is used for pilot scoring.
erd_set_seed(seed);
classes = unique(y(L));
trainIndices = zeros(0,1);
validationIndices = zeros(0,1);
for k=1:length(classes)
    ids = L(y(L)==classes(k));
    order = randperm(length(ids));
    validationCount = max(1,floor(length(ids)/2));
    validationCount = min(validationCount,length(ids)-1);
    validationIndices = [validationIndices; ids(order(1:validationCount))]; %#ok<AGROW>
    trainIndices = [trainIndices; ids(order(validationCount+1:end))]; %#ok<AGROW>
end
trainIndices = double(trainIndices(:));
validationIndices = double(validationIndices(:));
end


function [U0,V0] = erd_initial_factors(m,n,c,seed,epsilon)
erd_set_seed(seed+500001);
U0 = max(rand(m,c),epsilon);
V0 = max(rand(n,c),epsilon);
end


function [neighborIndex,neighborDistance] = erd_knn(X,p,blockSize)
n = size(X,2);
neighborIndex = zeros(n,p);
neighborDistance = zeros(n,p);
for first=1:blockSize:n
    last=min(first+blockSize-1,n);
    ids=first:last;
    similarities=X(:,ids)'*X;
    for q=1:length(ids), similarities(q,ids(q))=-Inf; end
    [sortedSimilarity,sortedIndex]=sort(similarities,2,'descend');
    neighborIndex(ids,:)=sortedIndex(:,1:p);
    neighborDistance(ids,:)=max(0,1-sortedSimilarity(:,1:p));
end
end


function W = erd_local_graph(neighborIndex,neighborDistance,n,epsilon)
p=size(neighborIndex,2);
rows=repmat((1:n)',1,p);
cols=neighborIndex;
sigma=max(neighborDistance(:,p),epsilon);
denominator=sqrt(sigma(rows(:)).*sigma(cols(:)))+epsilon;
weights=exp(-(neighborDistance(:)./denominator).^2);
A=sparse(rows(:),cols(:),weights,n,n);
W=max(A,A');
W=W-spdiags(diag(W),0,n,n);
values=nonzeros(W);
if ~isempty(values), W=W/mean(values); end
W=sparse(W);
end


function [V,obj0,obj1,nonincrease,iterationsUsed,stopResidual,converged] = ...
    snm_train_public_updates(X,y,L,LS,alpha,beta,U,V,maxIterations, ...
    minIterations,stopTolerance,stopPatience,epsilon)
% The two multiplicative updates below are transcribed from the authors'
% public SNMFWLP.m. Only stopping and numerical guards are supplied here.
[n,c] = size(V);
Y = zeros(n,c);
for t=1:length(L)
    Y(L(t),y(L(t))) = 1;
end
mask = zeros(n,1); mask(L)=1;
degree = full(sum(LS,2));

obj0 = snm_objective(X,U,V,LS,degree,Y,mask,alpha,beta);
previous = obj0;
nonincrease = 1;
stableCount = 0; converged = 0; stopResidual = Inf;
for iter=1:maxIterations
    oldPrediction = snm_row_normalize(V,epsilon);

    XV = X*V; XV(XV<0)=0;
    % Algebraically identical to (U*U')*XV, but avoids the dense m-by-m
    % intermediate used by the public expression and is essential on the
    % 8-GB MATLAB R2019a test machine.
    denominatorU = U*(U'*XV);
    U = U.*sqrt(XV./max(denominatorU,epsilon));
    U = max(U,epsilon);

    numeratorV = X'*U + alpha*(LS*V) + ...
        beta*bsxfun(@times,mask,Y);
    numeratorV(numeratorV<0)=0;
    denominatorV = V*(U'*U) + ...
        alpha*bsxfun(@times,degree,V) + ...
        beta*bsxfun(@times,mask,V);
    V = V.*(numeratorV./max(denominatorV,epsilon));
    V = max(V,epsilon);

    newPrediction = snm_row_normalize(V,epsilon);
    stopResidual = norm(newPrediction-oldPrediction,'fro')/ ...
        max(norm(oldPrediction,'fro'),epsilon);
    if iter==1 || mod(iter,10)==0 || iter==maxIterations
        current = snm_objective(X,U,V,LS,degree,Y,mask,alpha,beta);
        if current > previous + 1e-8*max(1,abs(previous))
            nonincrease = 0;
        end
        previous = current;
    end
    if iter>=minIterations && stopResidual<=stopTolerance
        stableCount=stableCount+1;
    else
        stableCount=0;
    end
    if stableCount>=stopPatience
        converged=1;
        break;
    end
end
iterationsUsed=iter;
obj1 = snm_objective(X,U,V,LS,degree,Y,mask,alpha,beta);
end


function value = snm_objective(X,U,V,LS,degree,Y,mask,alpha,beta)
XV = X*V;
reconstruction = norm(X,'fro')^2-2*sum(sum(U.*XV))+ ...
    sum(sum((U'*U).*(V'*V)));
reconstruction = max(reconstruction,0);
graphTerm = sum(sum(bsxfun(@times,degree,V).*V))-sum(sum((LS*V).*V));
labelDifference = bsxfun(@times,mask,V-Y);
value = reconstruction+alpha*graphTerm+ ...
    beta*sum(labelDifference(:).^2);
end


function P = snm_row_normalize(V,epsilon)
P = bsxfun(@rdivide,V,max(sum(V,2),epsilon));
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
fprintf(fid,'dataset,method,alpha,beta,parameter_source,max_iterations,label_fraction,seeds,direct_assignment,evaluation\n');
for i=1:length(D)
    fprintf(fid,'%s,SNMFWLP,%.15g,%.15g,%s,%d,0.10,20260617-20260636,row-wise-argmax,unlabeled-only\n', ...
        D(i).name,D(i).alpha,D(i).beta,D(i).alphaSource,maxIterations);
end
fclose(fid);
end


function erd_write_raw_csv(path,rows,completed,runs)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,seed,method,data_file,samples,features,classes,' ...
    'labeled_count,alpha,beta,parameter_source,iterations,unlabeled_ACC,' ...
    'unlabeled_NMI,all_ACC,all_NMI,labeled_ACC,validation_ACC,objective_initial,' ...
    'objective_final,objective_nonincrease,stop_residual,converged,' ...
    'graph_seconds,run_seconds\n']);
for d=1:size(completed,1)
    for r=1:runs
        if completed(d,r)
            idx = (d-1)*runs+r;
            S = rows(idx);
            fprintf(fid,['%s,%d,%s,%s,%d,%d,%d,%d,%.15g,%.15g,%s,%d,' ...
                '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d,' ...
                '%.15g,%d,%.6f,%.6f\n'], ...
                S.dataset,S.seed,S.method,S.data_file,S.samples,S.features, ...
                S.classes,S.labeled_count,S.alpha,S.beta,S.parameter_source, ...
                S.iterations,S.unlabeled_ACC,S.unlabeled_NMI,S.all_ACC, ...
                S.all_NMI,S.labeled_ACC,S.validation_ACC,S.objective_initial, ...
                S.objective_final,S.objective_nonincrease,S.stop_residual, ...
                S.converged,S.graph_seconds,S.run_seconds);
        end
    end
end
fclose(fid);
end


function rows = erd_summarize(rawRows,D,runs)
template = struct('dataset','','method','SNMFWLP','runs',runs, ...
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


function erd_write_decision(path,runs,rawRows,summaryRows,D,maxIterations)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'SNMFWLP COMMON-STOP AUDIT\n');
fprintf(fid,'runs=%d\n',runs);
fprintf(fid,'seeds=20260617-%d\n',20260617+runs-1);
fprintf(fid,'label_fraction=0.10 with at least two labeled samples per class\n');
fprintf(fid,'evaluation=unlabeled samples only\n');
fprintf(fid,'assignment=row-wise argmax; no K-means\n');
fprintf(fid,'maximum_iterations=%d\n',maxIterations);
fprintf(fid,'implementation=public SNMFWLP updates plus disclosed integration wrapper\n\n');
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
