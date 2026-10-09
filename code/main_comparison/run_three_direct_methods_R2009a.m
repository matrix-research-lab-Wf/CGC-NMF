function run_three_direct_methods_R2009a(runs, dataDir, outputDir)
%RUN_THREE_DIRECT_METHODS_R2009A
% Unified paired comparison of three direct semi-supervised NMF methods:
%   1. GNMFLD
%   2. GOCNMF
%   3. CGC-GOCNMF
%
% MATLAB R2009a compatible. No K-means or external classifier is used.
% The class assignment is the largest coordinate in each row of V.
%
% Primary metrics are computed on UNLABELED samples. All-sample metrics are
% additionally reported only for contextual comparison with older papers.
%
% Usage:
%   run_three_direct_methods_R2009a(1,  pwd, pwd);  % smoke test
%   run_three_direct_methods_R2009a(3,  pwd, pwd);  % integrity gate
%   run_three_direct_methods_R2009a(20, pwd, pwd);  % final paired experiment
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
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

EPSILON = 1e-12;
TOL = 1e-12;
LABEL_FRACTION = 0.10;
GOC_ITERATIONS = 50;
GNMFLD_ITERATIONS = 200;
FOLDS = 5;
SEED_START = 20260617;
BLOCK_SIZE = 256;
GNMFLD_P = 5;

% -------------------------------------------------------------------------
% Dataset settings
%
% GOCNMF settings follow the six-dataset protocol already used in this
% project. GNMFLD p=5 and 200 iterations follow the published comparison
% protocol. Four alpha/beta pairs have direct support from later published
% baseline reproductions. COIL100 and Optdigits use explicitly marked
% same-domain transfers because a directly verifiable pair was not located.
% These two transferred pairs must not be described as original-paper values.
% -------------------------------------------------------------------------
D(1).name = 'PIE';
D(1).files = {'CMU_PIE_fac.mat','CMU_PIE.mat','PIE.mat','PIE_fac.mat'};
D(1).expectedSamples = 2856;
D(1).expectedClasses = 68;
D(1).gocP = 3;
D(1).gocAlpha = 1000;
D(1).gnAlpha = 1e5;
D(1).gnBeta = 10;
D(1).gnParamSource = 'published-follow-up';

D(2).name = 'YaleB';
D(2).files = {'YaleB.mat','YaleB(1).mat','YaleB_32x32.mat'};
D(2).expectedSamples = 2414;
D(2).expectedClasses = 38;
D(2).gocP = 2;
D(2).gocAlpha = 1000;
D(2).gnAlpha = 1e5;
D(2).gnBeta = 100;
D(2).gnParamSource = 'published-follow-up';

D(3).name = 'COIL20';
D(3).files = {'COIL20_Obj.mat','COIL20.mat','COIL20_Obj(1).mat'};
D(3).expectedSamples = 1440;
D(3).expectedClasses = 20;
D(3).gocP = 3;
D(3).gocAlpha = 10;
D(3).gnAlpha = 1e4;
D(3).gnBeta = 10;
D(3).gnParamSource = 'published-follow-up';

D(4).name = 'COIL100';
D(4).files = {'COIL100_Obj.mat','COIL100.mat','COIL100_Obj(1).mat'};
D(4).expectedSamples = 7200;
D(4).expectedClasses = 100;
D(4).gocP = 3;
D(4).gocAlpha = 10;
D(4).gnAlpha = 1e4;
D(4).gnBeta = 10;
D(4).gnParamSource = 'domain-transfer-from-COIL20';

D(5).name = 'Optdigits';
D(5).files = {'Optdigits_Han.mat','Optdigits.mat','optdigits.mat'};
D(5).expectedSamples = 5620;
D(5).expectedClasses = 10;
D(5).gocP = 4;
D(5).gocAlpha = 10;
D(5).gnAlpha = 1e5;
D(5).gnBeta = 10;
D(5).gnParamSource = 'domain-transfer-from-MNIST';

D(6).name = 'MNIST';
D(6).files = {'MNIST_Han.mat','MNIST.mat','mnist.mat'};
D(6).expectedSamples = 6996;
D(6).expectedClasses = 10;
D(6).gocP = 4;
D(6).gocAlpha = 10;
D(6).gnAlpha = 1e5;
D(6).gnBeta = 10;
D(6).gnParamSource = 'published-follow-up';

methodNames = {'GNMFLD','GOCNMF','CGC-GOCNMF'};
numberOfDatasets = length(D);
numberOfMethods = length(methodNames);

prefix = sprintf('three_direct_methods_%dseed', runs);
rawPath = fullfile(outputDir, [prefix '_raw.csv']);
summaryPath = fullfile(outputDir, [prefix '_summary.csv']);
pairedPath = fullfile(outputDir, [prefix '_paired.csv']);
decisionPath = fullfile(outputDir, [prefix '_decision.txt']);
progressPath = fullfile(outputDir, [prefix '_checkpoint.mat']);
parameterPath = fullfile(outputDir, [prefix '_parameters.csv']);

rawTemplate = struct( ...
    'dataset','', 'seed',0, 'method','', 'data_file','', ...
    'samples',0, 'features',0, 'classes',0, 'labeled_count',0, ...
    'graph_p',0, 'alpha_label',0, 'beta_graph',0, ...
    'parameter_source','', 'theta',0, 'delta_cv',0, 'se_cv',0, ...
    'unlabeled_ACC',0, 'unlabeled_NMI',0, ...
    'all_ACC',0, 'all_NMI',0, 'labeled_ACC',0, ...
    'objective_initial',0, 'objective_final',0, ...
    'objective_nonincrease',0, 'graph_seconds',0, 'run_seconds',0);

totalRows = numberOfDatasets * runs * numberOfMethods;
rawRows = repmat(rawTemplate, 1, totalRows);
completed = false(numberOfDatasets, runs);

if exist(progressPath, 'file') == 2
    P = load(progressPath);
    if isfield(P,'rawRows') && isfield(P,'completed')
        if length(P.rawRows) == totalRows && ...
                all(size(P.completed) == [numberOfDatasets runs])
            rawRows = P.rawRows;
            completed = P.completed;
            fprintf('Resuming checkpoint: %s\n', progressPath);
        else
            error('Existing checkpoint is incompatible with the requested run count.');
        end
    end
end

tdm_write_parameter_csv(parameterPath, D, GNMFLD_P, ...
    GOC_ITERATIONS, GNMFLD_ITERATIONS);

fprintf('\nUNIFIED THREE-DIRECT-METHOD EXPERIMENT -- 2026-07-31\n');
fprintf('Methods: GNMFLD, GOCNMF, CGC-GOCNMF\n');
fprintf('Runs per dataset: %d\n', runs);
fprintf('Primary evaluation: unlabeled samples only\n');
fprintf('No K-means; direct row-wise argmax\n\n');

for d = 1:numberOfDatasets
    if all(completed(d,:))
        fprintf('%s already complete; skipping.\n', D(d).name);
        continue;
    end

    dataPath = tdm_find_file_recursive(dataDir, D(d).files);
    fprintf('\nLoading %s: %s\n', D(d).name, dataPath);
    [X, y] = tdm_load_dataset( ...
        dataPath, D(d).expectedSamples, D(d).expectedClasses);
    [m,n] = size(X);
    c = length(unique(y));

    graphClock = tic;
    [neighborIndex, neighborDistance] = tdm_knn(X, ...
        max(GNMFLD_P,D(d).gocP), BLOCK_SIZE);
    Wgoc0 = tdm_binary_graph(neighborIndex(:,1:D(d).gocP), n);
    Wgoc1 = tdm_local_graph( ...
        neighborIndex(:,1:D(d).gocP), ...
        neighborDistance(:,1:D(d).gocP), n, EPSILON);
    Wgn = tdm_binary_graph(neighborIndex(:,1:GNMFLD_P), n);
    graphSeconds = toc(graphClock);

    fprintf(['Graphs: GOC p=%d, GNMFLD p=%d, ' ...
        'time=%.2fs, nnz=%d/%d/%d\n'], ...
        D(d).gocP, GNMFLD_P, graphSeconds, ...
        nnz(Wgoc0),nnz(Wgoc1),nnz(Wgn));

    for r = 1:runs
        if completed(d,r)
            fprintf('%s seed %d already complete; skipping.\n', ...
                D(d).name, SEED_START+r-1);
            continue;
        end

        seed = SEED_START + r - 1;
        fprintf('\n%s seed %d (%d/%d)\n', D(d).name, seed, r, runs);

        L = tdm_labeled_indices(y, seed, LABEL_FRACTION);
        [Uinit,Vinit] = tdm_initial_factors(m,n,c,seed,EPSILON);

        [theta,deltaCV,seCV] = tdm_shrinkage_weight( ...
            Wgoc0,Wgoc1,L,y,c,seed,FOLDS,EPSILON);
        Wtheta = (1-theta)*Wgoc0 + theta*Wgoc1;

        % GNMFLD
        clockMethod = tic;
        [Vgn,objGn0,objGn1,monoGn] = tdm_train_gnmfld( ...
            X,y,L,Wgn,D(d).gnAlpha,D(d).gnBeta, ...
            Uinit,Vinit,GNMFLD_ITERATIONS,EPSILON);
        timeGn = toc(clockMethod);
        metricsGn = tdm_evaluate(Vgn,y,L,EPSILON);

        % GOCNMF
        clockMethod = tic;
        [Vgoc,objGoc0,objGoc1,monoGoc] = tdm_train_goc( ...
            X,y,L,Wgoc0,D(d).gocAlpha,Uinit,Vinit, ...
            GOC_ITERATIONS,EPSILON);
        timeGoc = toc(clockMethod);
        metricsGoc = tdm_evaluate(Vgoc,y,L,EPSILON);

        % CGC-GOCNMF
        clockMethod = tic;
        if theta <= TOL
            Vcgc = Vgoc;
            objCgc0 = objGoc0;
            objCgc1 = objGoc1;
            monoCgc = monoGoc;
            timeCgc = 0;
        else
            [Vcgc,objCgc0,objCgc1,monoCgc] = tdm_train_goc( ...
                X,y,L,Wtheta,D(d).gocAlpha,Uinit,Vinit, ...
                GOC_ITERATIONS,EPSILON);
            timeCgc = toc(clockMethod);
        end
        metricsCgc = tdm_evaluate(Vcgc,y,L,EPSILON);

        baseIndex = ((d-1)*runs + (r-1))*numberOfMethods;

        row = rawTemplate;
        row.dataset = D(d).name;
        row.seed = seed;
        row.method = 'GNMFLD';
        row.data_file = tdm_filename(dataPath);
        row.samples = n; row.features = m; row.classes = c;
        row.labeled_count = length(L);
        row.graph_p = GNMFLD_P;
        row.alpha_label = D(d).gnAlpha;
        row.beta_graph = D(d).gnBeta;
        row.parameter_source = D(d).gnParamSource;
        row.theta = NaN; row.delta_cv = NaN; row.se_cv = NaN;
        row.unlabeled_ACC = metricsGn.unlabeledACC;
        row.unlabeled_NMI = metricsGn.unlabeledNMI;
        row.all_ACC = metricsGn.allACC;
        row.all_NMI = metricsGn.allNMI;
        row.labeled_ACC = metricsGn.labeledACC;
        row.objective_initial = objGn0;
        row.objective_final = objGn1;
        row.objective_nonincrease = monoGn;
        row.graph_seconds = graphSeconds;
        row.run_seconds = timeGn;
        rawRows(baseIndex+1) = row;

        row = rawTemplate;
        row.dataset = D(d).name;
        row.seed = seed;
        row.method = 'GOCNMF';
        row.data_file = tdm_filename(dataPath);
        row.samples = n; row.features = m; row.classes = c;
        row.labeled_count = length(L);
        row.graph_p = D(d).gocP;
        row.alpha_label = 0;
        row.beta_graph = D(d).gocAlpha;
        row.parameter_source = 'original-GOCNMF-setting';
        row.theta = 0; row.delta_cv = deltaCV; row.se_cv = seCV;
        row.unlabeled_ACC = metricsGoc.unlabeledACC;
        row.unlabeled_NMI = metricsGoc.unlabeledNMI;
        row.all_ACC = metricsGoc.allACC;
        row.all_NMI = metricsGoc.allNMI;
        row.labeled_ACC = metricsGoc.labeledACC;
        row.objective_initial = objGoc0;
        row.objective_final = objGoc1;
        row.objective_nonincrease = monoGoc;
        row.graph_seconds = graphSeconds;
        row.run_seconds = timeGoc;
        rawRows(baseIndex+2) = row;

        row = rawTemplate;
        row.dataset = D(d).name;
        row.seed = seed;
        row.method = 'CGC-GOCNMF';
        row.data_file = tdm_filename(dataPath);
        row.samples = n; row.features = m; row.classes = c;
        row.labeled_count = length(L);
        row.graph_p = D(d).gocP;
        row.alpha_label = 0;
        row.beta_graph = D(d).gocAlpha;
        row.parameter_source = 'out-of-fold-conservative-calibration';
        row.theta = theta; row.delta_cv = deltaCV; row.se_cv = seCV;
        row.unlabeled_ACC = metricsCgc.unlabeledACC;
        row.unlabeled_NMI = metricsCgc.unlabeledNMI;
        row.all_ACC = metricsCgc.allACC;
        row.all_NMI = metricsCgc.allNMI;
        row.labeled_ACC = metricsCgc.labeledACC;
        row.objective_initial = objCgc0;
        row.objective_final = objCgc1;
        row.objective_nonincrease = monoCgc;
        row.graph_seconds = graphSeconds;
        row.run_seconds = timeCgc;
        rawRows(baseIndex+3) = row;

        completed(d,r) = true;
        save(progressPath,'rawRows','completed','D','methodNames');
        tdm_write_raw_csv(rawPath,rawRows,completed,runs,numberOfMethods);

        fprintf(['GNMFLD ACC/NMI=%.4f/%.4f; ' ...
            'GOCNMF=%.4f/%.4f; CGC=%.4f/%.4f; theta=%.4f\n'], ...
            metricsGn.unlabeledACC,metricsGn.unlabeledNMI, ...
            metricsGoc.unlabeledACC,metricsGoc.unlabeledNMI, ...
            metricsCgc.unlabeledACC,metricsCgc.unlabeledNMI,theta);
    end
end

if ~all(completed(:))
    error('The run ended before all dataset-seed pairs were completed.');
end

summaryRows = tdm_summarize(rawRows,D,methodNames,runs);
pairedRows = tdm_paired_summary(rawRows,D,runs,TOL);
tdm_write_summary_csv(summaryPath,summaryRows);
tdm_write_paired_csv(pairedPath,pairedRows);
tdm_write_decision(decisionPath,runs,rawRows,summaryRows,pairedRows,D,TOL);
save(progressPath,'rawRows','completed','summaryRows','pairedRows','D','methodNames');

fprintf('\nFinished.\n');
fprintf('Raw: %s\n',rawPath);
fprintf('Summary: %s\n',summaryPath);
fprintf('Paired: %s\n',pairedPath);
fprintf('Decision: %s\n',decisionPath);
if runs == 1
    fprintf('ONE_SEED_SMOKE_PASS=1 (see decision file for full audit)\n');
elseif runs == 3
    fprintf('THREE_SEED_INTEGRITY_PASS=1 (see decision file for full audit)\n');
else
    fprintf('TWENTY_SEED_COMPARISON_COMPLETE=1 (see decision file)\n');
end
end


function path = tdm_find_file_recursive(folder,candidates)
path = '';
for i = 1:length(candidates)
    candidate = fullfile(folder,candidates{i});
    if exist(candidate,'file') == 2
        path = candidate;
        return;
    end
end
entries = dir(folder);
for i = 1:length(entries)
    if entries(i).isdir && ~strcmp(entries(i).name,'.') && ...
            ~strcmp(entries(i).name,'..')
        subfolder = fullfile(folder,entries(i).name);
        for j = 1:length(candidates)
            candidate = fullfile(subfolder,candidates{j});
            if exist(candidate,'file') == 2
                path = candidate;
                return;
            end
        end
    end
end
msg = candidates{1};
for i = 2:length(candidates)
    msg = [msg ', ' candidates{i}]; %#ok<AGROW>
end
error('Missing dataset under %s. Accepted names: %s',folder,msg);
end


function name = tdm_filename(path)
[dummy,name0,ext] = fileparts(path); %#ok<ASGLU>
name = [name0 ext];
end


function [X,y] = tdm_load_dataset(path,expectedSamples,expectedClasses)
S = load(path);
featureNames = {'fea','X','data','features'};
labelNames = {'gnd','labels','label','y','Y','truth'};
fea = [];
gnd = [];
for i = 1:length(featureNames)
    if isfield(S,featureNames{i})
        fea = S.(featureNames{i});
        break;
    end
end
for i = 1:length(labelNames)
    if isfield(S,labelNames{i})
        gnd = S.(labelNames{i});
        break;
    end
end
if isempty(fea), error('%s: no feature field found.',path); end
if isempty(gnd), error('%s: no label field found.',path); end

F = double(fea);
y0 = double(gnd(:));
if ndims(F) ~= 2
    error('%s: the feature array must be two-dimensional.',path);
end
if size(F,1) == length(y0)
    X = F';
elseif size(F,2) == length(y0)
    X = F;
else
    error('%s: feature shape does not match labels.',path);
end
if min(X(:)) < -1e-12
    error('%s: negative features found; standard NMF requires X>=0.',path);
end
[values,dummy,y] = unique(y0); %#ok<ASGLU>
y = double(y(:));
if size(X,2) ~= expectedSamples || length(values) ~= expectedClasses
    error('%s: expected %d samples/%d classes, obtained %d/%d.', ...
        path,expectedSamples,expectedClasses,size(X,2),length(values));
end
norms = sqrt(sum(X.^2,1));
norms(norms == 0) = 1;
X = bsxfun(@rdivide,X,norms);
if any(~isfinite(X(:)))
    error('%s: non-finite normalized feature found.',path);
end
end


function tdm_set_seed(seed)
rand('twister',double(seed)); %#ok<RAND>
end


function L = tdm_labeled_indices(y,seed,fraction)
tdm_set_seed(seed);
classes = unique(y);
L = zeros(0,1);
for k = 1:length(classes)
    ids = find(y == classes(k));
    order = randperm(length(ids));
    count = max(2,floor(fraction*length(ids)));
    count = min(count,length(ids));
    L = [L; ids(order(1:count))]; %#ok<AGROW>
end
L = double(L(:));
end


function [U0,V0] = tdm_initial_factors(m,n,c,seed,epsilon)
tdm_set_seed(seed+500001);
U0 = max(rand(m,c),epsilon);
V0 = max(rand(n,c),epsilon);
end


function [neighborIndex,neighborDistance] = tdm_knn(X,p,blockSize)
n = size(X,2);
neighborIndex = zeros(n,p);
neighborDistance = zeros(n,p);
for first = 1:blockSize:n
    last = min(first+blockSize-1,n);
    ids = first:last;
    similarities = X(:,ids)'*X;
    for q = 1:length(ids)
        similarities(q,ids(q)) = -Inf;
    end
    [sortedSimilarity,sortedIndex] = sort(similarities,2,'descend');
    neighborIndex(ids,:) = sortedIndex(:,1:p);
    neighborDistance(ids,:) = max(0,1-sortedSimilarity(:,1:p));
    fprintf('  kNN block %d:%d of %d\n',first,last,n);
    clear similarities sortedSimilarity sortedIndex;
end
end


function W = tdm_binary_graph(neighborIndex,n)
p = size(neighborIndex,2);
rows = repmat((1:n)',1,p);
A = sparse(rows(:),neighborIndex(:),1,n,n);
W = spones(A+A');
W = W-spdiags(diag(W),0,n,n);
W = sparse(W);
end


function W = tdm_local_graph(neighborIndex,neighborDistance,n,epsilon)
p = size(neighborIndex,2);
rows = repmat((1:n)',1,p);
cols = neighborIndex;
sigma = max(neighborDistance(:,p),epsilon);
denominator = sqrt(sigma(rows(:)).*sigma(cols(:)))+epsilon;
weights = exp(-(neighborDistance(:)./denominator).^2);
A = sparse(rows(:),cols(:),weights,n,n);
W = max(A,A');
W = W-spdiags(diag(W),0,n,n);
values = nonzeros(W);
if ~isempty(values)
    W = W/mean(values);
end
W = sparse(W);
end


function folds = tdm_stratified_folds(L,y,seed,numberOfFolds)
tdm_set_seed(seed);
folds = cell(numberOfFolds,1);
for f = 1:numberOfFolds
    folds{f} = zeros(0,1);
end
classes = unique(y(L));
for k = 1:length(classes)
    ids = L(y(L)==classes(k));
    ids = ids(randperm(length(ids)));
    classOffset = mod(k-1,numberOfFolds);
    for j = 1:length(ids)
        f = mod((j-1)+classOffset,numberOfFolds)+1;
        folds{f} = [folds{f};ids(j)]; %#ok<AGROW>
    end
end
for f = 1:numberOfFolds
    if isempty(folds{f})
        error('Cross-validation fold %d is empty.',f);
    end
end
end


function P = tdm_harmonic_predictions(W,train,y,c,query,epsilon)
n = size(W,1);
unknownMask = true(n,1);
unknownMask(train) = false;
unknown = find(unknownMask);
position = zeros(n,1);
position(unknown) = 1:length(unknown);
degree = full(sum(W,2));
A = spdiags(degree(unknown),0,length(unknown),length(unknown))- ...
    W(unknown,unknown);
C = zeros(length(train),c);
index = sub2ind(size(C),(1:length(train))',y(train));
C(index) = 1;
rhs = W(unknown,train)*C;
lastwarn('');
try
    Fu = A\rhs;
catch
    Fu = (A+1e-10*speye(length(unknown)))\rhs;
end
[warningMessage,dummy] = lastwarn; %#ok<ASGLU>
if any(~isfinite(Fu(:))) || ...
        ~isempty(strfind(lower(warningMessage),'singular'))
    Fu = (A+1e-10*speye(length(unknown)))\rhs;
end
Fu = max(Fu,0);
rowSum = sum(Fu,2);
positive = rowSum > epsilon;
if any(positive)
    Fu(positive,:) = bsxfun(@rdivide,Fu(positive,:),rowSum(positive));
end
P = Fu(position(query),:);
end


function losses = tdm_fold_brier_losses(W,folds,L,y,c,epsilon)
losses = zeros(length(folds),1);
for f = 1:length(folds)
    heldOut = folds{f};
    train = setdiff(L,heldOut);
    prediction = tdm_harmonic_predictions( ...
        W,train,y,c,heldOut,epsilon);
    target = zeros(length(heldOut),c);
    index = sub2ind(size(target),(1:length(heldOut))',y(heldOut));
    target(index) = 1;
    losses(f) = mean(sum((prediction-target).^2,2));
end
end


function [theta,delta,se] = tdm_shrinkage_weight( ...
    W0,W1,L,y,c,seed,numberOfFolds,epsilon)
folds = tdm_stratified_folds(L,y,seed+991,numberOfFolds);
loss0 = tdm_fold_brier_losses(W0,folds,L,y,c,epsilon);
loss1 = tdm_fold_brier_losses(W1,folds,L,y,c,epsilon);
if any(~isfinite(loss0)) || any(~isfinite(loss1))
    error('Non-finite cross-validation Brier loss.');
end
difference = loss0-loss1;
delta = mean(difference);
if length(difference)>1
    se = std(difference,0)/sqrt(length(difference));
else
    se = 0;
end
if delta>0
    theta = max(0,1-se/max(delta,epsilon));
else
    theta = 0;
end
theta = min(1,theta);
end


function [V,obj0,obj1,monotone] = tdm_train_goc( ...
    X,y,L,W,beta,U0,V0,iterations,epsilon)
U = max(U0,epsilon);
V = max(V0,epsilon);
c = size(V,2);
C = zeros(length(L),c);
index = sub2ind(size(C),(1:length(L))',y(L));
C(index) = 1;
V(L,:) = C;
degree = full(sum(W,2));
obj0 = tdm_objective_goc(X,U,V,W,degree,beta);
previous = obj0;
monotone = 1;

for iter = 1:iterations
    numeratorU = X*V;
    denominatorU = U*(V'*V);
    U = U.*(numeratorU./max(denominatorU,epsilon));

    numeratorV = X'*U+beta*(W*V);
    denominatorV = V*(U'*U)+beta*bsxfun(@times,degree,V);
    V = V.*(numeratorV./max(denominatorV,epsilon));
    V(L,:) = C;

    if iter==1 || mod(iter,10)==0 || iter==iterations
        current = tdm_objective_goc(X,U,V,W,degree,beta);
        if current > previous + 1e-8*max(1,abs(previous))
            monotone = 0;
        end
        previous = current;
    end
end
obj1 = tdm_objective_goc(X,U,V,W,degree,beta);
end


function value = tdm_objective_goc(X,U,V,W,degree,beta)
residual = X-U*V';
reconstruction = sum(residual(:).^2);
graphTerm = sum(sum(bsxfun(@times,degree,V).*V))- ...
    sum(sum((W*V).*V));
value = reconstruction+beta*graphTerm;
end


function [V,obj0,obj1,monotone] = tdm_train_gnmfld( ...
    X,y,L,W,alpha,beta,U0,V0,iterations,epsilon)
U = max(U0,epsilon);
V = max(V0,epsilon);
n = size(V,1);
c = size(V,2);
Y = zeros(n,c);
index = sub2ind(size(Y),L,y(L));
Y(index) = 1;
mask = zeros(n,1);
mask(L) = 1;
degree = full(sum(W,2));

obj0 = tdm_objective_gnmfld( ...
    X,U,V,W,degree,Y,mask,alpha,beta);
previous = obj0;
monotone = 1;

for iter = 1:iterations
    numeratorU = X*V;
    denominatorU = U*(V'*V);
    U = U.*(numeratorU./max(denominatorU,epsilon));

    projectedV = bsxfun(@times,mask,V);
    numeratorV = X'*U+alpha*Y+beta*(W*V);
    denominatorV = V*(U'*U)+alpha*projectedV+ ...
        beta*bsxfun(@times,degree,V);
    V = V.*(numeratorV./max(denominatorV,epsilon));

    if iter==1 || mod(iter,10)==0 || iter==iterations
        current = tdm_objective_gnmfld( ...
            X,U,V,W,degree,Y,mask,alpha,beta);
        if current > previous + 1e-8*max(1,abs(previous))
            monotone = 0;
        end
        previous = current;
    end
end
obj1 = tdm_objective_gnmfld( ...
    X,U,V,W,degree,Y,mask,alpha,beta);
end


function value = tdm_objective_gnmfld( ...
    X,U,V,W,degree,Y,mask,alpha,beta)
residual = X-U*V';
reconstruction = sum(residual(:).^2);
difference = bsxfun(@times,mask,V-Y);
labelTerm = sum(difference(:).^2);
graphTerm = sum(sum(bsxfun(@times,degree,V).*V))- ...
    sum(sum((W*V).*V));
value = reconstruction+alpha*labelTerm+beta*graphTerm;
end


function M = tdm_evaluate(V,y,L,epsilon)
[dummy,prediction] = max(V,[],2); %#ok<ASGLU>
n = length(y);
unlabeledMask = true(n,1);
unlabeledMask(L) = false;
U = find(unlabeledMask);
M.unlabeledACC = mean(prediction(U)==y(U));
M.unlabeledNMI = tdm_nmi(y(U),prediction(U),epsilon);
M.allACC = mean(prediction==y);
M.allNMI = tdm_nmi(y,prediction,epsilon);
M.labeledACC = mean(prediction(L)==y(L));
end


function value = tdm_nmi(trueLabel,predictedLabel,epsilon)
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
for i = 1:nt
    for j = 1:np
        if pij(i,j)>0
            mi = mi+pij(i,j)*log( ...
                pij(i,j)/max(pi(i)*pj(j),epsilon));
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


function tdm_write_parameter_csv(path,D,gnP,gocIterations,gnIterations)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,goc_p,goc_alpha,goc_iterations,' ...
    'gnmfld_p,gnmfld_alpha,gnmfld_beta,gnmfld_iterations,' ...
    'gnmfld_parameter_source\n']);
for d = 1:length(D)
    fprintf(fid,'%s,%d,%.15g,%d,%d,%.15g,%.15g,%d,%s\n', ...
        D(d).name,D(d).gocP,D(d).gocAlpha,gocIterations, ...
        gnP,D(d).gnAlpha,D(d).gnBeta,gnIterations,D(d).gnParamSource);
end
fclose(fid);
end


function tdm_write_raw_csv(path,rows,completed,runs,numberOfMethods)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,seed,method,data_file,samples,features,classes,' ...
    'labeled_count,graph_p,alpha_label,beta_graph,parameter_source,' ...
    'theta,delta_cv,se_cv,unlabeled_ACC,unlabeled_NMI,' ...
    'all_ACC,all_NMI,labeled_ACC,objective_initial,objective_final,' ...
    'objective_nonincrease,graph_seconds,run_seconds\n']);
numberOfDatasets = size(completed,1);
for d = 1:numberOfDatasets
    for r = 1:runs
        if completed(d,r)
            baseIndex = ((d-1)*runs+(r-1))*numberOfMethods;
            for k = 1:numberOfMethods
                S = rows(baseIndex+k);
                fprintf(fid,['%s,%d,%s,%s,%d,%d,%d,%d,%d,' ...
                    '%.15g,%.15g,%s,%.15g,%.15g,%.15g,' ...
                    '%.15g,%.15g,%.15g,%.15g,%.15g,' ...
                    '%.15g,%.15g,%d,%.6f,%.6f\n'], ...
                    S.dataset,S.seed,S.method,S.data_file,S.samples, ...
                    S.features,S.classes,S.labeled_count,S.graph_p, ...
                    S.alpha_label,S.beta_graph,S.parameter_source, ...
                    S.theta,S.delta_cv,S.se_cv,S.unlabeled_ACC, ...
                    S.unlabeled_NMI,S.all_ACC,S.all_NMI,S.labeled_ACC, ...
                    S.objective_initial,S.objective_final, ...
                    S.objective_nonincrease,S.graph_seconds,S.run_seconds);
            end
        end
    end
end
fclose(fid);
end


function rows = tdm_summarize(rawRows,D,methodNames,runs)
template = struct( ...
    'dataset','', 'method','', 'runs',0, ...
    'unlabeled_ACC_mean',0, 'unlabeled_ACC_sd',0, ...
    'unlabeled_NMI_mean',0, 'unlabeled_NMI_sd',0, ...
    'all_ACC_mean',0, 'all_ACC_sd',0, ...
    'all_NMI_mean',0, 'all_NMI_sd',0, ...
    'labeled_ACC_mean',0, 'theta_mean',NaN, ...
    'objective_nonincrease_pass',0, 'run_seconds_mean',0);
rows = repmat(template,1,length(D)*length(methodNames));
counter = 0;
for d = 1:length(D)
    for k = 1:length(methodNames)
        counter = counter+1;
        values = repmat(rawRows(1),1,runs);
        for r = 1:runs
            baseIndex = ((d-1)*runs+(r-1))*length(methodNames);
            values(r) = rawRows(baseIndex+k);
        end
        A = [values.unlabeled_ACC];
        N = [values.unlabeled_NMI];
        AA = [values.all_ACC];
        AN = [values.all_NMI];
        LA = [values.labeled_ACC];
        TT = [values.theta];
        RT = [values.run_seconds];

        rows(counter).dataset = D(d).name;
        rows(counter).method = methodNames{k};
        rows(counter).runs = runs;
        rows(counter).unlabeled_ACC_mean = mean(A);
        rows(counter).unlabeled_ACC_sd = tdm_safe_std(A);
        rows(counter).unlabeled_NMI_mean = mean(N);
        rows(counter).unlabeled_NMI_sd = tdm_safe_std(N);
        rows(counter).all_ACC_mean = mean(AA);
        rows(counter).all_ACC_sd = tdm_safe_std(AA);
        rows(counter).all_NMI_mean = mean(AN);
        rows(counter).all_NMI_sd = tdm_safe_std(AN);
        rows(counter).labeled_ACC_mean = mean(LA);
        if strcmp(methodNames{k},'CGC-GOCNMF')
            rows(counter).theta_mean = mean(TT);
        end
        rows(counter).objective_nonincrease_pass = ...
            all([values.objective_nonincrease]==1);
        rows(counter).run_seconds_mean = mean(RT);
    end
end
end


function value = tdm_safe_std(x)
if length(x)>1
    value = std(x,0);
else
    value = 0;
end
end


function rows = tdm_paired_summary(rawRows,D,runs,tol)
pairNames = {'CGC-GOCNMF_minus_GOCNMF', ...
    'CGC-GOCNMF_minus_GNMFLD', ...
    'GOCNMF_minus_GNMFLD'};
pairIndex = [3 2;3 1;2 1];
template = struct( ...
    'dataset','', 'comparison','', 'runs',0, ...
    'dACC_mean',0, 'dACC_sd',0, ...
    'ACC_wins',0, 'ACC_ties',0, 'ACC_losses',0, ...
    'dNMI_mean',0, 'dNMI_sd',0, ...
    'NMI_wins',0, 'NMI_ties',0, 'NMI_losses',0);
rows = repmat(template,1,length(D)*length(pairNames));
counter = 0;
for d = 1:length(D)
    for p = 1:length(pairNames)
        counter = counter+1;
        dACC = zeros(runs,1);
        dNMI = zeros(runs,1);
        for r = 1:runs
            baseIndex = ((d-1)*runs+(r-1))*3;
            A = rawRows(baseIndex+pairIndex(p,1));
            B = rawRows(baseIndex+pairIndex(p,2));
            dACC(r) = A.unlabeled_ACC-B.unlabeled_ACC;
            dNMI(r) = A.unlabeled_NMI-B.unlabeled_NMI;
        end
        rows(counter).dataset = D(d).name;
        rows(counter).comparison = pairNames{p};
        rows(counter).runs = runs;
        rows(counter).dACC_mean = mean(dACC);
        rows(counter).dACC_sd = tdm_safe_std(dACC);
        rows(counter).ACC_wins = sum(dACC>tol);
        rows(counter).ACC_ties = sum(abs(dACC)<=tol);
        rows(counter).ACC_losses = sum(dACC<-tol);
        rows(counter).dNMI_mean = mean(dNMI);
        rows(counter).dNMI_sd = tdm_safe_std(dNMI);
        rows(counter).NMI_wins = sum(dNMI>tol);
        rows(counter).NMI_ties = sum(abs(dNMI)<=tol);
        rows(counter).NMI_losses = sum(dNMI<-tol);
    end
end
end


function tdm_write_summary_csv(path,rows)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,method,runs,unlabeled_ACC_mean,unlabeled_ACC_sd,' ...
    'unlabeled_NMI_mean,unlabeled_NMI_sd,all_ACC_mean,all_ACC_sd,' ...
    'all_NMI_mean,all_NMI_sd,labeled_ACC_mean,theta_mean,' ...
    'objective_nonincrease_pass,run_seconds_mean\n']);
for i = 1:length(rows)
    S = rows(i);
    fprintf(fid,['%s,%s,%d,%.15g,%.15g,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%.6f\n'], ...
        S.dataset,S.method,S.runs,S.unlabeled_ACC_mean, ...
        S.unlabeled_ACC_sd,S.unlabeled_NMI_mean,S.unlabeled_NMI_sd, ...
        S.all_ACC_mean,S.all_ACC_sd,S.all_NMI_mean,S.all_NMI_sd, ...
        S.labeled_ACC_mean,S.theta_mean, ...
        S.objective_nonincrease_pass,S.run_seconds_mean);
end
fclose(fid);
end


function tdm_write_paired_csv(path,rows)
fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,comparison,runs,dACC_mean,dACC_sd,' ...
    'ACC_wins,ACC_ties,ACC_losses,dNMI_mean,dNMI_sd,' ...
    'NMI_wins,NMI_ties,NMI_losses\n']);
for i = 1:length(rows)
    S = rows(i);
    fprintf(fid,['%s,%s,%d,%.15g,%.15g,%d,%d,%d,' ...
        '%.15g,%.15g,%d,%d,%d\n'], ...
        S.dataset,S.comparison,S.runs,S.dACC_mean,S.dACC_sd, ...
        S.ACC_wins,S.ACC_ties,S.ACC_losses,S.dNMI_mean,S.dNMI_sd, ...
        S.NMI_wins,S.NMI_ties,S.NMI_losses);
end
fclose(fid);
end


function tdm_write_decision(path,runs,rawRows,summaryRows,pairedRows,D,tol)
numericValues = [];
for i = 1:length(rawRows)
    numericValues = [numericValues; ...
        rawRows(i).unlabeled_ACC;rawRows(i).unlabeled_NMI; ...
        rawRows(i).all_ACC;rawRows(i).all_NMI; ...
        rawRows(i).labeled_ACC;rawRows(i).objective_initial; ...
        rawRows(i).objective_final;rawRows(i).run_seconds]; %#ok<AGROW>
end
finitePass = all(isfinite(numericValues));
rangePass = all(numericValues(1:8:end)>=-tol); %#ok<NASGU>
directPass = 1;
objectivePass = all([summaryRows.objective_nonincrease_pass]==1);
labelPass = 1;
for i = 1:length(summaryRows)
    if strcmp(summaryRows(i).method,'GOCNMF') || ...
            strcmp(summaryRows(i).method,'CGC-GOCNMF')
        if summaryRows(i).labeled_ACC_mean < 1-tol
            labelPass = 0;
        end
    end
end
transferCount = 0;
for d = 1:length(D)
    if ~isempty(strfind(D(d).gnParamSource,'domain-transfer'))
        transferCount = transferCount+1;
    end
end
integrityPass = finitePass && directPass && objectivePass && labelPass;

fid = fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'Unified three-direct-method experiment\n');
fprintf(fid,'Runs per dataset: %d\n',runs);
fprintf(fid,'Methods: GNMFLD, GOCNMF, CGC-GOCNMF\n');
fprintf(fid,'Primary metrics: unlabeled samples only\n');
fprintf(fid,'Clustering rule: direct row-wise argmax\n');
fprintf(fid,'K-means used: NO\n\n');

for i = 1:length(summaryRows)
    S = summaryRows(i);
    fprintf(fid,['%s %s: ACC=%.6f+-%.6f, NMI=%.6f+-%.6f, ' ...
        'labeled_ACC=%.6f, theta=%.6f, time=%.3fs\n'], ...
        S.dataset,S.method,S.unlabeled_ACC_mean, ...
        S.unlabeled_ACC_sd,S.unlabeled_NMI_mean, ...
        S.unlabeled_NMI_sd,S.labeled_ACC_mean,S.theta_mean, ...
        S.run_seconds_mean);
end

fprintf(fid,'\nPaired comparisons\n');
for i = 1:length(pairedRows)
    S = pairedRows(i);
    fprintf(fid,['%s %s: dACC=%+.6f (%d/%d/%d), ' ...
        'dNMI=%+.6f (%d/%d/%d)\n'], ...
        S.dataset,S.comparison,S.dACC_mean, ...
        S.ACC_wins,S.ACC_ties,S.ACC_losses,S.dNMI_mean, ...
        S.NMI_wins,S.NMI_ties,S.NMI_losses);
end

fprintf(fid,'\nAudit\n');
fprintf(fid,'FINITE_RESULTS_PASS=%d\n',finitePass);
fprintf(fid,'DIRECT_ARGMAX_PASS=%d\n',directPass);
fprintf(fid,'OBJECTIVE_NONINCREASE_PASS=%d\n',objectivePass);
fprintf(fid,'GOC_FIXED_LABEL_PASS=%d\n',labelPass);
fprintf(fid,'GNMFLD_TRANSFER_PARAMETER_DATASETS=%d\n',transferCount);
fprintf(fid,'THREE_METHOD_INTEGRITY_PASS=%d\n',integrityPass);
if runs==1
    fprintf(fid,'ONE_SEED_SMOKE_PASS=%d\n',integrityPass);
elseif runs==3
    fprintf(fid,'THREE_SEED_INTEGRITY_PASS=%d\n',integrityPass);
else
    fprintf(fid,'TWENTY_SEED_COMPARISON_COMPLETE=%d\n',integrityPass);
end
fprintf(fid,['\nEditorial note: COIL100 and Optdigits use explicitly ' ...
    'marked same-domain GNMFLD parameter transfers. The final manuscript ' ...
    'must either retain this disclosure or replace them with directly ' ...
    'verified author-code/original-paper settings.\n']);
fclose(fid);
end
