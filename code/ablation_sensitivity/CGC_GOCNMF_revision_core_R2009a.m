function CGC_GOCNMF_revision_core_R2009a(mode,runs,dataDir,outputDir)
%CGC_GOCNMF_REVISION_CORE_R2009A
% Shared core for:
%   1. six-dataset BASE/LOCAL/HARD/SHRINK ablation;
%   2. shrinkage-rule sensitivity over kappa and F.
%
% Frozen protocol:
%   - 10% labeled samples per class, with at least two labels;
%   - evaluation on unlabeled samples only;
%   - direct row-wise argmax, no K-means;
%   - 50 GOCNMF iterations;
%   - identical labeled splits and initial factors within every paired run;
%   - epsilon = 1e-12;
%   - formal seeds begin at 20260617.
%
% MATLAB R2009a compatible.

if nargin < 1 || isempty(mode), error('mode is required.'); end
if nargin < 2 || isempty(runs), runs = 3; end
if nargin < 3 || isempty(dataDir), dataDir = pwd; end
if nargin < 4 || isempty(outputDir), outputDir = pwd; end

if ~(runs == 3 || runs == 20)
    error('runs must be 3 or 20.');
end
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

if strcmpi(mode,'ablation6')
    cgc_run_ablation6(runs,dataDir,outputDir);
elseif strcmpi(mode,'sensitivity3')
    cgc_run_sensitivity3(runs,dataDir,outputDir);
else
    error('Unknown mode: %s',mode);
end
end


function cgc_run_ablation6(runs,dataDir,outputDir)
EPSILON = 1e-12;
LABEL_FRACTION = 0.10;
ITERATIONS = 50;
FOLDS = 5;
SEED_START = 20260617;
BLOCK_SIZE = 256;
TOL = 1e-12;

D = cgc_dataset_spec();
numberOfDatasets = length(D);

prefix = sprintf('CGC_GOCNMF_ablation6_%dseed',runs);
rawPath = fullfile(outputDir,[prefix '_raw.csv']);
summaryPath = fullfile(outputDir,[prefix '_summary.csv']);
decisionPath = fullfile(outputDir,[prefix '_decision.txt']);
latexPath = fullfile(outputDir,[prefix '_table.tex']);
calibrationLatexPath = fullfile(outputDir,[prefix '_calibration_table.tex']);
checkpointPath = fullfile(outputDir,[prefix '_checkpoint.mat']);

template = struct( ...
    'dataset','', 'seed',0, 'p',0, 'alpha',0, 'labeled_count',0, ...
    'theta',0, 'delta_cv',0, 'se_cv',0, ...
    'binary_cv_brier',0, 'local_cv_brier',0, ...
    'hard_choice',0, ...
    'base_ACC',0, 'local_ACC',0, 'hard_ACC',0, 'shrink_ACC',0, ...
    'base_NMI',0, 'local_NMI',0, 'hard_NMI',0, 'shrink_NMI',0);

rows = repmat(template,1,numberOfDatasets*runs);
completed = false(numberOfDatasets,runs);

if exist(checkpointPath,'file') == 2
    P = load(checkpointPath);
    if isfield(P,'rows') && isfield(P,'completed') && ...
            length(P.rows) == numberOfDatasets*runs && ...
            all(size(P.completed) == [numberOfDatasets runs])
        rows = P.rows;
        completed = P.completed;
        fprintf('Resuming checkpoint: %s\n',checkpointPath);
    else
        error('Existing ablation checkpoint is incompatible.');
    end
end

fprintf('\nSIX-DATASET CGC-GOCNMF ABLATION\n');
fprintf('Runs per dataset: %d\n',runs);
fprintf('Variants: BASE, LOCAL, HARD, SHRINK\n\n');

for d = 1:numberOfDatasets
    dataPath = cgc_find_file_recursive(dataDir,D(d).files);
    fprintf('Loading %s: %s\n',D(d).name,dataPath);
    [X,y] = cgc_load_dataset( ...
        dataPath,D(d).expectedSamples,D(d).expectedClasses);
    c = length(unique(y));

    fprintf('Constructing graphs for %s.\n',D(d).name);
    [W0,W1] = cgc_construct_graphs(X,D(d).p,BLOCK_SIZE,EPSILON);

    for r = 1:runs
        if completed(d,r)
            fprintf('%s seed %d already complete; skipping.\n', ...
                D(d).name,SEED_START+r-1);
            continue;
        end

        seed = SEED_START+r-1;
        L = cgc_labeled_indices(y,seed,LABEL_FRACTION);

        [deltaCV,seCV,loss0,loss1] = cgc_calibration_statistics( ...
            W0,W1,L,y,c,seed,FOLDS,EPSILON);
        theta = cgc_theta(deltaCV,seCV,1,EPSILON);
        hardChoice = double(deltaCV > 0);

        Vbase = cgc_train(X,y,L,W0,D(d).alpha,seed,ITERATIONS,EPSILON);
        Vlocal = cgc_train(X,y,L,W1,D(d).alpha,seed,ITERATIONS,EPSILON);
        [baseACC,baseNMI] = cgc_evaluate(Vbase,y,L,EPSILON);
        [localACC,localNMI] = cgc_evaluate(Vlocal,y,L,EPSILON);

        if hardChoice == 1
            hardACC = localACC;
            hardNMI = localNMI;
        else
            hardACC = baseACC;
            hardNMI = baseNMI;
        end

        if theta <= TOL
            shrinkACC = baseACC;
            shrinkNMI = baseNMI;
        elseif theta >= 1-TOL
            shrinkACC = localACC;
            shrinkNMI = localNMI;
        else
            Wtheta = (1-theta)*W0 + theta*W1;
            Vshrink = cgc_train( ...
                X,y,L,Wtheta,D(d).alpha,seed,ITERATIONS,EPSILON);
            [shrinkACC,shrinkNMI] = cgc_evaluate( ...
                Vshrink,y,L,EPSILON);
        end

        index = (d-1)*runs+r;
        S = template;
        S.dataset = D(d).name;
        S.seed = seed;
        S.p = D(d).p;
        S.alpha = D(d).alpha;
        S.labeled_count = length(L);
        S.theta = theta;
        S.delta_cv = deltaCV;
        S.se_cv = seCV;
        S.binary_cv_brier = loss0;
        S.local_cv_brier = loss1;
        S.hard_choice = hardChoice;
        S.base_ACC = baseACC;
        S.local_ACC = localACC;
        S.hard_ACC = hardACC;
        S.shrink_ACC = shrinkACC;
        S.base_NMI = baseNMI;
        S.local_NMI = localNMI;
        S.hard_NMI = hardNMI;
        S.shrink_NMI = shrinkNMI;
        rows(index) = S;
        completed(d,r) = true;

        save(checkpointPath,'rows','completed','D','runs');
        cgc_write_ablation_raw(rawPath,rows,completed,runs);

        fprintf(['%s seed=%d theta=%.4f dACC=' ...
            '%+.5f/%+.5f/%+.5f dNMI=%+.5f/%+.5f/%+.5f\n'], ...
            D(d).name,seed,theta, ...
            localACC-baseACC,hardACC-baseACC,shrinkACC-baseACC, ...
            localNMI-baseNMI,hardNMI-baseNMI,shrinkNMI-baseNMI);
    end
end

if ~all(completed(:))
    error('Ablation experiment is incomplete.');
end

summary = cgc_summarize_ablation(rows,D,runs,TOL);
cgc_write_ablation_raw(rawPath,rows,completed,runs);
cgc_write_ablation_summary(summaryPath,summary);
cgc_write_ablation_latex(latexPath,summary);
cgc_write_calibration_latex(calibrationLatexPath,summary);
cgc_write_ablation_decision(decisionPath,runs,rows,summary,TOL);
save(checkpointPath,'rows','completed','summary','D','runs');

fprintf('\nAblation complete.\n');
fprintf('RAW=%s\n',rawPath);
fprintf('SUMMARY=%s\n',summaryPath);
fprintf('TABLE=%s\n',latexPath);
fprintf('CALIBRATION_TABLE=%s\n',calibrationLatexPath);
fprintf('DECISION=%s\n',decisionPath);
end


function cgc_run_sensitivity3(runs,dataDir,outputDir)
EPSILON = 1e-12;
LABEL_FRACTION = 0.10;
ITERATIONS = 50;
SEED_START = 20260617;
BLOCK_SIZE = 256;
TOL = 1e-12;
FGRID = [3 5 10];
KGRID = [0 0.5 1 1.5];

allD = cgc_dataset_spec();
% Pre-specified regimes: fallback, borderline, and clear benefit.
D = allD([1 2 4]); % PIE, YaleB, COIL100

numberOfDatasets = length(D);
numberOfF = length(FGRID);
numberOfK = length(KGRID);
rowsPerRun = numberOfF*numberOfK;

prefix = sprintf('CGC_GOCNMF_sensitivity3_%dseed',runs);
rawPath = fullfile(outputDir,[prefix '_raw.csv']);
summaryPath = fullfile(outputDir,[prefix '_summary.csv']);
decisionPath = fullfile(outputDir,[prefix '_decision.txt']);
kappaLatexPath = fullfile(outputDir,[prefix '_kappa_table.tex']);
foldLatexPath = fullfile(outputDir,[prefix '_fold_table.tex']);
checkpointPath = fullfile(outputDir,[prefix '_checkpoint.mat']);

template = struct( ...
    'dataset','', 'seed',0, 'F',0, 'kappa',0, ...
    'p',0, 'alpha',0, 'labeled_count',0, ...
    'delta_cv',0, 'se_cv',0, 'theta',0, ...
    'binary_cv_brier',0, 'local_cv_brier',0, ...
    'base_ACC',0, 'model_ACC',0, ...
    'base_NMI',0, 'model_NMI',0);

totalRows = numberOfDatasets*runs*rowsPerRun;
rows = repmat(template,1,totalRows);
completed = false(numberOfDatasets,runs);

if exist(checkpointPath,'file') == 2
    P = load(checkpointPath);
    if isfield(P,'rows') && isfield(P,'completed') && ...
            length(P.rows) == totalRows && ...
            all(size(P.completed) == [numberOfDatasets runs])
        rows = P.rows;
        completed = P.completed;
        fprintf('Resuming checkpoint: %s\n',checkpointPath);
    else
        error('Existing sensitivity checkpoint is incompatible.');
    end
end

fprintf('\nCGC-GOCNMF SHRINKAGE SENSITIVITY\n');
fprintf('Runs per dataset: %d\n',runs);
fprintf('Datasets: PIE, YaleB, COIL100\n');
fprintf('F values: 3, 5, 10\n');
fprintf('kappa values: 0, 0.5, 1, 1.5\n\n');

for d = 1:numberOfDatasets
    dataPath = cgc_find_file_recursive(dataDir,D(d).files);
    fprintf('Loading %s: %s\n',D(d).name,dataPath);
    [X,y] = cgc_load_dataset( ...
        dataPath,D(d).expectedSamples,D(d).expectedClasses);
    c = length(unique(y));

    fprintf('Constructing graphs for %s.\n',D(d).name);
    [W0,W1] = cgc_construct_graphs(X,D(d).p,BLOCK_SIZE,EPSILON);

    for r = 1:runs
        if completed(d,r)
            fprintf('%s seed %d already complete; skipping.\n', ...
                D(d).name,SEED_START+r-1);
            continue;
        end

        seed = SEED_START+r-1;
        L = cgc_labeled_indices(y,seed,LABEL_FRACTION);

        Vbase = cgc_train(X,y,L,W0,D(d).alpha,seed,ITERATIONS,EPSILON);
        Vlocal = cgc_train(X,y,L,W1,D(d).alpha,seed,ITERATIONS,EPSILON);
        [baseACC,baseNMI] = cgc_evaluate(Vbase,y,L,EPSILON);
        [localACC,localNMI] = cgc_evaluate(Vlocal,y,L,EPSILON);

        localCounter = 0;
        for fi = 1:numberOfF
            F = FGRID(fi);
            [deltaCV,seCV,loss0,loss1] = cgc_calibration_statistics( ...
                W0,W1,L,y,c,seed,F,EPSILON);

            for ki = 1:numberOfK
                kappa = KGRID(ki);
                theta = cgc_theta(deltaCV,seCV,kappa,EPSILON);

                if theta <= TOL
                    modelACC = baseACC;
                    modelNMI = baseNMI;
                elseif theta >= 1-TOL
                    modelACC = localACC;
                    modelNMI = localNMI;
                else
                    Wtheta = (1-theta)*W0 + theta*W1;
                    Vmodel = cgc_train( ...
                        X,y,L,Wtheta,D(d).alpha,seed,ITERATIONS,EPSILON);
                    [modelACC,modelNMI] = cgc_evaluate( ...
                        Vmodel,y,L,EPSILON);
                end

                localCounter = localCounter+1;
                index = ((d-1)*runs+(r-1))*rowsPerRun+localCounter;
                S = template;
                S.dataset = D(d).name;
                S.seed = seed;
                S.F = F;
                S.kappa = kappa;
                S.p = D(d).p;
                S.alpha = D(d).alpha;
                S.labeled_count = length(L);
                S.delta_cv = deltaCV;
                S.se_cv = seCV;
                S.theta = theta;
                S.binary_cv_brier = loss0;
                S.local_cv_brier = loss1;
                S.base_ACC = baseACC;
                S.model_ACC = modelACC;
                S.base_NMI = baseNMI;
                S.model_NMI = modelNMI;
                rows(index) = S;

                fprintf(['%s seed=%d F=%d kappa=%.1f theta=%.4f ' ...
                    'dACC=%+.5f dNMI=%+.5f\n'], ...
                    D(d).name,seed,F,kappa,theta, ...
                    modelACC-baseACC,modelNMI-baseNMI);
            end
        end

        completed(d,r) = true;
        save(checkpointPath,'rows','completed','D','runs','FGRID','KGRID');
        cgc_write_sensitivity_raw(rawPath,rows,completed, ...
            runs,rowsPerRun);
    end
end

if ~all(completed(:))
    error('Sensitivity experiment is incomplete.');
end

summary = cgc_summarize_sensitivity( ...
    rows,D,runs,FGRID,KGRID,TOL);
cgc_write_sensitivity_raw(rawPath,rows,completed,runs,rowsPerRun);
cgc_write_sensitivity_summary(summaryPath,summary);
cgc_write_kappa_latex(kappaLatexPath,summary,KGRID);
cgc_write_fold_latex(foldLatexPath,summary,FGRID);
cgc_write_sensitivity_decision( ...
    decisionPath,runs,rows,summary,FGRID,KGRID,TOL);
save(checkpointPath,'rows','completed','summary','D','runs','FGRID','KGRID');

fprintf('\nSensitivity experiment complete.\n');
fprintf('RAW=%s\n',rawPath);
fprintf('SUMMARY=%s\n',summaryPath);
fprintf('KAPPA_TABLE=%s\n',kappaLatexPath);
fprintf('FOLD_TABLE=%s\n',foldLatexPath);
fprintf('DECISION=%s\n',decisionPath);
end


function D = cgc_dataset_spec()
D(1).name = 'PIE';
D(1).files = {'CMU_PIE_fac.mat','CMU_PIE.mat','PIE.mat','PIE_fac.mat'};
D(1).p = 3; D(1).alpha = 1000;
D(1).expectedSamples = 2856; D(1).expectedClasses = 68;

D(2).name = 'YaleB';
D(2).files = {'YaleB.mat','YaleB(1).mat','YaleB_32x32.mat'};
D(2).p = 2; D(2).alpha = 1000;
D(2).expectedSamples = 2414; D(2).expectedClasses = 38;

D(3).name = 'COIL20';
D(3).files = {'COIL20_Obj.mat','COIL20.mat','COIL20_Obj(1).mat'};
D(3).p = 3; D(3).alpha = 10;
D(3).expectedSamples = 1440; D(3).expectedClasses = 20;

D(4).name = 'COIL100';
D(4).files = {'COIL100_Obj.mat','COIL100_Obj(1).mat','COIL100.mat'};
D(4).p = 3; D(4).alpha = 10;
D(4).expectedSamples = 7200; D(4).expectedClasses = 100;

D(5).name = 'Optdigits';
D(5).files = {'Optdigits_Han.mat','Optdigits.mat','optdigits.mat'};
D(5).p = 4; D(5).alpha = 10;
D(5).expectedSamples = 5620; D(5).expectedClasses = 10;

D(6).name = 'MNIST';
D(6).files = {'MNIST_Han.mat','MNIST.mat','mnist.mat'};
D(6).p = 4; D(6).alpha = 10;
D(6).expectedSamples = 6996; D(6).expectedClasses = 10;
end


function path = cgc_find_file_recursive(folder,candidates)
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
    if entries(i).isdir && ...
            ~strcmp(entries(i).name,'.') && ...
            ~strcmp(entries(i).name,'..')
        child = fullfile(folder,entries(i).name);
        try
            path = cgc_find_file_recursive(child,candidates);
        catch
            path = '';
        end
        if ~isempty(path)
            return;
        end
    end
end

error('Required dataset was not found under %s.',folder);
end


function [X,y] = cgc_load_dataset(path,expectedSamples,expectedClasses)
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
    error('%s: feature array must be two-dimensional.',path);
end
if size(F,1) == length(y0)
    X = F';
elseif size(F,2) == length(y0)
    X = F;
else
    error('%s: feature shape does not match labels.',path);
end
if min(X(:)) < -1e-12
    error('%s: negative feature found.',path);
end

[classValues,dummy,y] = unique(y0); %#ok<ASGLU>
y = double(y(:));

if size(X,2) ~= expectedSamples || ...
        length(classValues) ~= expectedClasses
    error('%s: expected %d samples/%d classes, obtained %d/%d.', ...
        path,expectedSamples,expectedClasses, ...
        size(X,2),length(classValues));
end

norms = sqrt(sum(X.^2,1));
norms(norms == 0) = 1;
X = bsxfun(@rdivide,X,norms);

if any(~isfinite(X(:)))
    error('%s: nonfinite normalized feature found.',path);
end
end


function cgc_set_seed(seed)
rand('twister',double(seed)); %#ok<RAND>
end


function L = cgc_labeled_indices(y,seed,fraction)
cgc_set_seed(seed);
classes = unique(y);
L = zeros(0,1);

for k = 1:length(classes)
    ids = find(y == classes(k));
    order = randperm(length(ids));
    count = max(2,floor(fraction*length(ids)));
    count = min(count,length(ids));
    L = [L;ids(order(1:count))]; %#ok<AGROW>
end

L = double(L(:));
end


function [W0,W1] = cgc_construct_graphs(X,p,blockSize,epsilon)
n = size(X,2);
rows = zeros(n*p,1);
cols = zeros(n*p,1);
dvals = zeros(n*p,1);
sigma = zeros(n,1);
position = 0;

for first = 1:blockSize:n
    last = min(first+blockSize-1,n);
    ids = first:last;
    similarities = X(:,ids)'*X;

    for q = 1:length(ids)
        similarities(q,ids(q)) = -Inf;
    end

    [sortedSimilarity,sortedIndex] = sort( ...
        similarities,2,'descend');
    neighborIndex = sortedIndex(:,1:p);
    neighborDistance = max(0,1-sortedSimilarity(:,1:p));

    rowBlock = repmat(ids(:),1,p);
    count = length(ids)*p;
    range = position+(1:count);
    rows(range) = rowBlock(:);
    cols(range) = neighborIndex(:);
    dvals(range) = neighborDistance(:);
    sigma(ids) = max(neighborDistance(:,p),epsilon);
    position = position+count;

    fprintf('  kNN block %d:%d of %d\n',first,last,n);
end

A = sparse(rows,cols,ones(length(rows),1),n,n);
W0 = spones(A+A');
W0 = W0-spdiags(diag(W0),0,n,n);
W0 = sparse(W0);

denominator = sqrt(sigma(rows).*sigma(cols))+epsilon;
weights = exp(-(dvals./denominator).^2);
A1 = sparse(rows,cols,weights,n,n);
W1 = max(A1,A1');
W1 = W1-spdiags(diag(W1),0,n,n);

values = nonzeros(W1);
if isempty(values)
    error('The locally weighted graph is empty.');
end
W1 = W1/mean(values);
W1 = sparse(W1);

if nnz(spones(W0)-spones(W1)) ~= 0
    error('Candidate graph supports do not match.');
end
end


function folds = cgc_stratified_folds(L,y,seed,numberOfFolds)
cgc_set_seed(seed);
folds = cell(numberOfFolds,1);
for f = 1:numberOfFolds
    folds{f} = zeros(0,1);
end

classes = unique(y(L));
for k = 1:length(classes)
    ids = L(y(L) == classes(k));
    order = randperm(length(ids));
    ids = ids(order);
    offset = mod(k-1,numberOfFolds);

    for j = 1:length(ids)
        f = mod((j-1)+offset,numberOfFolds)+1;
        folds{f} = [folds{f};ids(j)]; %#ok<AGROW>
    end
end

for f = 1:numberOfFolds
    if isempty(folds{f})
        error('Cross-validation fold %d is empty.',f);
    end
end
end


function [delta,se,loss0Mean,loss1Mean] = ...
    cgc_calibration_statistics(W0,W1,L,y,c,seed,F,epsilon)

folds = cgc_stratified_folds(L,y,seed+991,F);
loss0 = cgc_fold_brier_losses(W0,folds,L,y,c,epsilon);
loss1 = cgc_fold_brier_losses(W1,folds,L,y,c,epsilon);

if any(~isfinite(loss0)) || any(~isfinite(loss1))
    error('Nonfinite cross-validation Brier loss detected.');
end

difference = loss0-loss1;
delta = mean(difference);

if length(difference) > 1
    se = std(difference,0)/sqrt(length(difference));
else
    se = 0;
end

loss0Mean = mean(loss0);
loss1Mean = mean(loss1);
end


function theta = cgc_theta(delta,se,kappa,epsilon)
if delta > 0
    theta = max(0,(delta-kappa*se)/delta);
else
    theta = 0;
end
theta = min(1,max(0,theta));
end


function losses = cgc_fold_brier_losses(W,folds,L,y,c,epsilon)
numberOfFolds = length(folds);
losses = zeros(numberOfFolds,1);

for f = 1:numberOfFolds
    heldOut = folds{f};
    train = setdiff(L,heldOut);
    prediction = cgc_harmonic_predictions( ...
        W,train,y,c,heldOut,epsilon);

    target = zeros(length(heldOut),c);
    index = sub2ind(size(target),(1:length(heldOut))',y(heldOut));
    target(index) = 1;
    losses(f) = mean(sum((prediction-target).^2,2));
end
end


function P = cgc_harmonic_predictions(W,train,y,c,query,epsilon)
n = size(W,1);
unknownMask = true(n,1);
unknownMask(train) = false;
unknown = find(unknownMask);

position = zeros(n,1);
position(unknown) = 1:length(unknown);

degree = full(sum(W,2));
A = spdiags(degree(unknown),0,length(unknown),length(unknown)) ...
    -W(unknown,unknown);

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
[warningMessage,dummyWarningId] = lastwarn; %#ok<ASGLU>

if any(~isfinite(Fu(:))) || ...
        ~isempty(strfind(lower(warningMessage),'singular'))
    Fu = (A+1e-10*speye(length(unknown)))\rhs;
end

Fu = max(Fu,0);
rowSum = sum(Fu,2);
positive = rowSum > epsilon;
if any(positive)
    Fu(positive,:) = bsxfun(@rdivide, ...
        Fu(positive,:),rowSum(positive));
end

P = Fu(position(query),:);
end


function V = cgc_train(X,y,L,W,alpha,seed,iterations,epsilon)
[m,n] = size(X);
c = length(unique(y));

cgc_set_seed(seed);
U = max(rand(m,c),epsilon);
V = max(rand(n,c),epsilon);

C = zeros(length(L),c);
index = sub2ind(size(C),(1:length(L))',y(L));
C(index) = 1;
V(L,:) = C;

degree = full(sum(W,2));

for iter = 1:iterations
    numeratorU = X*V;
    denominatorU = U*(V'*V);
    U = U.*(numeratorU./max(denominatorU,epsilon));
    U = max(U,epsilon);

    numeratorV = X'*U+alpha*(W*V);
    denominatorV = V*(U'*U)+ ...
        alpha*bsxfun(@times,degree,V);
    V = V.*(numeratorV./max(denominatorV,epsilon));
    V = max(V,epsilon);
    V(L,:) = C;
end
end


function [acc,nmi] = cgc_evaluate(V,y,L,epsilon)
mask = true(length(y),1);
mask(L) = false;
U = find(mask);

[dummy,prediction] = max(V,[],2); %#ok<ASGLU>
acc = mean(prediction(U) == y(U));
nmi = cgc_nmi(y(U),prediction(U),epsilon);
end


function value = cgc_nmi(trueLabel,predictedLabel,epsilon)
trueLabel = trueLabel(:);
predictedLabel = predictedLabel(:);

[trueValues,dummy1,trueIndex] = unique(trueLabel); %#ok<ASGLU>
[predValues,dummy2,predIndex] = unique(predictedLabel); %#ok<ASGLU>

contingency = accumarray( ...
    [trueIndex predIndex],1,[length(trueValues) length(predValues)]);
P = contingency/length(trueLabel);
pi = sum(P,2);
pj = sum(P,1);
mi = 0;

for i = 1:size(P,1)
    for j = 1:size(P,2)
        if P(i,j) > 0
            mi = mi+P(i,j)*log(P(i,j)/max(pi(i)*pj(j),epsilon));
        end
    end
end

ht = -sum(pi(pi>0).*log(pi(pi>0)));
hp = -sum(pj(pj>0).*log(pj(pj>0)));
denominator = 0.5*(ht+hp);

if denominator <= epsilon
    value = 1;
else
    value = mi/denominator;
end
value = max(0,min(1,value));
end


function cgc_write_ablation_raw(path,rows,completed,runs)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,['dataset,seed,p,alpha,labeled_count,theta,delta_cv,se_cv,' ...
    'binary_cv_brier,local_cv_brier,hard_choice,' ...
    'base_ACC,local_ACC,hard_ACC,shrink_ACC,' ...
    'local_dACC,hard_dACC,shrink_dACC,' ...
    'base_NMI,local_NMI,hard_NMI,shrink_NMI,' ...
    'local_dNMI,hard_dNMI,shrink_dNMI\n']);

for d = 1:size(completed,1)
    for r = 1:runs
        if completed(d,r)
            S = rows((d-1)*runs+r);
            fprintf(fid,['%s,%d,%d,%.15g,%d,%.15g,%.15g,%.15g,' ...
                '%.15g,%.15g,%d,' ...
                '%.15g,%.15g,%.15g,%.15g,' ...
                '%.15g,%.15g,%.15g,' ...
                '%.15g,%.15g,%.15g,%.15g,' ...
                '%.15g,%.15g,%.15g\n'], ...
                S.dataset,S.seed,S.p,S.alpha,S.labeled_count, ...
                S.theta,S.delta_cv,S.se_cv, ...
                S.binary_cv_brier,S.local_cv_brier,S.hard_choice, ...
                S.base_ACC,S.local_ACC,S.hard_ACC,S.shrink_ACC, ...
                S.local_ACC-S.base_ACC, ...
                S.hard_ACC-S.base_ACC, ...
                S.shrink_ACC-S.base_ACC, ...
                S.base_NMI,S.local_NMI,S.hard_NMI,S.shrink_NMI, ...
                S.local_NMI-S.base_NMI, ...
                S.hard_NMI-S.base_NMI, ...
                S.shrink_NMI-S.base_NMI);
        end
    end
end

fclose(fid);
end


function summary = cgc_summarize_ablation(rows,D,runs,tol)
template = struct( ...
    'dataset','', 'runs',runs, ...
    'theta_mean',0, 'theta_sd',0, 'fallbacks',0, ...
    'delta_mean',0, 'se_mean',0, 'hard_local_count',0, ...
    'base_ACC_mean',0, 'base_NMI_mean',0, ...
    'local_dACC_mean',0, 'hard_dACC_mean',0, ...
    'shrink_dACC_mean',0, ...
    'local_dNMI_mean',0, 'hard_dNMI_mean',0, ...
    'shrink_dNMI_mean',0, ...
    'local_ACC_losses',0, 'hard_ACC_losses',0, ...
    'shrink_ACC_losses',0, ...
    'local_NMI_losses',0, 'hard_NMI_losses',0, ...
    'shrink_NMI_losses',0);

summary = repmat(template,1,length(D));

for d = 1:length(D)
    ids = (d-1)*runs+(1:runs);
    R = rows(ids);

    theta = [R.theta];
    delta = [R.delta_cv];
    se = [R.se_cv];

    baseACC = [R.base_ACC];
    localACC = [R.local_ACC];
    hardACC = [R.hard_ACC];
    shrinkACC = [R.shrink_ACC];

    baseNMI = [R.base_NMI];
    localNMI = [R.local_NMI];
    hardNMI = [R.hard_NMI];
    shrinkNMI = [R.shrink_NMI];

    S = template;
    S.dataset = D(d).name;
    S.theta_mean = mean(theta);
    S.theta_sd = cgc_std(theta);
    S.fallbacks = sum(theta <= tol);
    S.delta_mean = mean(delta);
    S.se_mean = mean(se);
    S.hard_local_count = sum([R.hard_choice] == 1);
    S.base_ACC_mean = mean(baseACC);
    S.base_NMI_mean = mean(baseNMI);

    S.local_dACC_mean = mean(localACC-baseACC);
    S.hard_dACC_mean = mean(hardACC-baseACC);
    S.shrink_dACC_mean = mean(shrinkACC-baseACC);

    S.local_dNMI_mean = mean(localNMI-baseNMI);
    S.hard_dNMI_mean = mean(hardNMI-baseNMI);
    S.shrink_dNMI_mean = mean(shrinkNMI-baseNMI);

    S.local_ACC_losses = sum(localACC-baseACC < -tol);
    S.hard_ACC_losses = sum(hardACC-baseACC < -tol);
    S.shrink_ACC_losses = sum(shrinkACC-baseACC < -tol);

    S.local_NMI_losses = sum(localNMI-baseNMI < -tol);
    S.hard_NMI_losses = sum(hardNMI-baseNMI < -tol);
    S.shrink_NMI_losses = sum(shrinkNMI-baseNMI < -tol);

    summary(d) = S;
end
end


function cgc_write_ablation_summary(path,summary)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,['dataset,runs,theta_mean,theta_sd,fallbacks,' ...
    'delta_mean,se_mean,hard_local_count,base_ACC_mean,base_NMI_mean,' ...
    'local_dACC_mean,hard_dACC_mean,shrink_dACC_mean,' ...
    'local_dNMI_mean,hard_dNMI_mean,shrink_dNMI_mean,' ...
    'local_ACC_losses,hard_ACC_losses,shrink_ACC_losses,' ...
    'local_NMI_losses,hard_NMI_losses,shrink_NMI_losses\n']);

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,['%s,%d,%.15g,%.15g,%d,%.15g,%.15g,%d,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
        '%d,%d,%d,%d,%d,%d\n'], ...
        S.dataset,S.runs,S.theta_mean,S.theta_sd,S.fallbacks, ...
        S.delta_mean,S.se_mean,S.hard_local_count, ...
        S.base_ACC_mean,S.base_NMI_mean, ...
        S.local_dACC_mean,S.hard_dACC_mean,S.shrink_dACC_mean, ...
        S.local_dNMI_mean,S.hard_dNMI_mean,S.shrink_dNMI_mean, ...
        S.local_ACC_losses,S.hard_ACC_losses,S.shrink_ACC_losses, ...
        S.local_NMI_losses,S.hard_NMI_losses,S.shrink_NMI_losses);
end

fclose(fid);
end


function cgc_write_ablation_latex(path,summary)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'\\begin{table*}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,['\\caption{Six-dataset ablation of local weighting and ' ...
    'conservative calibration.}\\label{tab:ablation6}\n']);
fprintf(fid,'\\small\n');
fprintf(fid,'\\begin{tabular}{llrrc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,['Dataset & Variant & $\\Delta$ACC & $\\Delta$NMI ' ...
    '& ACC/NMI losses \\\\\n']);
fprintf(fid,'\\midrule\n');

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,'\\multirow{3}{*}{%s} ',S.dataset);
    fprintf(fid,'& LOCAL & %+.3f & %+.3f & %d/%d \\\\\n', ...
        100*S.local_dACC_mean,100*S.local_dNMI_mean, ...
        S.local_ACC_losses,S.local_NMI_losses);
    fprintf(fid,'& HARD & %+.3f & %+.3f & %d/%d \\\\\n', ...
        100*S.hard_dACC_mean,100*S.hard_dNMI_mean, ...
        S.hard_ACC_losses,S.hard_NMI_losses);
    fprintf(fid,'& SHRINK & %+.3f & %+.3f & %d/%d \\\\\n', ...
        100*S.shrink_dACC_mean,100*S.shrink_dNMI_mean, ...
        S.shrink_ACC_losses,S.shrink_NMI_losses);
    if i < length(summary)
        fprintf(fid,'\\midrule\n');
    end
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n\n');
fprintf(fid,'\\vspace{1mm}\n');
fprintf(fid,'\\parbox{0.98\\textwidth}{\\footnotesize\n');
fprintf(fid,['Changes are absolute percentage points relative to BASE ' ...
    sprintf('over %d paired trials on unlabeled samples. ',summary(1).runs) ...
    'Losses are counted relative to BASE.}\n']);
fprintf(fid,'\\end{table*}\n');

fclose(fid);
end


function cgc_write_calibration_latex(path,summary)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'\\begin{table}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,['\\caption{Calibration behavior in the six-dataset ' ...
    'ablation.}\\label{tab:ablation6_calibration}\n']);
fprintf(fid,'\\small\n');
fprintf(fid,'\\begin{tabular}{lrrrrr}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,['Dataset & $\\bar\\Delta$ & $\\bar s$ & $\\bar\\theta$ ' ...
    '& Fallbacks & HARD local \\\\\n']);
fprintf(fid,'\\midrule\n');

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,'%s & %.5f & %.5f & %.3f & %d/%d & %d/%d \\\\\n', ...
        S.dataset,S.delta_mean,S.se_mean,S.theta_mean, ...
        S.fallbacks,S.runs,S.hard_local_count,S.runs);
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);
end


function cgc_write_ablation_decision(path,runs,rows,summary,tol)
numeric = [];
for i = 1:length(rows)
    numeric = [numeric; ...
        rows(i).theta; rows(i).delta_cv; rows(i).se_cv; ...
        rows(i).base_ACC; rows(i).local_ACC; ...
        rows(i).hard_ACC; rows(i).shrink_ACC; ...
        rows(i).base_NMI; rows(i).local_NMI; ...
        rows(i).hard_NMI; rows(i).shrink_NMI]; %#ok<AGROW>
end

finitePass = all(isfinite(numeric));
thetaPass = all([rows.theta] >= -tol & [rows.theta] <= 1+tol);
shrinkNoLoss = 1;

for i = 1:length(summary)
    if summary(i).shrink_ACC_losses > 0 || ...
            summary(i).shrink_NMI_losses > 0
        shrinkNoLoss = 0;
    end
end

fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'Six-dataset CGC-GOCNMF ablation audit\n');
fprintf(fid,'Runs per dataset: %d\n',runs);
fprintf(fid,'Variants: BASE, LOCAL, HARD, SHRINK\n\n');

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,['%s: theta=%.6f, fallbacks=%d/%d, ' ...
        'dACC local/hard/shrink=%+.6f/%+.6f/%+.6f, ' ...
        'ACC losses=%d/%d/%d, ' ...
        'dNMI local/hard/shrink=%+.6f/%+.6f/%+.6f, ' ...
        'NMI losses=%d/%d/%d\n'], ...
        S.dataset,S.theta_mean,S.fallbacks,S.runs, ...
        S.local_dACC_mean,S.hard_dACC_mean,S.shrink_dACC_mean, ...
        S.local_ACC_losses,S.hard_ACC_losses,S.shrink_ACC_losses, ...
        S.local_dNMI_mean,S.hard_dNMI_mean,S.shrink_dNMI_mean, ...
        S.local_NMI_losses,S.hard_NMI_losses,S.shrink_NMI_losses);
end

fprintf(fid,'\nFINITE_OUTPUT_PASS=%d\n',finitePass);
fprintf(fid,'THETA_RANGE_PASS=%d\n',thetaPass);
fprintf(fid,'SHRINK_NO_LOSS_OBSERVED=%d\n',shrinkNoLoss);

if runs == 3
    fprintf(fid,'THREE_SEED_ABLATION_INTEGRITY_PASS=%d\n', ...
        finitePass && thetaPass);
else
    fprintf(fid,'TWENTY_SEED_ABLATION_COMPLETE=%d\n', ...
        finitePass && thetaPass);
end

fclose(fid);
end


function cgc_write_sensitivity_raw(path,rows,completed,runs,rowsPerRun)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,['dataset,seed,F,kappa,p,alpha,labeled_count,' ...
    'delta_cv,se_cv,theta,binary_cv_brier,local_cv_brier,' ...
    'base_ACC,model_ACC,dACC,base_NMI,model_NMI,dNMI\n']);

for d = 1:size(completed,1)
    for r = 1:runs
        if completed(d,r)
            first = ((d-1)*runs+(r-1))*rowsPerRun+1;
            last = first+rowsPerRun-1;
            for index = first:last
                S = rows(index);
                fprintf(fid,['%s,%d,%d,%.15g,%d,%.15g,%d,' ...
                    '%.15g,%.15g,%.15g,%.15g,%.15g,' ...
                    '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n'], ...
                    S.dataset,S.seed,S.F,S.kappa,S.p,S.alpha, ...
                    S.labeled_count,S.delta_cv,S.se_cv,S.theta, ...
                    S.binary_cv_brier,S.local_cv_brier, ...
                    S.base_ACC,S.model_ACC,S.model_ACC-S.base_ACC, ...
                    S.base_NMI,S.model_NMI,S.model_NMI-S.base_NMI);
            end
        end
    end
end

fclose(fid);
end


function summary = cgc_summarize_sensitivity( ...
    rows,D,runs,FGRID,KGRID,tol)

template = struct( ...
    'dataset','', 'F',0, 'kappa',0, 'runs',runs, ...
    'theta_mean',0, 'theta_sd',0, 'fallbacks',0, ...
    'delta_mean',0, 'se_mean',0, ...
    'base_ACC_mean',0, 'model_ACC_mean',0, ...
    'dACC_mean',0, 'ACC_losses',0, ...
    'base_NMI_mean',0, 'model_NMI_mean',0, ...
    'dNMI_mean',0, 'NMI_losses',0);

summary = repmat(template,1,length(D)*length(FGRID)*length(KGRID));
counter = 0;

for d = 1:length(D)
    for fi = 1:length(FGRID)
        for ki = 1:length(KGRID)
            selected = false(1,length(rows));
            for j = 1:length(rows)
                selected(j) = strcmp(rows(j).dataset,D(d).name) && ...
                    rows(j).F == FGRID(fi) && ...
                    abs(rows(j).kappa-KGRID(ki)) <= 1e-14;
            end
            R = rows(selected);

            if length(R) ~= runs
                error('Sensitivity summary selection failed.');
            end

            counter = counter+1;
            S = template;
            S.dataset = D(d).name;
            S.F = FGRID(fi);
            S.kappa = KGRID(ki);

            theta = [R.theta];
            baseACC = [R.base_ACC];
            modelACC = [R.model_ACC];
            baseNMI = [R.base_NMI];
            modelNMI = [R.model_NMI];

            S.theta_mean = mean(theta);
            S.theta_sd = cgc_std(theta);
            S.fallbacks = sum(theta <= tol);
            S.delta_mean = mean([R.delta_cv]);
            S.se_mean = mean([R.se_cv]);

            S.base_ACC_mean = mean(baseACC);
            S.model_ACC_mean = mean(modelACC);
            S.dACC_mean = mean(modelACC-baseACC);
            S.ACC_losses = sum(modelACC-baseACC < -tol);

            S.base_NMI_mean = mean(baseNMI);
            S.model_NMI_mean = mean(modelNMI);
            S.dNMI_mean = mean(modelNMI-baseNMI);
            S.NMI_losses = sum(modelNMI-baseNMI < -tol);

            summary(counter) = S;
        end
    end
end
end


function cgc_write_sensitivity_summary(path,summary)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,['dataset,F,kappa,runs,theta_mean,theta_sd,fallbacks,' ...
    'delta_mean,se_mean,base_ACC_mean,model_ACC_mean,dACC_mean,' ...
    'ACC_losses,base_NMI_mean,model_NMI_mean,dNMI_mean,NMI_losses\n']);

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,['%s,%d,%.15g,%d,%.15g,%.15g,%d,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%d,%.15g,%.15g,%.15g,%d\n'], ...
        S.dataset,S.F,S.kappa,S.runs,S.theta_mean,S.theta_sd, ...
        S.fallbacks,S.delta_mean,S.se_mean, ...
        S.base_ACC_mean,S.model_ACC_mean,S.dACC_mean,S.ACC_losses, ...
        S.base_NMI_mean,S.model_NMI_mean,S.dNMI_mean,S.NMI_losses);
end

fclose(fid);
end


function cgc_write_kappa_latex(path,summary,KGRID)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'\\begin{table*}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,['\\caption{Sensitivity to the shrinkage strength ' ...
    '$\\kappa$ with $F=5$.}\\label{tab:kappa_sensitivity}\n']);
fprintf(fid,'\\small\n');
fprintf(fid,'\\begin{tabular}{lrrrrrr}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,['Dataset & $\\kappa$ & $\\bar\\theta$ & Fallbacks ' ...
    '& $\\Delta$ACC & $\\Delta$NMI & ACC/NMI losses \\\\\n']);
fprintf(fid,'\\midrule\n');

datasets = {'PIE','YaleB','COIL100'};
for d = 1:length(datasets)
    for k = 1:length(KGRID)
        S = cgc_find_summary(summary,datasets{d},5,KGRID(k));
        if k == 1
            fprintf(fid,'\\multirow{%d}{*}{%s} ',length(KGRID),datasets{d});
        else
            fprintf(fid,' ');
        end
        fprintf(fid,'& %.1f & %.3f & %d/%d & %+.3f & %+.3f & %d/%d \\\\\n', ...
            S.kappa,S.theta_mean,S.fallbacks,S.runs, ...
            100*S.dACC_mean,100*S.dNMI_mean, ...
            S.ACC_losses,S.NMI_losses);
    end
    if d < length(datasets), fprintf(fid,'\\midrule\n'); end
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{table*}\n');

fclose(fid);
end


function cgc_write_fold_latex(path,summary,FGRID)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'\\begin{table*}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,['\\caption{Sensitivity to the number of folds with ' ...
    '$\\kappa=1$.}\\label{tab:fold_sensitivity}\n']);
fprintf(fid,'\\small\n');
fprintf(fid,'\\begin{tabular}{lrrrrrr}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,['Dataset & $F$ & $\\bar\\theta$ & Fallbacks ' ...
    '& $\\Delta$ACC & $\\Delta$NMI & ACC/NMI losses \\\\\n']);
fprintf(fid,'\\midrule\n');

datasets = {'PIE','YaleB','COIL100'};
for d = 1:length(datasets)
    for f = 1:length(FGRID)
        S = cgc_find_summary(summary,datasets{d},FGRID(f),1);
        if f == 1
            fprintf(fid,'\\multirow{%d}{*}{%s} ',length(FGRID),datasets{d});
        else
            fprintf(fid,' ');
        end
        fprintf(fid,'& %d & %.3f & %d/%d & %+.3f & %+.3f & %d/%d \\\\\n', ...
            S.F,S.theta_mean,S.fallbacks,S.runs, ...
            100*S.dACC_mean,100*S.dNMI_mean, ...
            S.ACC_losses,S.NMI_losses);
    end
    if d < length(datasets), fprintf(fid,'\\midrule\n'); end
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{table*}\n');

fclose(fid);
end


function S = cgc_find_summary(summary,dataset,F,kappa)
found = 0;
S = summary(1);

for i = 1:length(summary)
    if strcmp(summary(i).dataset,dataset) && ...
            summary(i).F == F && ...
            abs(summary(i).kappa-kappa) <= 1e-14
        S = summary(i);
        found = 1;
        break;
    end
end

if ~found
    error('Requested sensitivity summary was not found.');
end
end


function cgc_write_sensitivity_decision( ...
    path,runs,rows,summary,FGRID,KGRID,tol)

numeric = [];
for i = 1:length(rows)
    numeric = [numeric; ...
        rows(i).delta_cv; rows(i).se_cv; rows(i).theta; ...
        rows(i).base_ACC; rows(i).model_ACC; ...
        rows(i).base_NMI; rows(i).model_NMI]; %#ok<AGROW>
end

finitePass = all(isfinite(numeric));
thetaPass = all([rows.theta] >= -tol & [rows.theta] <= 1+tol);

referenceCount = 0;
for i = 1:length(summary)
    if summary(i).F == 5 && abs(summary(i).kappa-1) <= 1e-14
        referenceCount = referenceCount+1;
    end
end
referencePass = (referenceCount == 3);

fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'CGC-GOCNMF shrinkage sensitivity audit\n');
fprintf(fid,'Runs per dataset: %d\n',runs);
fprintf(fid,'Datasets: PIE, YaleB, COIL100\n');
fprintf(fid,'F grid:');
for i = 1:length(FGRID), fprintf(fid,' %d',FGRID(i)); end
fprintf(fid,'\nkappa grid:');
for i = 1:length(KGRID), fprintf(fid,' %.1f',KGRID(i)); end
fprintf(fid,'\n\n');

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,['%s F=%d kappa=%.1f: theta=%.6f, fallback=%d/%d, ' ...
        'dACC=%+.6f, ACC losses=%d, dNMI=%+.6f, NMI losses=%d\n'], ...
        S.dataset,S.F,S.kappa,S.theta_mean,S.fallbacks,S.runs, ...
        S.dACC_mean,S.ACC_losses,S.dNMI_mean,S.NMI_losses);
end

fprintf(fid,'\nFINITE_OUTPUT_PASS=%d\n',finitePass);
fprintf(fid,'THETA_RANGE_PASS=%d\n',thetaPass);
fprintf(fid,'REFERENCE_F5_KAPPA1_PRESENT_PASS=%d\n',referencePass);

if runs == 3
    fprintf(fid,'THREE_SEED_SENSITIVITY_INTEGRITY_PASS=%d\n', ...
        finitePass && thetaPass && referencePass);
else
    fprintf(fid,'TWENTY_SEED_SENSITIVITY_COMPLETE=%d\n', ...
        finitePass && thetaPass && referencePass);
end

fclose(fid);
end


function value = cgc_std(x)
if length(x) > 1
    value = std(x,0);
else
    value = 0;
end
end
