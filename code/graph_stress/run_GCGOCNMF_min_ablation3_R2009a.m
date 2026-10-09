function run_GCGOCNMF_min_ablation3_R2009a(runs, dataDir, outputDir)
%RUN_GCGOCNMF_MIN_ABLATION3_R2009A
% MATLAB R2009a-compatible minimal three-dataset ablation for conservative
% graph-calibrated GOCNMF.
%
% Four graph variants are compared under identical labels and initialization:
%   BASE   : original binary p-nearest-neighbor graph W0;
%   LOCAL  : local self-tuning weighted graph W1;
%   HARD   : W1 if mean cross-fitted Brier improvement Delta>0, else W0;
%   SHRINK : (1-theta)W0+theta W1, theta=max(0,1-SE/Delta) for Delta>0.
%
% Usage:
%   run_GCGOCNMF_min_ablation3_R2009a(3,  pwd, pwd);
%   run_GCGOCNMF_min_ablation3_R2009a(20, pwd, pwd);
%
% The model, graph construction, label rate, folds, alpha, p, iteration count,
% and random seeds are frozen. No K-means or external classifier is used.

if nargin < 1 || isempty(runs)
    runs = 3;
end
if nargin < 2 || isempty(dataDir)
    dataDir = pwd;
end
if nargin < 3 || isempty(outputDir)
    outputDir = pwd;
end
if ~(runs == 3 || runs == 20)
    error('runs must be 3 or 20.');
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

EPSILON = 1e-12;
LABEL_FRACTION = 0.10;
ITERATIONS = 50;
FOLDS = 5;
SEED_START = 20260617;
BLOCK_SIZE = 256;
TOL = 1e-12;

D(1).name = 'COIL100';
D(1).files = {'COIL100_Obj(1).mat', 'COIL100_Obj.mat', 'COIL100.mat'};
D(1).p = 3; D(1).alpha = 10; D(1).expectedSamples = 7200; D(1).expectedClasses = 100;
D(2).name = 'PIE';
D(2).files = {'PIE.mat'};
D(2).p = 3; D(2).alpha = 1000; D(2).expectedSamples = 2856; D(2).expectedClasses = 68;
D(3).name = 'YaleB';
D(3).files = {'YaleB.mat', 'YaleB(1).mat'};
D(3).p = 2; D(3).alpha = 1000; D(3).expectedSamples = 2414; D(3).expectedClasses = 38;

rawPath = fullfile(outputDir, sprintf('GCGOCNMF_min_ablation3_%dseed_raw.csv', runs));
summaryPath = fullfile(outputDir, sprintf('GCGOCNMF_min_ablation3_%dseed_summary.csv', runs));
decisionPath = fullfile(outputDir, sprintf('GCGOCNMF_min_ablation3_%dseed_decision.txt', runs));
matPath = fullfile(outputDir, sprintf('GCGOCNMF_min_ablation3_%dseed_results.mat', runs));

fid = fopen(rawPath, 'w');
if fid < 0, error('Cannot create %s.', rawPath); end
fprintf(fid, ['dataset,seed,file,samples,features,classes,p,alpha,labeled_count,' ...
    'theta,delta_cv,se_cv,binary_cv_brier,local_cv_brier,hard_choice,' ...
    'base_ACC,local_ACC,hard_ACC,shrink_ACC,' ...
    'local_dACC,hard_dACC,shrink_dACC,' ...
    'base_NMI,local_NMI,hard_NMI,shrink_NMI,' ...
    'local_dNMI,hard_dNMI,shrink_dNMI,' ...
    'graph_seconds,cv_seconds,base_seconds,local_seconds,shrink_seconds,total_seconds\n']);
fclose(fid);

summaryTemplate = struct( ...
    'dataset','', 'runs',0, 'theta_mean',0, 'theta_min',0, 'theta_max',0, ...
    'hard_local_count',0, 'hard_binary_count',0, ...
    'base_ACC_mean',0, 'local_ACC_mean',0, 'hard_ACC_mean',0, 'shrink_ACC_mean',0, ...
    'local_dACC_mean',0, 'hard_dACC_mean',0, 'shrink_dACC_mean',0, ...
    'local_ACC_losses',0, 'hard_ACC_losses',0, 'shrink_ACC_losses',0, ...
    'base_NMI_mean',0, 'local_NMI_mean',0, 'hard_NMI_mean',0, 'shrink_NMI_mean',0, ...
    'local_dNMI_mean',0, 'hard_dNMI_mean',0, 'shrink_dNMI_mean',0, ...
    'local_NMI_losses',0, 'hard_NMI_losses',0, 'shrink_NMI_losses',0, ...
    'shrink_dual_wins',0, 'shrink_dual_losses',0, ...
    'cv_delta_mean',0, 'cv_se_mean',0, ...
    'graph_seconds',0, 'cv_seconds_mean',0, 'base_seconds_mean',0, ...
    'local_seconds_mean',0, 'shrink_seconds_mean',0);
summaryRows = repmat(summaryTemplate, 1, length(D));

fprintf('\nRUNNING MINIMAL THREE-DATASET GCGOCNMF ABLATION -- 2026-07-30\n');
fprintf('Runs per dataset: %d\n', runs);
fprintf('Data directory: %s\n', dataDir);
fprintf('Output directory: %s\n\n', outputDir);

for d = 1:length(D)
    dataPath = gcgoc_find_file(dataDir, D(d).files);
    fprintf('Loading %s: %s\n', D(d).name, dataPath);
    [X, y] = gcgoc_load_dataset(dataPath, D(d).expectedSamples, D(d).expectedClasses);
    n = size(X,2); m = size(X,1); c = length(unique(y));

    graphClock = tic;
    [W0, W1] = gcgoc_construct_graphs(X, D(d).p, BLOCK_SIZE, EPSILON);
    graphSeconds = toc(graphClock);

    thetaVec=zeros(runs,1); deltaVec=zeros(runs,1); seVec=zeros(runs,1);
    hardChoiceVec=zeros(runs,1);
    baseACC=zeros(runs,1); localACC=zeros(runs,1); hardACC=zeros(runs,1); shrinkACC=zeros(runs,1);
    baseNMI=zeros(runs,1); localNMI=zeros(runs,1); hardNMI=zeros(runs,1); shrinkNMI=zeros(runs,1);
    cvTime=zeros(runs,1); baseTime=zeros(runs,1); localTime=zeros(runs,1); shrinkTime=zeros(runs,1);

    for r = 1:runs
        seed = SEED_START + r - 1;
        totalClock = tic;
        L = gcgoc_labeled_indices(y, seed, LABEL_FRACTION);

        cvClock = tic;
        [theta, deltaCV, seCV, loss0, loss1] = gcgoc_shrinkage_weight( ...
            W0, W1, L, y, c, seed, FOLDS, EPSILON);
        cvTime(r) = toc(cvClock);
        if deltaCV > 0
            hardChoice = 1;
        else
            hardChoice = 0;
        end
        Wshrink = (1-theta)*W0 + theta*W1;

        t0=tic; V0=gcgoc_train(X,y,L,W0,D(d).alpha,seed,ITERATIONS,EPSILON); baseTime(r)=toc(t0);
        t1=tic; VL=gcgoc_train(X,y,L,W1,D(d).alpha,seed,ITERATIONS,EPSILON); localTime(r)=toc(t1);
        [baseACC(r),baseNMI(r)] = gcgoc_evaluate(V0,y,L,EPSILON);
        [localACC(r),localNMI(r)] = gcgoc_evaluate(VL,y,L,EPSILON);
        if hardChoice == 1
            hardACC(r)=localACC(r); hardNMI(r)=localNMI(r);
        else
            hardACC(r)=baseACC(r); hardNMI(r)=baseNMI(r);
        end
        if theta <= TOL
            shrinkACC(r)=baseACC(r); shrinkNMI(r)=baseNMI(r); shrinkTime(r)=0;
        elseif abs(theta-1) <= TOL
            shrinkACC(r)=localACC(r); shrinkNMI(r)=localNMI(r); shrinkTime(r)=0;
        else
            ts=tic; VS=gcgoc_train(X,y,L,Wshrink,D(d).alpha,seed,ITERATIONS,EPSILON); shrinkTime(r)=toc(ts);
            [shrinkACC(r),shrinkNMI(r)] = gcgoc_evaluate(VS,y,L,EPSILON);
        end

        thetaVec(r)=theta; deltaVec(r)=deltaCV; seVec(r)=seCV; hardChoiceVec(r)=hardChoice;
        totalSeconds=toc(totalClock);
        [pathOnly,fileOnly,extOnly]=fileparts(dataPath); %#ok<ASGLU>
        fileNameOnly=[fileOnly extOnly];
        fid=fopen(rawPath,'a'); if fid<0, error('Cannot append to %s.',rawPath); end
        fprintf(fid, ['%s,%d,%s,%d,%d,%d,%d,%.15g,%d,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,%d,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
            '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n'], ...
            D(d).name,seed,fileNameOnly,n,m,c,D(d).p,D(d).alpha,length(L), ...
            theta,deltaCV,seCV,loss0,loss1,hardChoice, ...
            baseACC(r),localACC(r),hardACC(r),shrinkACC(r), ...
            localACC(r)-baseACC(r),hardACC(r)-baseACC(r),shrinkACC(r)-baseACC(r), ...
            baseNMI(r),localNMI(r),hardNMI(r),shrinkNMI(r), ...
            localNMI(r)-baseNMI(r),hardNMI(r)-baseNMI(r),shrinkNMI(r)-baseNMI(r), ...
            graphSeconds,cvTime(r),baseTime(r),localTime(r),shrinkTime(r),totalSeconds);
        fclose(fid);
        fprintf('%s seed=%d theta=%.4f local dACC=%+.5f hard=%+.5f shrink=%+.5f\n', ...
            D(d).name,seed,theta,localACC(r)-baseACC(r),hardACC(r)-baseACC(r),shrinkACC(r)-baseACC(r));
    end

    summaryRows(d).dataset=D(d).name; summaryRows(d).runs=runs;
    summaryRows(d).theta_mean=mean(thetaVec); summaryRows(d).theta_min=min(thetaVec); summaryRows(d).theta_max=max(thetaVec);
    summaryRows(d).hard_local_count=sum(hardChoiceVec==1); summaryRows(d).hard_binary_count=sum(hardChoiceVec==0);
    summaryRows(d).base_ACC_mean=mean(baseACC); summaryRows(d).local_ACC_mean=mean(localACC);
    summaryRows(d).hard_ACC_mean=mean(hardACC); summaryRows(d).shrink_ACC_mean=mean(shrinkACC);
    summaryRows(d).local_dACC_mean=mean(localACC-baseACC); summaryRows(d).hard_dACC_mean=mean(hardACC-baseACC);
    summaryRows(d).shrink_dACC_mean=mean(shrinkACC-baseACC);
    summaryRows(d).local_ACC_losses=sum(localACC-baseACC < -TOL); summaryRows(d).hard_ACC_losses=sum(hardACC-baseACC < -TOL);
    summaryRows(d).shrink_ACC_losses=sum(shrinkACC-baseACC < -TOL);
    summaryRows(d).base_NMI_mean=mean(baseNMI); summaryRows(d).local_NMI_mean=mean(localNMI);
    summaryRows(d).hard_NMI_mean=mean(hardNMI); summaryRows(d).shrink_NMI_mean=mean(shrinkNMI);
    summaryRows(d).local_dNMI_mean=mean(localNMI-baseNMI); summaryRows(d).hard_dNMI_mean=mean(hardNMI-baseNMI);
    summaryRows(d).shrink_dNMI_mean=mean(shrinkNMI-baseNMI);
    summaryRows(d).local_NMI_losses=sum(localNMI-baseNMI < -TOL); summaryRows(d).hard_NMI_losses=sum(hardNMI-baseNMI < -TOL);
    summaryRows(d).shrink_NMI_losses=sum(shrinkNMI-baseNMI < -TOL);
    summaryRows(d).shrink_dual_wins=sum((shrinkACC-baseACC>TOL)&(shrinkNMI-baseNMI>TOL));
    summaryRows(d).shrink_dual_losses=sum((shrinkACC-baseACC<-TOL)|(shrinkNMI-baseNMI<-TOL));
    summaryRows(d).cv_delta_mean=mean(deltaVec); summaryRows(d).cv_se_mean=mean(seVec);
    summaryRows(d).graph_seconds=graphSeconds; summaryRows(d).cv_seconds_mean=mean(cvTime);
    summaryRows(d).base_seconds_mean=mean(baseTime); summaryRows(d).local_seconds_mean=mean(localTime);
    summaryRows(d).shrink_seconds_mean=mean(shrinkTime);
    save(matPath,'summaryRows');
end

fid=fopen(summaryPath,'w'); if fid<0, error('Cannot create %s.',summaryPath); end
fprintf(fid, ['dataset,runs,theta_mean,theta_min,theta_max,hard_local_count,hard_binary_count,' ...
    'base_ACC_mean,local_ACC_mean,hard_ACC_mean,shrink_ACC_mean,' ...
    'local_dACC_mean,hard_dACC_mean,shrink_dACC_mean,local_ACC_losses,hard_ACC_losses,shrink_ACC_losses,' ...
    'base_NMI_mean,local_NMI_mean,hard_NMI_mean,shrink_NMI_mean,' ...
    'local_dNMI_mean,hard_dNMI_mean,shrink_dNMI_mean,local_NMI_losses,hard_NMI_losses,shrink_NMI_losses,' ...
    'shrink_dual_wins,shrink_dual_losses,cv_delta_mean,cv_se_mean,' ...
    'graph_seconds,cv_seconds_mean,base_seconds_mean,local_seconds_mean,shrink_seconds_mean\n']);
for d=1:length(summaryRows)
    S=summaryRows(d);
    fprintf(fid, ['%s,%d,%.15g,%.15g,%.15g,%d,%d,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d,%d,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d,%d,%d,%d,' ...
        '%.15g,%.15g,%.6f,%.6f,%.6f,%.6f,%.6f\n'], ...
        S.dataset,S.runs,S.theta_mean,S.theta_min,S.theta_max,S.hard_local_count,S.hard_binary_count, ...
        S.base_ACC_mean,S.local_ACC_mean,S.hard_ACC_mean,S.shrink_ACC_mean, ...
        S.local_dACC_mean,S.hard_dACC_mean,S.shrink_dACC_mean,S.local_ACC_losses,S.hard_ACC_losses,S.shrink_ACC_losses, ...
        S.base_NMI_mean,S.local_NMI_mean,S.hard_NMI_mean,S.shrink_NMI_mean, ...
        S.local_dNMI_mean,S.hard_dNMI_mean,S.shrink_dNMI_mean,S.local_NMI_losses,S.hard_NMI_losses,S.shrink_NMI_losses, ...
        S.shrink_dual_wins,S.shrink_dual_losses,S.cv_delta_mean,S.cv_se_mean, ...
        S.graph_seconds,S.cv_seconds_mean,S.base_seconds_mean,S.local_seconds_mean,S.shrink_seconds_mean);
end
fclose(fid);

finitePass=1; shrinkSafetyPass=1;
for d=1:length(summaryRows)
    S=summaryRows(d);
    values=[S.theta_mean S.theta_min S.theta_max S.cv_delta_mean S.cv_se_mean ...
        S.base_ACC_mean S.local_ACC_mean S.hard_ACC_mean S.shrink_ACC_mean ...
        S.base_NMI_mean S.local_NMI_mean S.hard_NMI_mean S.shrink_NMI_mean];
    if any(~isfinite(values)), finitePass=0; end
    if S.shrink_ACC_losses>0 || S.shrink_NMI_losses>0, shrinkSafetyPass=0; end
end

% The mechanism gate is descriptive and was fixed before viewing results:
% (i) SHRINK must improve both metrics on COIL100;
% (ii) LOCAL must damage at least one metric on PIE or YaleB;
% (iii) SHRINK must have no ACC/NMI losses on any of the three datasets.
mechanismPass=0;
coilIndex=find(strcmp({summaryRows.dataset},'COIL100'));
pieIndex=find(strcmp({summaryRows.dataset},'PIE'));
yaleIndex=find(strcmp({summaryRows.dataset},'YaleB'));
if length(coilIndex)==1 && length(pieIndex)==1 && length(yaleIndex)==1
    coilGain=(summaryRows(coilIndex).shrink_dACC_mean>0) && ...
        (summaryRows(coilIndex).shrink_dNMI_mean>0);
    localDamage=(summaryRows(pieIndex).local_dACC_mean<0) || ...
        (summaryRows(pieIndex).local_dNMI_mean<0) || ...
        (summaryRows(yaleIndex).local_dACC_mean<0) || ...
        (summaryRows(yaleIndex).local_dNMI_mean<0);
    mechanismPass=finitePass && shrinkSafetyPass && coilGain && localDamage;
end
fid=fopen(decisionPath,'w'); if fid<0, error('Cannot create %s.',decisionPath); end
fprintf(fid,'Minimal three-dataset GCGOCNMF ablation\nRuns per dataset: %d\n\n',runs);
for d=1:length(summaryRows)
    S=summaryRows(d);
    fprintf(fid,['%s: theta=%.6f, dACC local/hard/shrink=%+.6f/%+.6f/%+.6f, ' ...
        'ACC losses=%d/%d/%d, dNMI local/hard/shrink=%+.6f/%+.6f/%+.6f, NMI losses=%d/%d/%d\n'], ...
        S.dataset,S.theta_mean,S.local_dACC_mean,S.hard_dACC_mean,S.shrink_dACC_mean, ...
        S.local_ACC_losses,S.hard_ACC_losses,S.shrink_ACC_losses, ...
        S.local_dNMI_mean,S.hard_dNMI_mean,S.shrink_dNMI_mean, ...
        S.local_NMI_losses,S.hard_NMI_losses,S.shrink_NMI_losses);
end
fprintf(fid,'\nCV_FINITE_PASS=%d\nSHRINK_SAFETY_PASS=%d\nABLATION_MECHANISM_PASS=%d\n',finitePass,shrinkSafetyPass,mechanismPass);
fclose(fid);
save(matPath,'summaryRows','finitePass','shrinkSafetyPass','mechanismPass');

fprintf('\nFinished.\nRaw: %s\nSummary: %s\nDecision: %s\n',rawPath,summaryPath,decisionPath);
fprintf('CV_FINITE_PASS=%d\nSHRINK_SAFETY_PASS=%d\nABLATION_MECHANISM_PASS=%d\n',finitePass,shrinkSafetyPass,mechanismPass);
end


function path = gcgoc_find_file(folder, candidates)
path = '';
for i = 1:length(candidates)
    candidate = fullfile(folder, candidates{i});
    if exist(candidate, 'file') == 2
        path = candidate;
        return;
    end
end
msg = candidates{1};
for i = 2:length(candidates)
    msg = [msg ', ' candidates{i}]; %#ok<AGROW>
end
error('Missing dataset in %s. Accepted names: %s', folder, msg);
end


function [X, y] = gcgoc_load_dataset(path, expectedSamples, expectedClasses)
S = load(path);
featureNames = {'fea', 'X', 'data', 'features'};
labelNames = {'gnd', 'labels', 'label', 'y', 'Y', 'truth'};
fea = [];
gnd = [];
for i = 1:length(featureNames)
    if isfield(S, featureNames{i})
        fea = S.(featureNames{i});
        break;
    end
end
for i = 1:length(labelNames)
    if isfield(S, labelNames{i})
        gnd = S.(labelNames{i});
        break;
    end
end
if isempty(fea)
    error('%s: no feature field found.', path);
end
if isempty(gnd)
    error('%s: no label field found.', path);
end

y0 = double(gnd(:));
F = double(fea);
if ndims(F) ~= 2
    error('%s: the feature array must be two-dimensional.', path);
end
if size(F,1) == length(y0)
    X = F';
elseif size(F,2) == length(y0)
    X = F;
else
    error('%s: feature shape is incompatible with the label count.', path);
end
if min(X(:)) < -1e-12
    error('%s: negative feature values were found; standard NMF requires nonnegative data.', path);
end

[values, dummyIndex, y] = unique(y0); %#ok<ASGLU>
y = double(y(:));
if size(X,2) ~= expectedSamples || length(values) ~= expectedClasses
    error('%s: expected %d samples/%d classes, obtained %d/%d.', ...
        path, expectedSamples, expectedClasses, size(X,2), length(values));
end

norms = sqrt(sum(X.^2, 1));
norms(norms == 0) = 1;
X = bsxfun(@rdivide, X, norms);
end


function gcgoc_set_seed(seed)
% Legacy syntax is used for MATLAB R2009a compatibility.
rand('twister', double(seed)); %#ok<RAND>
end


function L = gcgoc_labeled_indices(y, seed, fraction)
gcgoc_set_seed(seed);
classes = unique(y);
L = zeros(0,1);
for k = 1:length(classes)
    ids = find(y == classes(k));
    order = randperm(length(ids));
    count = max(2, floor(fraction * length(ids)));
    count = min(count, length(ids));
    L = [L; ids(order(1:count))]; %#ok<AGROW>
end
L = double(L(:));
end


function [W0, W1] = gcgoc_construct_graphs(X, p, blockSize, epsilon)
% X is features-by-samples and has unit-norm columns.
n = size(X,2);
rows = zeros(n*p,1);
cols = zeros(n*p,1);
dvals = zeros(n*p,1);
sigma = zeros(n,1);
position = 0;

for first = 1:blockSize:n
    last = min(first + blockSize - 1, n);
    ids = first:last;
    similarities = X(:,ids)' * X;
    for q = 1:length(ids)
        similarities(q, ids(q)) = -Inf;
    end
    [sortedSimilarity, sortedIndex] = sort(similarities, 2, 'descend');
    neighborIndex = sortedIndex(:,1:p);
    neighborDistance = max(0, 1 - sortedSimilarity(:,1:p));

    rowBlock = repmat(ids(:), 1, p);
    count = length(ids) * p;
    range = position + (1:count);
    rows(range) = rowBlock(:);
    cols(range) = neighborIndex(:);
    dvals(range) = neighborDistance(:);
    sigma(ids) = max(neighborDistance(:,p), epsilon);
    position = position + count;

    fprintf('  kNN block %d:%d of %d\n', first, last, n);
    clear similarities sortedSimilarity sortedIndex neighborIndex neighborDistance rowBlock;
end

A = sparse(rows, cols, ones(length(rows),1), n, n);
W0 = spones(A + A');
W0 = W0 - spdiags(diag(W0), 0, n, n);
W0 = sparse(W0);

denominator = sqrt(sigma(rows) .* sigma(cols)) + epsilon;
weights = exp(-(dvals ./ denominator).^2);
Als = sparse(rows, cols, weights, n, n);
W1 = max(Als, Als');
W1 = W1 - spdiags(diag(W1), 0, n, n);
nonzeroValues = nonzeros(W1);
if ~isempty(nonzeroValues)
    W1 = W1 / mean(nonzeroValues);
end
W1 = sparse(W1);
end


function folds = gcgoc_stratified_folds(L, y, seed, numberOfFolds)
gcgoc_set_seed(seed);
folds = cell(numberOfFolds,1);
for f = 1:numberOfFolds
    folds{f} = zeros(0,1);
end
classes = unique(y(L));
for k = 1:length(classes)
    ids = L(y(L) == classes(k));
    order = randperm(length(ids));
    ids = ids(order);
    % Rotate the starting fold across classes. This is essential when a
    % class has fewer labeled samples than the requested number of folds
    % Offset class-wise folds to avoid an empty global fold when a class has fewer than five labeled samples.
    % Without the rotation, the fifth fold is globally empty and its
    % Brier loss becomes NaN.
    classOffset = mod(k-1, numberOfFolds);
    for j = 1:length(ids)
        f = mod((j-1) + classOffset, numberOfFolds) + 1;
        folds{f} = [folds{f}; ids(j)]; %#ok<AGROW>
    end
end
for f = 1:numberOfFolds
    if isempty(folds{f})
        error('Cross-validation fold %d is empty. Reduce FOLDS or inspect the labeled split.', f);
    end
end
end


function P = gcgoc_harmonic_predictions(W, train, y, c, query, epsilon)
n = size(W,1);
unknownMask = true(n,1);
unknownMask(train) = false;
unknown = find(unknownMask);
position = zeros(n,1);
position(unknown) = 1:length(unknown);
degree = full(sum(W,2));
A = spdiags(degree(unknown), 0, length(unknown), length(unknown)) - W(unknown,unknown);
C = zeros(length(train), c);
index = sub2ind(size(C), (1:length(train))', y(train));
C(index) = 1;
rhs = W(unknown,train) * C;

lastwarn('');
try
    Fu = A \ rhs;
catch
    Fu = (A + 1e-10 * speye(length(unknown))) \ rhs;
end
[warningMessage, warningId] = lastwarn; %#ok<ASGLU>
if any(~isfinite(Fu(:))) || ~isempty(strfind(lower(warningMessage), 'singular'))
    Fu = (A + 1e-10 * speye(length(unknown))) \ rhs;
end
Fu = max(Fu, 0);
rowSum = sum(Fu,2);
positive = rowSum > epsilon;
if any(positive)
    Fu(positive,:) = bsxfun(@rdivide, Fu(positive,:), rowSum(positive));
end
P = Fu(position(query),:);
end


function losses = gcgoc_fold_brier_losses(W, folds, L, y, c, epsilon)
numberOfFolds = length(folds);
losses = zeros(numberOfFolds,1);
for f = 1:numberOfFolds
    heldOut = folds{f};
    if isempty(heldOut)
        error('Cross-validation fold %d is empty.', f);
    end
    train = setdiff(L, heldOut);
    prediction = gcgoc_harmonic_predictions(W, train, y, c, heldOut, epsilon);
    if any(~isfinite(prediction(:)))
        error('Non-finite harmonic prediction detected in fold %d.', f);
    end
    target = zeros(length(heldOut), c);
    index = sub2ind(size(target), (1:length(heldOut))', y(heldOut));
    target(index) = 1;
    losses(f) = mean(sum((prediction - target).^2, 2));
end
end


function [theta, delta, se, loss0Mean, loss1Mean] = gcgoc_shrinkage_weight( ...
    W0, W1, L, y, c, seed, numberOfFolds, epsilon)
folds = gcgoc_stratified_folds(L, y, seed + 991, numberOfFolds);
loss0 = gcgoc_fold_brier_losses(W0, folds, L, y, c, epsilon);
loss1 = gcgoc_fold_brier_losses(W1, folds, L, y, c, epsilon);
if any(~isfinite(loss0)) || any(~isfinite(loss1))
    error('Non-finite cross-validation Brier loss detected.');
end
difference = loss0 - loss1;
delta = mean(difference);
if length(difference) > 1
    se = std(difference, 0) / sqrt(length(difference));
else
    se = 0;
end
if delta > 0
    theta = max(0, 1 - se / max(delta, epsilon));
else
    theta = 0;
end
theta = min(1, theta);
loss0Mean = mean(loss0);
loss1Mean = mean(loss1);
end


function V = gcgoc_train(X, y, L, W, alpha, seed, iterations, epsilon)
[m, n] = size(X);
c = length(unique(y));
gcgoc_set_seed(seed);
U = max(rand(m,c), epsilon);
V = max(rand(n,c), epsilon);
C = zeros(length(L), c);
index = sub2ind(size(C), (1:length(L))', y(L));
C(index) = 1;
V(L,:) = C;
degree = full(sum(W,2));

for iter = 1:iterations
    numeratorU = X * V;
    denominatorU = U * (V' * V);
    U = U .* (numeratorU ./ max(denominatorU, epsilon));

    numeratorV = X' * U + alpha * (W * V);
    denominatorV = V * (U' * U) + alpha * bsxfun(@times, degree, V);
    V = V .* (numeratorV ./ max(denominatorV, epsilon));
    V(L,:) = C;
end
end


function [acc, nmi] = gcgoc_evaluate(V, y, L, epsilon)
n = length(y);
mask = true(n,1);
mask(L) = false;
U = find(mask);
[dummy, prediction] = max(V, [], 2); %#ok<ASGLU>
acc = mean(prediction(U) == y(U));
nmi = gcgoc_nmi(y(U), prediction(U), epsilon);
end


function value = gcgoc_nmi(trueLabel, predictedLabel, epsilon)
trueLabel = trueLabel(:);
predictedLabel = predictedLabel(:);
[trueValues, dummy1, trueIndex] = unique(trueLabel); %#ok<ASGLU>
[predValues, dummy2, predIndex] = unique(predictedLabel); %#ok<ASGLU>
nt = length(trueValues);
np = length(predValues);
N = length(trueLabel);
contingency = accumarray([trueIndex predIndex], 1, [nt np]);
pij = contingency / N;
pi = sum(pij, 2);
pj = sum(pij, 1);

mi = 0;
for i = 1:nt
    for j = 1:np
        if pij(i,j) > 0
            mi = mi + pij(i,j) * log(pij(i,j) / max(pi(i)*pj(j), epsilon));
        end
    end
end
ht = -sum(pi(pi>0) .* log(pi(pi>0)));
hp = -sum(pj(pj>0) .* log(pj(pj>0)));
denominator = 0.5 * (ht + hp);
if denominator <= epsilon
    value = 1;
else
    value = mi / denominator;
end
value = max(0, min(1, value));
end


function gcgoc_write_summary_csv(path, rows)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot create %s.', path);
end
fprintf(fid, ['dataset,runs,theta_mean,theta_min,theta_max,' ...
    'base_ACC_mean,model_ACC_mean,dACC_mean,dACC_sd,ACC_wins,ACC_ties,ACC_losses,' ...
    'base_NMI_mean,model_NMI_mean,dNMI_mean,dNMI_sd,NMI_wins,NMI_ties,NMI_losses,' ...
    'dual_wins,cv_delta_mean,cv_se_mean\n']);
for d = 1:length(rows)
    S = rows(d);
    fprintf(fid, ['%s,%d,%.15g,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%.15g,%d,%d,%d,' ...
        '%.15g,%.15g,%.15g,%.15g,%d,%d,%d,%d,%.15g,%.15g\n'], ...
        S.dataset, S.runs, S.theta_mean, S.theta_min, S.theta_max, ...
        S.base_ACC_mean, S.model_ACC_mean, S.dACC_mean, S.dACC_sd, ...
        S.ACC_wins, S.ACC_ties, S.ACC_losses, ...
        S.base_NMI_mean, S.model_NMI_mean, S.dNMI_mean, S.dNMI_sd, ...
        S.NMI_wins, S.NMI_ties, S.NMI_losses, S.dual_wins, ...
        S.cv_delta_mean, S.cv_se_mean);
end
fclose(fid);
end


function gcgoc_write_decision(path, runs, rows, passFlag, cvFiniteFlag)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot create %s.', path);
end
fprintf(fid, 'Conservative graph-calibrated GOCNMF: first-three MATLAB replication\n');
fprintf(fid, 'Runs per dataset: %d\n', runs);
fprintf(fid, 'Formula and protocol are frozen.\n\n');
for d = 1:length(rows)
    S = rows(d);
    fprintf(fid, ['%s: theta_mean=%.8f, base_ACC=%.8f, model_ACC=%.8f, dACC=%.8f, ' ...
        'ACC_wins/ties/losses=%d/%d/%d, base_NMI=%.8f, model_NMI=%.8f, ' ...
        'dNMI=%.8f, NMI_wins/ties/losses=%d/%d/%d, dual_wins=%d\n'], ...
        S.dataset, S.theta_mean, S.base_ACC_mean, S.model_ACC_mean, S.dACC_mean, ...
        S.ACC_wins, S.ACC_ties, S.ACC_losses, S.base_NMI_mean, S.model_NMI_mean, ...
        S.dNMI_mean, S.NMI_wins, S.NMI_ties, S.NMI_losses, S.dual_wins);
end
fprintf(fid, '\nCV_FINITE_PASS=%d\n', cvFiniteFlag);
if runs == 3
    fprintf(fid, 'THREE_SEED_INTEGRITY_PASS=%d\n', passFlag);
    fprintf(fid, 'Run 20 seeds only when this value is 1.\n');
else
    fprintf(fid, 'FIRST_THREE_20SEED_PASS=%d\n', passFlag);
end
fclose(fid);
end
