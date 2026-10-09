function run_GCGOCNMF_first3_R2009a_FIX4(runs, dataDir, outputDir)
%RUN_GCGOCNMF_FIRST3_R2009A
% MATLAB R2009a-compatible validation program for the conservative
% graph-calibrated GOCNMF candidate on COIL20, Optdigits, and YaleB.
%
% Usage:
%   run_GCGOCNMF_first3_R2009a_FIX4(3,  pwd, pwd);
%   run_GCGOCNMF_first3_R2009a_FIX4(20, pwd, pwd);
%
% Run the 3-seed gate first. Run 20 seeds only when the decision file says
% THREE_SEED_INTEGRITY_PASS=1.
%
% Required data files (accepted alternatives are listed below):
%   COIL20    : COIL20_Obj.mat / COIL20.mat / COIL20_Obj(1).mat
%   Optdigits : Optdigits_Han.mat / Optdigits.mat
%   YaleB     : YaleB.mat / YaleB(1).mat
%
% Supported feature fields: fea, X, data, features
% Supported label fields  : gnd, labels, label, y, Y, truth
%
% The program keeps the protocol fixed:
%   label rate = 10%, five folds, 50 NMF iterations,
%   seeds = 20260617 onward,
%   COIL20    p=3, alpha=10,
%   Optdigits p=4, alpha=10,
%   YaleB     p=2, alpha=1000.
%
% No K-means or external classifier is used. Prediction is direct argmax.

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

D(1).name = 'COIL20';
D(1).files = {'COIL20_Obj.mat', 'COIL20.mat', 'COIL20_Obj(1).mat'};
D(1).p = 3;
D(1).alpha = 10;
D(1).expectedSamples = 1440;
D(1).expectedClasses = 20;

D(2).name = 'Optdigits';
D(2).files = {'Optdigits_Han.mat', 'Optdigits.mat'};
D(2).p = 4;
D(2).alpha = 10;
D(2).expectedSamples = 5620;
D(2).expectedClasses = 10;

D(3).name = 'YaleB';
D(3).files = {'YaleB.mat', 'YaleB(1).mat'};
D(3).p = 2;
D(3).alpha = 1000;
D(3).expectedSamples = 2414;
D(3).expectedClasses = 38;

rawPath = fullfile(outputDir, sprintf('graph_calibrated_first3_FIX4_%dseed_raw.csv', runs));
summaryPath = fullfile(outputDir, sprintf('graph_calibrated_first3_FIX4_%dseed_summary.csv', runs));
decisionPath = fullfile(outputDir, sprintf('graph_calibrated_first3_FIX4_%dseed_decision.txt', runs));
matPath = fullfile(outputDir, sprintf('graph_calibrated_first3_FIX4_%dseed_results.mat', runs));

fid = fopen(rawPath, 'w');
if fid < 0
    error('Cannot create %s.', rawPath);
end
fprintf(fid, ['dataset,seed,file,samples,features,classes,p,alpha,labeled_count,' ...
    'theta,delta_cv,se_cv,binary_cv_brier,local_cv_brier,' ...
    'base_ACC,ACC,dACC,base_NMI,NMI,dNMI,graph_seconds,run_seconds\n']);
fclose(fid);

rawRows = struct([]);
summaryTemplate = struct( ...
    'dataset', '', 'runs', 0, ...
    'theta_mean', 0, 'theta_min', 0, 'theta_max', 0, ...
    'base_ACC_mean', 0, 'model_ACC_mean', 0, ...
    'dACC_mean', 0, 'dACC_sd', 0, ...
    'ACC_wins', 0, 'ACC_ties', 0, 'ACC_losses', 0, ...
    'base_NMI_mean', 0, 'model_NMI_mean', 0, ...
    'dNMI_mean', 0, 'dNMI_sd', 0, ...
    'NMI_wins', 0, 'NMI_ties', 0, 'NMI_losses', 0, ...
    'dual_wins', 0, 'cv_delta_mean', 0, 'cv_se_mean', 0);
summaryRows = repmat(summaryTemplate, 1, length(D));
rowCounter = 0;

fprintf('\nRUNNING FIRST3 FIX4 MATLAB REPLICATION -- 2026-07-30\n');
fprintf('\nConservative graph-calibrated GOCNMF validation\n');
fprintf('Runs per dataset: %d\n', runs);
fprintf('Data directory: %s\n', dataDir);
fprintf('Output directory: %s\n\n', outputDir);

for d = 1:length(D)
    dataPath = gcgoc_find_file(dataDir, D(d).files);
    fprintf('Loading %s: %s\n', D(d).name, dataPath);
    [X, y] = gcgoc_load_dataset(dataPath, D(d).expectedSamples, D(d).expectedClasses);
    n = size(X, 2);
    m = size(X, 1);
    c = length(unique(y));

    graphClock = tic;
    [W0, W1] = gcgoc_construct_graphs(X, D(d).p, BLOCK_SIZE, EPSILON);
    graphSeconds = toc(graphClock);
    fprintf('%s graph built in %.2f seconds; nnz(W0)=%d, nnz(W1)=%d\n', ...
        D(d).name, graphSeconds, nnz(W0), nnz(W1));

    thetaVec = zeros(runs,1);
    deltaVec = zeros(runs,1);
    seVec = zeros(runs,1);
    baseACCVec = zeros(runs,1);
    accVec = zeros(runs,1);
    dACCVec = zeros(runs,1);
    baseNMIVec = zeros(runs,1);
    nmiVec = zeros(runs,1);
    dNMIVec = zeros(runs,1);

    for r = 1:runs
        seed = SEED_START + r - 1;
        runClock = tic;

        L = gcgoc_labeled_indices(y, seed, LABEL_FRACTION);
        [theta, deltaCV, seCV, loss0, loss1] = gcgoc_shrinkage_weight( ...
            W0, W1, L, y, c, seed, FOLDS, EPSILON);
        Wtheta = (1-theta) * W0 + theta * W1;

        V0 = gcgoc_train(X, y, L, W0, D(d).alpha, seed, ITERATIONS, EPSILON);
        V1 = gcgoc_train(X, y, L, Wtheta, D(d).alpha, seed, ITERATIONS, EPSILON);

        [baseACC, baseNMI] = gcgoc_evaluate(V0, y, L, EPSILON);
        [modelACC, modelNMI] = gcgoc_evaluate(V1, y, L, EPSILON);
        dACC = modelACC - baseACC;
        dNMI = modelNMI - baseNMI;
        runSeconds = toc(runClock);

        thetaVec(r) = theta;
        deltaVec(r) = deltaCV;
        seVec(r) = seCV;
        baseACCVec(r) = baseACC;
        accVec(r) = modelACC;
        dACCVec(r) = dACC;
        baseNMIVec(r) = baseNMI;
        nmiVec(r) = modelNMI;
        dNMIVec(r) = dNMI;

        rowCounter = rowCounter + 1;
        rawRows(rowCounter).dataset = D(d).name;
        rawRows(rowCounter).seed = seed;
        rawRows(rowCounter).file = dataPath;
        rawRows(rowCounter).samples = n;
        rawRows(rowCounter).features = m;
        rawRows(rowCounter).classes = c;
        rawRows(rowCounter).p = D(d).p;
        rawRows(rowCounter).alpha = D(d).alpha;
        rawRows(rowCounter).labeled_count = length(L);
        rawRows(rowCounter).theta = theta;
        rawRows(rowCounter).delta_cv = deltaCV;
        rawRows(rowCounter).se_cv = seCV;
        rawRows(rowCounter).binary_cv_brier = loss0;
        rawRows(rowCounter).local_cv_brier = loss1;
        rawRows(rowCounter).base_ACC = baseACC;
        rawRows(rowCounter).ACC = modelACC;
        rawRows(rowCounter).dACC = dACC;
        rawRows(rowCounter).base_NMI = baseNMI;
        rawRows(rowCounter).NMI = modelNMI;
        rawRows(rowCounter).dNMI = dNMI;
        rawRows(rowCounter).graph_seconds = graphSeconds;
        rawRows(rowCounter).run_seconds = runSeconds;

        fid = fopen(rawPath, 'a');
        if fid < 0
            error('Cannot append to %s.', rawPath);
        end
        [pathOnly, fileOnly, extOnly] = fileparts(dataPath);
        fileNameOnly = [fileOnly extOnly];
        fprintf(fid, ['%s,%d,%s,%d,%d,%d,%d,%.15g,%d,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.6f,%.6f\n'], ...
            D(d).name, seed, fileNameOnly, n, m, c, D(d).p, D(d).alpha, length(L), ...
            theta, deltaCV, seCV, loss0, loss1, ...
            baseACC, modelACC, dACC, baseNMI, modelNMI, dNMI, graphSeconds, runSeconds);
        fclose(fid);

        save(matPath, 'rawRows', 'summaryRows');
        fprintf('%s seed=%d theta=%.6f dACC=%+.6f dNMI=%+.6f time=%.1fs\n', ...
            D(d).name, seed, theta, dACC, dNMI, runSeconds);
    end

    % Assign fields directly.  In MATLAB R2009a, assigning a populated
    % structure to struct([]) can raise
    % "Subscripted assignment between dissimilar structures".
    summaryRows(d).dataset = D(d).name;
    summaryRows(d).runs = runs;
    summaryRows(d).theta_mean = mean(thetaVec);
    summaryRows(d).theta_min = min(thetaVec);
    summaryRows(d).theta_max = max(thetaVec);
    summaryRows(d).base_ACC_mean = mean(baseACCVec);
    summaryRows(d).model_ACC_mean = mean(accVec);
    summaryRows(d).dACC_mean = mean(dACCVec);
    if runs > 1
        summaryRows(d).dACC_sd = std(dACCVec, 0);
    else
        summaryRows(d).dACC_sd = 0;
    end
    summaryRows(d).ACC_wins = sum(dACCVec > TOL);
    summaryRows(d).ACC_ties = sum(abs(dACCVec) <= TOL);
    summaryRows(d).ACC_losses = sum(dACCVec < -TOL);
    summaryRows(d).base_NMI_mean = mean(baseNMIVec);
    summaryRows(d).model_NMI_mean = mean(nmiVec);
    summaryRows(d).dNMI_mean = mean(dNMIVec);
    if runs > 1
        summaryRows(d).dNMI_sd = std(dNMIVec, 0);
    else
        summaryRows(d).dNMI_sd = 0;
    end
    summaryRows(d).NMI_wins = sum(dNMIVec > TOL);
    summaryRows(d).NMI_ties = sum(abs(dNMIVec) <= TOL);
    summaryRows(d).NMI_losses = sum(dNMIVec < -TOL);
    summaryRows(d).dual_wins = sum((dACCVec > TOL) & (dNMIVec > TOL));
    summaryRows(d).cv_delta_mean = mean(deltaVec);
    summaryRows(d).cv_se_mean = mean(seVec);
    save(matPath, 'rawRows', 'summaryRows');
end

gcgoc_write_summary_csv(summaryPath, summaryRows);

cvFiniteFlag = 1;
passFlag = 1;
for d = 1:length(summaryRows)
    cvValues = [summaryRows(d).theta_mean, summaryRows(d).theta_min, ...
        summaryRows(d).theta_max, summaryRows(d).cv_delta_mean, ...
        summaryRows(d).cv_se_mean];
    if any(~isfinite(cvValues))
        cvFiniteFlag = 0;
        passFlag = 0;
    end
    if summaryRows(d).dACC_mean < -TOL || summaryRows(d).dNMI_mean < -TOL || ...
            summaryRows(d).ACC_losses > 0 || summaryRows(d).NMI_losses > 0
        passFlag = 0;
    end
end

gcgoc_write_decision(decisionPath, runs, summaryRows, passFlag, cvFiniteFlag);
save(matPath, 'rawRows', 'summaryRows', 'passFlag', 'cvFiniteFlag');

fprintf('\nFinished.\n');
fprintf('Raw results    : %s\n', rawPath);
fprintf('Summary results: %s\n', summaryPath);
fprintf('Decision report: %s\n', decisionPath);
fprintf('MAT results    : %s\n', matPath);
if runs == 3
    fprintf('THREE_SEED_INTEGRITY_PASS=%d\n', passFlag);
else
    fprintf('FIRST_THREE_20SEED_PASS=%d\n', passFlag);
end
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
