function run_GCGOCNMF_graph_stress_R2009a(runs, dataDir, outputDir)
%RUN_GCGOCNMF_GRAPH_STRESS_R2009A
% Controlled graph-reliability stress test for CGC-GOCNMF.
%
% The clean local graph W1 is corrupted without changing its support or its
% marginal edge-weight distribution. Let Wperm be obtained by randomly
% permuting the undirected nonzero weights of W1 over the same edge set:
%
%   Wrho = (1-rho) W1 + rho Wperm,
%
% where rho = 0,0.2,...,1. The graph remains symmetric, nonnegative, and
% supported on exactly the same p-nearest-neighbor edges. Only the association
% between edge location and local similarity is progressively destroyed.
%
% For every rho, Wrho is evaluated against the original binary graph W0 by
% the frozen five-fold out-of-fold Brier rule. The resulting theta calibrates
%
%   Wshrink = (1-theta) W0 + theta Wrho.
%
% Usage:
%   run_GCGOCNMF_graph_stress_R2009a(3,  pwd, pwd);
%   run_GCGOCNMF_graph_stress_R2009a(20, pwd, pwd);
%
% MATLAB R2009a compatible. No K-means or external classifier is used.

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
TOL = 1e-12;
LABEL_FRACTION = 0.10;
ITERATIONS = 50;
FOLDS = 5;
SEED_START = 20260617;
BLOCK_SIZE = 256;
RHO = [0 0.2 0.4 0.6 0.8 1.0];

DATASET_NAME = 'COIL100';
DATA_FILES = {'COIL100_Obj(1).mat', 'COIL100_Obj.mat', 'COIL100.mat'};
P = 3;
ALPHA = 10;
EXPECTED_SAMPLES = 7200;
EXPECTED_CLASSES = 100;

rawPath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_raw.csv', runs));
summaryPath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_summary.csv', runs));
decisionPath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_decision.txt', runs));
matPath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_results.mat', runs));

thetaFigurePath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_theta.png', runs));
cvFigurePath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_cv.png', runs));
accFigurePath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_dACC.png', runs));
nmiFigurePath = fullfile(outputDir, sprintf( ...
    'GCGOCNMF_graph_stress_%dseed_dNMI.png', runs));

fid = fopen(rawPath, 'w');
if fid < 0, error('Cannot create %s.', rawPath); end
fprintf(fid, ['dataset,seed,rho,samples,features,classes,p,alpha,labeled_count,' ...
    'edge_correlation,edge_distortion,theta,delta_cv,se_cv,' ...
    'binary_cv_brier,candidate_cv_brier,' ...
    'base_ACC,candidate_ACC,shrink_ACC,candidate_dACC,shrink_dACC,' ...
    'base_NMI,candidate_NMI,shrink_NMI,candidate_dNMI,shrink_dNMI,' ...
    'cv_seconds,candidate_seconds,shrink_seconds\n']);
fclose(fid);

fprintf('\nCONTROLLED GRAPH-RELIABILITY STRESS TEST -- 2026-07-31\n');
fprintf('Dataset: %s\n', DATASET_NAME);
fprintf('Runs: %d\n', runs);
fprintf('rho levels: 0, 0.2, 0.4, 0.6, 0.8, 1.0\n\n');

dataPath = gcgoc_find_file(dataDir, DATA_FILES);
[X, y] = gcgoc_load_dataset( ...
    dataPath, EXPECTED_SAMPLES, EXPECTED_CLASSES);
n = size(X,2);
m = size(X,1);
c = length(unique(y));

fprintf('Constructing clean binary and local graphs.\n');
[W0, W1] = gcgoc_construct_graphs(X, P, BLOCK_SIZE, EPSILON);
[edgeI, edgeJ, cleanWeight] = find(triu(W1, 1));
edgeCount = length(cleanWeight);
if edgeCount == 0
    error('The local graph has no undirected edge.');
end

numberOfLevels = length(RHO);
thetaValue = zeros(runs, numberOfLevels);
deltaValue = zeros(runs, numberOfLevels);
seValue = zeros(runs, numberOfLevels);
binaryLossValue = zeros(runs, numberOfLevels);
candidateLossValue = zeros(runs, numberOfLevels);
edgeCorrelation = zeros(runs, numberOfLevels);
edgeDistortion = zeros(runs, numberOfLevels);

baseACC = zeros(runs,1);
baseNMI = zeros(runs,1);
candidateACC = zeros(runs, numberOfLevels);
candidateNMI = zeros(runs, numberOfLevels);
shrinkACC = zeros(runs, numberOfLevels);
shrinkNMI = zeros(runs, numberOfLevels);

cvSeconds = zeros(runs, numberOfLevels);
candidateSeconds = zeros(runs, numberOfLevels);
shrinkSeconds = zeros(runs, numberOfLevels);

for r = 1:runs
    seed = SEED_START + r - 1;
    fprintf('\n%s seed %d (%d/%d)\n', DATASET_NAME, seed, r, runs);

    L = gcgoc_labeled_indices(y, seed, LABEL_FRACTION);
    folds = gcgoc_stratified_folds(L, y, seed + 991, FOLDS);

    % The same baseline folds and baseline Brier losses are reused at all rho.
    loss0 = gcgoc_fold_brier_losses(W0, folds, L, y, c, EPSILON);

    baseClock = tic;
    V0 = gcgoc_train(X, y, L, W0, ALPHA, seed, ITERATIONS, EPSILON);
    [baseACC(r), baseNMI(r)] = gcgoc_evaluate(V0, y, L, EPSILON);
    baseElapsed = toc(baseClock); %#ok<NASGU>

    % One fixed random permutation per run creates a nested corruption path.
    gcgoc_set_seed(seed + 700001);
    permutation = randperm(edgeCount);
    if all(permutation(:) == (1:edgeCount)')
        permutation = [2:edgeCount 1];
    end
    permutedWeight = cleanWeight(permutation);
    Wperm = sparse( ...
        [edgeI; edgeJ], [edgeJ; edgeI], ...
        [permutedWeight; permutedWeight], n, n);

    for q = 1:numberOfLevels
        rho = RHO(q);
        currentWeight = (1-rho) * cleanWeight + rho * permutedWeight;
        Wrho = (1-rho) * W1 + rho * Wperm;
        Wrho = sparse(Wrho);

        Ccorr = corrcoef(cleanWeight, currentWeight);
        if numel(Ccorr) < 4 || ~isfinite(Ccorr(1,2))
            if rho == 0
                edgeCorrelation(r,q) = 1;
            else
                edgeCorrelation(r,q) = 0;
            end
        else
            edgeCorrelation(r,q) = Ccorr(1,2);
        end
        edgeDistortion(r,q) = norm(currentWeight-cleanWeight) / ...
            max(norm(cleanWeight), EPSILON);

        cvClock = tic;
        lossRho = gcgoc_fold_brier_losses( ...
            Wrho, folds, L, y, c, EPSILON);
        difference = loss0 - lossRho;
        deltaCV = mean(difference);
        if length(difference) > 1
            seCV = std(difference, 0) / sqrt(length(difference));
        else
            seCV = 0;
        end
        if deltaCV > 0
            theta = max(0, 1-seCV/max(deltaCV,EPSILON));
        else
            theta = 0;
        end
        theta = min(1, theta);
        cvSeconds(r,q) = toc(cvClock);

        Wshrink = (1-theta)*W0 + theta*Wrho;

        candidateClock = tic;
        VC = gcgoc_train( ...
            X, y, L, Wrho, ALPHA, seed, ITERATIONS, EPSILON);
        candidateSeconds(r,q) = toc(candidateClock);
        [candidateACC(r,q), candidateNMI(r,q)] = ...
            gcgoc_evaluate(VC, y, L, EPSILON);

        if theta <= TOL
            shrinkACC(r,q) = baseACC(r);
            shrinkNMI(r,q) = baseNMI(r);
            shrinkSeconds(r,q) = 0;
        elseif abs(theta-1) <= TOL
            shrinkACC(r,q) = candidateACC(r,q);
            shrinkNMI(r,q) = candidateNMI(r,q);
            shrinkSeconds(r,q) = 0;
        else
            shrinkClock = tic;
            VS = gcgoc_train( ...
                X, y, L, Wshrink, ALPHA, seed, ITERATIONS, EPSILON);
            shrinkSeconds(r,q) = toc(shrinkClock);
            [shrinkACC(r,q), shrinkNMI(r,q)] = ...
                gcgoc_evaluate(VS, y, L, EPSILON);
        end

        thetaValue(r,q) = theta;
        deltaValue(r,q) = deltaCV;
        seValue(r,q) = seCV;
        binaryLossValue(r,q) = mean(loss0);
        candidateLossValue(r,q) = mean(lossRho);

        fid = fopen(rawPath, 'a');
        if fid < 0, error('Cannot append to %s.', rawPath); end
        fprintf(fid, ['%s,%d,%.4f,%d,%d,%d,%d,%.15g,%d,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,' ...
            '%.15g,%.15g,%.15g,%.15g,%.15g,' ...
            '%.6f,%.6f,%.6f\n'], ...
            DATASET_NAME, seed, rho, n, m, c, P, ALPHA, length(L), ...
            edgeCorrelation(r,q), edgeDistortion(r,q), ...
            theta, deltaCV, seCV, mean(loss0), mean(lossRho), ...
            baseACC(r), candidateACC(r,q), shrinkACC(r,q), ...
            candidateACC(r,q)-baseACC(r), shrinkACC(r,q)-baseACC(r), ...
            baseNMI(r), candidateNMI(r,q), shrinkNMI(r,q), ...
            candidateNMI(r,q)-baseNMI(r), shrinkNMI(r,q)-baseNMI(r), ...
            cvSeconds(r,q), candidateSeconds(r,q), shrinkSeconds(r,q));
        fclose(fid);

        fprintf(['rho=%.1f corr=%+.3f distortion=%.3f theta=%.3f ' ...
            'dACC candidate/shrink=%+.4f/%+.4f\n'], ...
            rho, edgeCorrelation(r,q), edgeDistortion(r,q), theta, ...
            candidateACC(r,q)-baseACC(r), shrinkACC(r,q)-baseACC(r));
    end
end

summaryTemplate = struct( ...
    'rho',0, 'runs',0, ...
    'edge_correlation_mean',0, 'edge_distortion_mean',0, ...
    'theta_mean',0, 'theta_sd',0, ...
    'delta_cv_mean',0, 'delta_cv_sd',0, 'se_cv_mean',0, ...
    'binary_cv_brier_mean',0, 'candidate_cv_brier_mean',0, ...
    'base_ACC_mean',0, 'candidate_ACC_mean',0, 'shrink_ACC_mean',0, ...
    'candidate_dACC_mean',0, 'shrink_dACC_mean',0, ...
    'candidate_ACC_losses',0, 'shrink_ACC_losses',0, ...
    'base_NMI_mean',0, 'candidate_NMI_mean',0, 'shrink_NMI_mean',0, ...
    'candidate_dNMI_mean',0, 'shrink_dNMI_mean',0, ...
    'candidate_NMI_losses',0, 'shrink_NMI_losses',0, ...
    'cv_seconds_mean',0, 'candidate_seconds_mean',0, ...
    'shrink_seconds_mean',0);
summaryRows = repmat(summaryTemplate, 1, numberOfLevels);

for q = 1:numberOfLevels
    S = summaryTemplate;
    S.rho = RHO(q);
    S.runs = runs;
    S.edge_correlation_mean = mean(edgeCorrelation(:,q));
    S.edge_distortion_mean = mean(edgeDistortion(:,q));
    S.theta_mean = mean(thetaValue(:,q));
    S.theta_sd = std(thetaValue(:,q),0);
    S.delta_cv_mean = mean(deltaValue(:,q));
    S.delta_cv_sd = std(deltaValue(:,q),0);
    S.se_cv_mean = mean(seValue(:,q));
    S.binary_cv_brier_mean = mean(binaryLossValue(:,q));
    S.candidate_cv_brier_mean = mean(candidateLossValue(:,q));

    S.base_ACC_mean = mean(baseACC);
    S.candidate_ACC_mean = mean(candidateACC(:,q));
    S.shrink_ACC_mean = mean(shrinkACC(:,q));
    S.candidate_dACC_mean = mean(candidateACC(:,q)-baseACC);
    S.shrink_dACC_mean = mean(shrinkACC(:,q)-baseACC);
    S.candidate_ACC_losses = sum(candidateACC(:,q)-baseACC < -TOL);
    S.shrink_ACC_losses = sum(shrinkACC(:,q)-baseACC < -TOL);

    S.base_NMI_mean = mean(baseNMI);
    S.candidate_NMI_mean = mean(candidateNMI(:,q));
    S.shrink_NMI_mean = mean(shrinkNMI(:,q));
    S.candidate_dNMI_mean = mean(candidateNMI(:,q)-baseNMI);
    S.shrink_dNMI_mean = mean(shrinkNMI(:,q)-baseNMI);
    S.candidate_NMI_losses = sum(candidateNMI(:,q)-baseNMI < -TOL);
    S.shrink_NMI_losses = sum(shrinkNMI(:,q)-baseNMI < -TOL);

    S.cv_seconds_mean = mean(cvSeconds(:,q));
    S.candidate_seconds_mean = mean(candidateSeconds(:,q));
    S.shrink_seconds_mean = mean(shrinkSeconds(:,q));
    summaryRows(q) = S;
end

fid = fopen(summaryPath, 'w');
if fid < 0, error('Cannot create %s.', summaryPath); end
fprintf(fid, ['rho,runs,edge_correlation_mean,edge_distortion_mean,' ...
    'theta_mean,theta_sd,delta_cv_mean,delta_cv_sd,se_cv_mean,' ...
    'binary_cv_brier_mean,candidate_cv_brier_mean,' ...
    'base_ACC_mean,candidate_ACC_mean,shrink_ACC_mean,' ...
    'candidate_dACC_mean,shrink_dACC_mean,' ...
    'candidate_ACC_losses,shrink_ACC_losses,' ...
    'base_NMI_mean,candidate_NMI_mean,shrink_NMI_mean,' ...
    'candidate_dNMI_mean,shrink_dNMI_mean,' ...
    'candidate_NMI_losses,shrink_NMI_losses,' ...
    'cv_seconds_mean,candidate_seconds_mean,shrink_seconds_mean\n']);
for q = 1:numberOfLevels
    S = summaryRows(q);
    fprintf(fid, ['%.4f,%d,%.15g,%.15g,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d,' ...
        '%.6f,%.6f,%.6f\n'], ...
        S.rho,S.runs,S.edge_correlation_mean,S.edge_distortion_mean, ...
        S.theta_mean,S.theta_sd,S.delta_cv_mean,S.delta_cv_sd,S.se_cv_mean, ...
        S.binary_cv_brier_mean,S.candidate_cv_brier_mean, ...
        S.base_ACC_mean,S.candidate_ACC_mean,S.shrink_ACC_mean, ...
        S.candidate_dACC_mean,S.shrink_dACC_mean, ...
        S.candidate_ACC_losses,S.shrink_ACC_losses, ...
        S.base_NMI_mean,S.candidate_NMI_mean,S.shrink_NMI_mean, ...
        S.candidate_dNMI_mean,S.shrink_dNMI_mean, ...
        S.candidate_NMI_losses,S.shrink_NMI_losses, ...
        S.cv_seconds_mean,S.candidate_seconds_mean,S.shrink_seconds_mean);
end
fclose(fid);

finitePass = 1;
allValues = [thetaValue(:); deltaValue(:); seValue(:); ...
    binaryLossValue(:); candidateLossValue(:); ...
    edgeCorrelation(:); edgeDistortion(:); ...
    baseACC(:); candidateACC(:); shrinkACC(:); ...
    baseNMI(:); candidateNMI(:); shrinkNMI(:)];
if any(~isfinite(allValues))
    finitePass = 0;
end

constructionPass = 1;
if abs(summaryRows(1).edge_correlation_mean-1) > 1e-10
    constructionPass = 0;
end
if abs(summaryRows(1).edge_distortion_mean) > 1e-10
    constructionPass = 0;
end
for q = 2:numberOfLevels
    if summaryRows(q).edge_distortion_mean + TOL < ...
            summaryRows(q-1).edge_distortion_mean
        constructionPass = 0;
    end
end

thetaDrop = summaryRows(numberOfLevels).theta_mean < ...
    summaryRows(1).theta_mean - TOL;
deltaDrop = summaryRows(numberOfLevels).delta_cv_mean < ...
    summaryRows(1).delta_cv_mean - TOL;
candidateACCDegrade = ...
    summaryRows(numberOfLevels).candidate_dACC_mean < ...
    summaryRows(1).candidate_dACC_mean - TOL;
candidateNMIDegrade = ...
    summaryRows(numberOfLevels).candidate_dNMI_mean < ...
    summaryRows(1).candidate_dNMI_mean - TOL;

trendResponsePass = finitePass && constructionPass && ...
    thetaDrop && deltaDrop && candidateACCDegrade && candidateNMIDegrade;

high = summaryRows(numberOfLevels);
shrinkNoMoreLosses = ...
    high.shrink_ACC_losses <= high.candidate_ACC_losses && ...
    high.shrink_NMI_losses <= high.candidate_NMI_losses;
shrinkCloserToBase = ...
    abs(high.shrink_dACC_mean) <= abs(high.candidate_dACC_mean) + TOL && ...
    abs(high.shrink_dNMI_mean) <= abs(high.candidate_dNMI_mean) + TOL;
riskControlPass = finitePass && shrinkNoMoreLosses && shrinkCloserToBase;
stressMechanismPass = trendResponsePass && riskControlPass;

fid = fopen(decisionPath, 'w');
if fid < 0, error('Cannot create %s.', decisionPath); end
fprintf(fid, 'Controlled graph-reliability stress test\n');
fprintf(fid, 'Dataset: %s\nRuns: %d\n', DATASET_NAME, runs);
fprintf(fid, ['Corruption preserves graph support, symmetry, nonnegativity, ' ...
    'and the edge-weight distribution.\n\n']);
for q = 1:numberOfLevels
    S = summaryRows(q);
    fprintf(fid, ['rho=%.1f: edge_corr=%+.6f, distortion=%.6f, ' ...
        'theta=%.6f, delta=%.6f, se=%.6f, ' ...
        'dACC candidate/shrink=%+.6f/%+.6f, losses=%d/%d, ' ...
        'dNMI candidate/shrink=%+.6f/%+.6f, losses=%d/%d\n'], ...
        S.rho,S.edge_correlation_mean,S.edge_distortion_mean, ...
        S.theta_mean,S.delta_cv_mean,S.se_cv_mean, ...
        S.candidate_dACC_mean,S.shrink_dACC_mean, ...
        S.candidate_ACC_losses,S.shrink_ACC_losses, ...
        S.candidate_dNMI_mean,S.shrink_dNMI_mean, ...
        S.candidate_NMI_losses,S.shrink_NMI_losses);
end
fprintf(fid, '\nSTRESS_FINITE_PASS=%d\n', finitePass);
fprintf(fid, 'CORRUPTION_CONSTRUCTION_PASS=%d\n', constructionPass);
fprintf(fid, 'TREND_RESPONSE_PASS=%d\n', trendResponsePass);
fprintf(fid, 'RISK_CONTROL_PASS=%d\n', riskControlPass);
fprintf(fid, 'STRESS_MECHANISM_PASS=%d\n', stressMechanismPass);
if runs == 3
    fprintf(fid, 'THREE_SEED_INTEGRITY_PASS=%d\n', ...
        finitePass && constructionPass);
else
    fprintf(fid, 'TWENTY_SEED_STRESS_PASS=%d\n', ...
        stressMechanismPass);
end
fclose(fid);

save(matPath, 'RHO', 'summaryRows', ...
    'thetaValue','deltaValue','seValue', ...
    'edgeCorrelation','edgeDistortion', ...
    'baseACC','candidateACC','shrinkACC', ...
    'baseNMI','candidateNMI','shrinkNMI', ...
    'finitePass','constructionPass','trendResponsePass', ...
    'riskControlPass','stressMechanismPass');

gcgoc_plot_stress(RHO, thetaValue, ...
    deltaValue, seValue, ...
    candidateACC-repmat(baseACC,1,numberOfLevels), ...
    shrinkACC-repmat(baseACC,1,numberOfLevels), ...
    candidateNMI-repmat(baseNMI,1,numberOfLevels), ...
    shrinkNMI-repmat(baseNMI,1,numberOfLevels), ...
    thetaFigurePath,cvFigurePath,accFigurePath,nmiFigurePath);

fprintf('\nFinished.\n');
fprintf('Raw: %s\nSummary: %s\nDecision: %s\n', ...
    rawPath, summaryPath, decisionPath);
fprintf('STRESS_FINITE_PASS=%d\n', finitePass);
fprintf('CORRUPTION_CONSTRUCTION_PASS=%d\n', constructionPass);
fprintf('TREND_RESPONSE_PASS=%d\n', trendResponsePass);
fprintf('RISK_CONTROL_PASS=%d\n', riskControlPass);
fprintf('STRESS_MECHANISM_PASS=%d\n', stressMechanismPass);
if runs == 3
    fprintf('THREE_SEED_INTEGRITY_PASS=%d\n', ...
        finitePass && constructionPass);
else
    fprintf('TWENTY_SEED_STRESS_PASS=%d\n', stressMechanismPass);
end
end


function gcgoc_plot_stress(rho, thetaValue, deltaValue, seValue, ...
    candidateDACC, shrinkDACC, candidateDNMI, shrinkDNMI, ...
    thetaPath, cvPath, accPath, nmiPath)

runs = size(thetaValue,1);
denominator = sqrt(max(runs,1));

thetaMean = mean(thetaValue,1);
thetaSE = std(thetaValue,0,1) / denominator;
deltaMean = mean(deltaValue,1);
seMean = mean(seValue,1);
candidateACCMean = mean(candidateDACC,1) * 100;
candidateACCSE = std(candidateDACC,0,1) / denominator * 100;
shrinkACCMean = mean(shrinkDACC,1) * 100;
shrinkACCSE = std(shrinkDACC,0,1) / denominator * 100;
candidateNMIMean = mean(candidateDNMI,1) * 100;
candidateNMISE = std(candidateDNMI,0,1) / denominator * 100;
shrinkNMIMean = mean(shrinkDNMI,1) * 100;
shrinkNMISE = std(shrinkDNMI,0,1) / denominator * 100;

f = figure('Visible','off');
errorbar(rho,thetaMean,thetaSE,'ko-','LineWidth',1.5,'MarkerSize',6);
xlabel('Corruption level \rho');
ylabel('Calibration coefficient \theta');
xlim([0 1]); ylim([0 1]);
grid on; set(gca,'FontSize',10);
set(f,'PaperPositionMode','auto');
print(f,'-dpng','-r300',thetaPath); close(f);

f = figure('Visible','off');
plot(rho,deltaMean,'b-o','LineWidth',1.5,'MarkerSize',6); hold on;
plot(rho,seMean,'r-s','LineWidth',1.5,'MarkerSize',6);
plot(rho,zeros(size(rho)),'k--','LineWidth',1);
xlabel('Corruption level \rho');
ylabel('Out-of-fold Brier difference');
legend('\Delta','s','zero','Location','Best');
xlim([0 1]); grid on; set(gca,'FontSize',10);
set(f,'PaperPositionMode','auto');
print(f,'-dpng','-r300',cvPath); close(f);

f = figure('Visible','off');
errorbar(rho,candidateACCMean,candidateACCSE, ...
    'b-o','LineWidth',1.5,'MarkerSize',6); hold on;
errorbar(rho,shrinkACCMean,shrinkACCSE, ...
    'r-s','LineWidth',1.5,'MarkerSize',6);
plot(rho,zeros(size(rho)),'k--','LineWidth',1);
xlabel('Corruption level \rho');
ylabel('\DeltaACC (percentage points)');
legend('Candidate graph','SHRINK','BASE','Location','Best');
xlim([0 1]); grid on; set(gca,'FontSize',10);
set(f,'PaperPositionMode','auto');
print(f,'-dpng','-r300',accPath); close(f);

f = figure('Visible','off');
errorbar(rho,candidateNMIMean,candidateNMISE, ...
    'b-o','LineWidth',1.5,'MarkerSize',6); hold on;
errorbar(rho,shrinkNMIMean,shrinkNMISE, ...
    'r-s','LineWidth',1.5,'MarkerSize',6);
plot(rho,zeros(size(rho)),'k--','LineWidth',1);
xlabel('Corruption level \rho');
ylabel('\DeltaNMI (percentage points)');
legend('Candidate graph','SHRINK','BASE','Location','Best');
xlim([0 1]); grid on; set(gca,'FontSize',10);
set(f,'PaperPositionMode','auto');
print(f,'-dpng','-r300',nmiPath); close(f);
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
