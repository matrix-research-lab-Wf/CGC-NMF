function regenerate_release_statistics_R2019a(repoRoot)
%REGENERATE_RELEASE_STATISTICS_R2019A Rebuild archived inferential statistics.
% Compatible with MATLAB R2019a and base MATLAB (no Statistics Toolbox).
%
% Usage:
%   regenerate_release_statistics_R2019a(pwd)
%
% Outputs:
%   results/cgc_gnmfld_transfer/CGC_GNMFLD_exact_wilcoxon_holm.csv
%   results/main_comparison/CGC_minus_GOCNMF_bootstrap_10000.csv

if nargin < 1 || isempty(repoRoot)
    repoRoot = pwd;
end

transferFile = fullfile(repoRoot,'results','cgc_gnmfld_transfer', ...
    'CGC_GNMFLD_transfer_raw_FIXED.csv');
mainFile = fullfile(repoRoot,'results','main_comparison', ...
    'three_direct_methods_20seed_raw.csv');

assert(exist(transferFile,'file') == 2,'Missing transfer raw CSV.');
assert(exist(mainFile,'file') == 2,'Missing main-comparison raw CSV.');

datasets = {'PIE','YaleB','COIL20','COIL100','Optdigits','MNIST'};
metrics = {'ACC','NMI'};

%% Exact paired Wilcoxon tests and baseline-wise Holm correction.
T = readtable(transferFile,'Delimiter',',');
nTests = numel(datasets) * numel(metrics);
outDataset = cell(nTests,1);
outMetric = cell(nTests,1);
nPairs = zeros(nTests,1);
wPlus = zeros(nTests,1);
pRaw = zeros(nTests,1);
wins = zeros(nTests,1);
ties = zeros(nTests,1);
losses = zeros(nTests,1);

row = 0;
for d = 1:numel(datasets)
    mask = strcmp(T.Dataset,datasets{d});
    for m = 1:numel(metrics)
        row = row + 1;
        if strcmp(metrics{m},'ACC')
            diff = T.ACC_CGC_GNMFLD(mask) - T.ACC_GNMFLD(mask);
        else
            diff = T.NMI_CGC_GNMFLD(mask) - T.NMI_GNMFLD(mask);
        end
        [pRaw(row),wPlus(row),nPairs(row)] = exact_signed_rank_two_sided(diff);
        tol = 1e-12;
        wins(row) = sum(diff > tol);
        ties(row) = sum(abs(diff) <= tol);
        losses(row) = sum(diff < -tol);
        outDataset{row} = datasets{d};
        outMetric{row} = metrics{m};
    end
end
pHolm = holm_adjust(pRaw);
wilcoxonTable = table(outDataset,outMetric,nPairs,wPlus,pRaw,pHolm, ...
    wins,ties,losses,'VariableNames',{'dataset','metric','n_nonzero', ...
    'W_plus','p_exact_two_sided','p_holm_12','wins','ties','losses'});
wilcoxonOut = fullfile(repoRoot,'results','cgc_gnmfld_transfer', ...
    'CGC_GNMFLD_exact_wilcoxon_holm.csv');
writetable(wilcoxonTable,wilcoxonOut);

%% Fixed-seed paired percentile bootstrap for CGC-GOCNMF minus GOCNMF.
M = readtable(mainFile,'Delimiter',',');
rng(20261009,'twister');
B = 10000;
nRows = numel(datasets) * numel(metrics);
bDataset = cell(nRows,1);
bMetric = cell(nRows,1);
bPairs = zeros(nRows,1);
bMean = zeros(nRows,1);
bLow = zeros(nRows,1);
bHigh = zeros(nRows,1);
bWins = zeros(nRows,1);
bTies = zeros(nRows,1);
bLosses = zeros(nRows,1);

row = 0;
for d = 1:numel(datasets)
    datasetMask = strcmp(M.dataset,datasets{d});
    G = M(datasetMask & strcmp(M.method,'GOCNMF'),:);
    C = M(datasetMask & strcmp(M.method,'CGC-GOCNMF'),:);
    [commonSeeds,ig,ic] = intersect(G.seed,C.seed,'stable');
    assert(numel(commonSeeds) == 20,'Expected 20 paired seeds for %s.',datasets{d});
    for m = 1:numel(metrics)
        row = row + 1;
        if strcmp(metrics{m},'ACC')
            diff = 100 * (C.unlabeled_ACC(ic) - G.unlabeled_ACC(ig));
        else
            diff = 100 * (C.unlabeled_NMI(ic) - G.unlabeled_NMI(ig));
        end
        n = numel(diff);
        bootMeans = zeros(B,1);
        for b = 1:B
            idx = randi(n,n,1);
            bootMeans(b) = mean(diff(idx));
        end
        ci = percentile_linear(bootMeans,[2.5 97.5]);
        tol = 1e-12;
        bDataset{row} = datasets{d};
        bMetric{row} = metrics{m};
        bPairs(row) = n;
        bMean(row) = mean(diff);
        bLow(row) = ci(1);
        bHigh(row) = ci(2);
        bWins(row) = sum(diff > tol);
        bTies(row) = sum(abs(diff) <= tol);
        bLosses(row) = sum(diff < -tol);
    end
end
bootstrapTable = table(bDataset,bMetric,bPairs,bMean,bLow,bHigh, ...
    bWins,bTies,bLosses,repmat(B,nRows,1),repmat(20261009,nRows,1), ...
    'VariableNames',{'dataset','metric','n_pairs','mean_difference_pp', ...
    'bootstrap95_low_pp','bootstrap95_high_pp','wins','ties','losses', ...
    'bootstrap_resamples','rng_seed'});
bootstrapOut = fullfile(repoRoot,'results','main_comparison', ...
    'CGC_minus_GOCNMF_bootstrap_10000.csv');
writetable(bootstrapTable,bootstrapOut);

fprintf('Wrote: %s\n',wilcoxonOut);
fprintf('Wrote: %s\n',bootstrapOut);
end

function [p,wPlus,n] = exact_signed_rank_two_sided(diff)
tol = 1e-12;
diff = diff(:);
diff = diff(abs(diff) > tol);
n = numel(diff);
if n == 0
    p = 1;
    wPlus = 0;
    return;
end
ranks = local_tied_rank(abs(diff));
weights = round(2*ranks); % integer half-rank units
wPlus2 = sum(weights(diff > 0));
total = sum(weights);
counts = zeros(1,total+1);
counts(1) = 1;
for i = 1:n
    w = weights(i);
    old = counts;
    counts((w+1):end) = counts((w+1):end) + old(1:(end-w));
end
lower = sum(counts(1:(wPlus2+1)));
upper = sum(counts((wPlus2+1):end));
p = min(1,2*min(lower,upper)/(2^n));
wPlus = wPlus2/2;
end

function ranks = local_tied_rank(x)
[sx,order] = sort(x(:));
n = numel(sx);
sortedRanks = zeros(n,1);
i = 1;
while i <= n
    j = i;
    while j < n && sx(j+1) == sx(i)
        j = j + 1;
    end
    sortedRanks(i:j) = (i+j)/2;
    i = j + 1;
end
ranks = zeros(n,1);
ranks(order) = sortedRanks;
end

function adjusted = holm_adjust(p)
p = p(:);
m = numel(p);
[sorted,index] = sort(p);
adjSorted = zeros(m,1);
running = 0;
for i = 1:m
    running = max(running,min(1,(m-i+1)*sorted(i)));
    adjSorted(i) = running;
end
adjusted = zeros(m,1);
adjusted(index) = adjSorted;
end

function q = percentile_linear(x,percent)
x = sort(x(:));
n = numel(x);
q = zeros(size(percent));
for i = 1:numel(percent)
    pos = 1 + (n-1)*percent(i)/100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q(i) = x(lo);
    else
        q(i) = x(lo) + (pos-lo)*(x(hi)-x(lo));
    end
end
end
