function analyze_NEU_statistics_R2009a(rawCsv, outputDir)
%ANALYZE_NEU_STATISTICS_R2009A
% Reproduce the formal paired statistics for the NEU-CLS application.
%
% Requires MATLAB Statistics Toolbox for SIGNRANK and TIEDRANK.
% Compatible with MATLAB R2009a.
%
% The automatic Wilcoxon rule used here is:
%   - exact calculation when there are no zero differences and no tied
%     absolute differences;
%   - normal approximation otherwise.
% This matches the exact/approximate distinction used in the paper audit.
%
% Holm correction is applied separately for each baseline over the nine
% combinations of three label rates and three metrics.

if nargin < 1 || isempty(rawCsv)
    rawCsv = fullfile(pwd,'NEU_application_20seed_raw.csv');
end
if nargin < 2 || isempty(outputDir)
    outputDir = pwd;
end
if exist(rawCsv,'file') ~= 2
    error('Raw CSV not found: %s',rawCsv);
end
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

fid = fopen(rawCsv,'r');
if fid < 0
    error('Cannot open %s.',rawCsv);
end
header = fgetl(fid); %#ok<NASGU>
% Columns:
% fraction, seed, method, labeled_count, theta, delta, se, fallback,
% ACC, NMI, ARI, labeled_ACC, obj0, obj1, monotone, time
C = textscan(fid,['%f%f%s' repmat('%f',1,13)], ...
    'Delimiter',',','CollectOutput',false);
fclose(fid);

fraction = C{1};
seed = C{2};
method = C{3};
numeric = zeros(length(fraction),13);
for j = 1:13
    numeric(:,j) = C{j+3};
end
ACC = numeric(:,6);
NMI = numeric(:,7);
ARI = numeric(:,8);

rates = [0.05 0.10 0.20];
baselineNames = {'GOCNMF','GNMFLD'};
metricNames = {'ACC','NMI','ARI'};
metricValues = {ACC,NMI,ARI};
numberOfTests = 9;

outputPath = fullfile(outputDir, ...
    'NEU_application_paired_statistics_MATLAB.csv');
fid = fopen(outputPath,'w');
fprintf(fid,['baseline,label_fraction,metric,n_pairs,mean_difference,' ...
    'median_difference,bootstrap95_low,bootstrap95_high,' ...
    'wilcoxon_W,p_raw,p_holm,rank_biserial,wins,ties,losses\n']);

rand('twister',20260801); %#ok<RAND>

for b = 1:length(baselineNames)
    localRows = cell(numberOfTests,1);
    rawP = zeros(numberOfTests,1);
    testIndex = 0;

    for q = 1:length(rates)
        rate = rates(q);
        cgcIndex = abs(fraction-rate)<1e-12 & ...
            strcmp(method,'CGC-GOCNMF');
        baseIndex = abs(fraction-rate)<1e-12 & ...
            strcmp(method,baselineNames{b});

        cgcSeed = seed(cgcIndex);
        baseSeed = seed(baseIndex);
        [commonSeed,ia,ib] = intersect(cgcSeed,baseSeed); %#ok<ASGLU>

        if length(commonSeed) ~= 20
            error('Expected 20 paired seeds.');
        end

        for k = 1:length(metricNames)
            testIndex = testIndex+1;
            cgcMetric = metricValues{k}(cgcIndex);
            baseMetric = metricValues{k}(baseIndex);
            d = cgcMetric(ia)-baseMetric(ib);

            nonzero = abs(d)>1e-12;
            dnz = d(nonzero);
            absValues = abs(dnz);
            tied = length(unique(absValues)) < length(absValues);

            if isempty(dnz)
                p = NaN;
                W = NaN;
            elseif ~tied && all(nonzero)
                [p,dummy,stats] = signrank(dnz,0, ...
                    'method','exact'); %#ok<ASGLU>
                W = min(stats.signedrank, ...
                    sum(tiedrank(abs(dnz)))-stats.signedrank);
            else
                [p,dummy,stats] = signrank(dnz,0, ...
                    'method','approximate'); %#ok<ASGLU>
                ranks = tiedrank(abs(dnz));
                Wplus = sum(ranks(dnz>0));
                Wminus = sum(ranks(dnz<0));
                W = min(Wplus,Wminus);
            end

            repetitions = 100000;
            bootMean = zeros(repetitions,1);
            n = length(d);
            for r = 1:repetitions
                ids = ceil(n*rand(n,1)); %#ok<RAND>
                bootMean(r) = mean(d(ids));
            end
            sortedBoot = sort(bootMean);
            lowIndex = max(1,round(0.025*repetitions));
            highIndex = min(repetitions,round(0.975*repetitions));
            ciLow = sortedBoot(lowIndex);
            ciHigh = sortedBoot(highIndex);

            ranks = tiedrank(abs(dnz));
            Wplus = sum(ranks(dnz>0));
            Wminus = sum(ranks(dnz<0));
            if Wplus+Wminus==0
                effect = NaN;
            else
                effect = (Wplus-Wminus)/(Wplus+Wminus);
            end

            row = struct;
            row.baseline = baselineNames{b};
            row.rate = rate;
            row.metric = metricNames{k};
            row.n = length(d);
            row.meanDifference = mean(d);
            row.medianDifference = median(d);
            row.ciLow = ciLow;
            row.ciHigh = ciHigh;
            row.W = W;
            row.pRaw = p;
            row.effect = effect;
            row.wins = sum(d>1e-12);
            row.ties = sum(abs(d)<=1e-12);
            row.losses = sum(d<-1e-12);
            localRows{testIndex} = row;
            rawP(testIndex) = p;
        end
    end

    adjustedP = local_holm(rawP);

    for t = 1:numberOfTests
        row = localRows{t};
        fprintf(fid,['%s,%.6f,%s,%d,%.15g,%.15g,%.15g,%.15g,' ...
            '%.15g,%.15g,%.15g,%.15g,%d,%d,%d\n'], ...
            row.baseline,row.rate,row.metric,row.n, ...
            row.meanDifference,row.medianDifference,row.ciLow,row.ciHigh, ...
            row.W,row.pRaw,adjustedP(t),row.effect, ...
            row.wins,row.ties,row.losses);
    end
end

fclose(fid);
fprintf('Output: %s\n',outputPath);
fprintf('NEU_FORMAL_STATISTICS_COMPLETE=1\n');
end


function adjusted = local_holm(p)
m = length(p);
[sortedP,order] = sort(p);
sortedAdjusted = zeros(m,1);
running = 0;
for i = 1:m
    value = (m-i+1)*sortedP(i);
    running = max(running,value);
    sortedAdjusted(i) = min(1,running);
end
adjusted = zeros(m,1);
adjusted(order) = sortedAdjusted;
end
