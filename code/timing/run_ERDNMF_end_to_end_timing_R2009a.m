function run_ERDNMF_end_to_end_timing_R2009a(runs,dataDir,outputDir)
%RUN_ERDNMF_END_TO_END_TIMING_R2009A
% Dedicated timing supplement for adding ERDNMF to Table 10.
%
% Included:
%   ERDNMF matrix factorization with the locked 200 iterations.
%
% Excluded:
%   data loading, common column normalization, and post-training evaluation.
%
% ERDNMF has no neighborhood-graph construction or calibration stage, so
%   T_total = T_factorization.
%
% The program uses the same datasets, labeled fraction, seeds, labeled
% subsets, initialization rule, and final ERDNMF coefficients as the unified
% six-dataset experiment.
%
% Usage:
%   run_ERDNMF_end_to_end_timing_R2009a(1,pwd,pwd);
%   run_ERDNMF_end_to_end_timing_R2009a(3,pwd,pwd);
%   run_ERDNMF_end_to_end_timing_R2009a(20,pwd,pwd);
%
% MATLAB R2009a compatible. Checkpointing occurs after every dataset-seed run.

if nargin < 1 || isempty(runs), runs = 1; end
if nargin < 2 || isempty(dataDir), dataDir = pwd; end
if nargin < 3 || isempty(outputDir), outputDir = pwd; end
if ~(runs == 1 || runs == 3 || runs == 20)
    error('runs must be 1, 3, or 20.');
end
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

EPSILON = 1e-12;
LABEL_FRACTION = 0.10;
ITERATIONS = 200;
SEED_START = 20260617;

D(1).name = 'PIE';
D(1).files = {'CMU_PIE_fac.mat','CMU_PIE.mat','PIE.mat','PIE_fac.mat'};
D(1).expectedSamples = 2856;
D(1).expectedClasses = 68;
D(1).alpha = 400;
D(1).source = 'original-authors-code';

D(2).name = 'YaleB';
D(2).files = {'YaleB.mat','YaleB(1).mat','YaleB_32x32.mat'};
D(2).expectedSamples = 2414;
D(2).expectedClasses = 38;
D(2).alpha = 1000;
D(2).source = 'labeled-only-CV';

D(3).name = 'COIL20';
D(3).files = {'COIL20_Obj.mat','COIL20.mat','COIL20_Obj(1).mat'};
D(3).expectedSamples = 1440;
D(3).expectedClasses = 20;
D(3).alpha = 800;
D(3).source = 'original-authors-code';

D(4).name = 'COIL100';
D(4).files = {'COIL100_Obj.mat','COIL100.mat','COIL100_Obj(1).mat'};
D(4).expectedSamples = 7200;
D(4).expectedClasses = 100;
D(4).alpha = 1000;
D(4).source = 'original-authors-code';

D(5).name = 'Optdigits';
D(5).files = {'Optdigits_Han.mat','Optdigits.mat','optdigits.mat'};
D(5).expectedSamples = 5620;
D(5).expectedClasses = 10;
D(5).alpha = 600;
D(5).source = 'original-authors-code';

D(6).name = 'MNIST';
D(6).files = {'MNIST_Han.mat','MNIST.mat','mnist.mat'};
D(6).expectedSamples = 6996;
D(6).expectedClasses = 10;
D(6).alpha = 1000;
D(6).source = 'labeled-only-CV';

numberOfDatasets = length(D);
prefix = sprintf('ERDNMF_end_to_end_timing_%dseed',runs);
rawPath = fullfile(outputDir,[prefix '_raw.csv']);
summaryPath = fullfile(outputDir,[prefix '_summary.csv']);
decisionPath = fullfile(outputDir,[prefix '_decision.txt']);
latexPath = fullfile(outputDir,[prefix '_table_column.tex']);
checkpointPath = fullfile(outputDir,[prefix '_checkpoint.mat']);

template = struct( ...
    'dataset','', ...
    'seed',0, ...
    'alpha',0, ...
    'parameter_source','', ...
    'factorization_seconds',0, ...
    'total_seconds',0, ...
    'unlabeled_ACC',0, ...
    'unlabeled_NMI',0, ...
    'objective_initial',0, ...
    'objective_final',0, ...
    'objective_nonincrease',0);

rows = repmat(template,1,numberOfDatasets*runs);
completed = false(numberOfDatasets,runs);
warmDone = false(numberOfDatasets,1);

if exist(checkpointPath,'file') == 2
    P = load(checkpointPath);
    if isfield(P,'rows') && isfield(P,'completed') && ...
            length(P.rows) == numberOfDatasets*runs && ...
            all(size(P.completed) == [numberOfDatasets runs])
        rows = P.rows;
        completed = P.completed;
        if isfield(P,'warmDone')
            warmDone = P.warmDone;
        end
        fprintf('Resuming checkpoint: %s\n',checkpointPath);
    else
        error('Existing checkpoint is incompatible with this run count.');
    end
end

fprintf('\nERDNMF DEDICATED TABLE-10 TIMING\n');
fprintf('Runs per dataset: %d\n',runs);
fprintf('Seeds: %d to %d\n',SEED_START,SEED_START+runs-1);
fprintf('Iterations: %d\n',ITERATIONS);
fprintf('Data loading, common normalization, and evaluation are excluded.\n');
fprintf('ERDNMF has no graph or calibration stage.\n\n');

for d = 1:numberOfDatasets
    dataPath = erdt_find_file_recursive(dataDir,D(d).files);
    [X,y] = erdt_load_dataset( ...
        dataPath,D(d).expectedSamples,D(d).expectedClasses);
    [m,n] = size(X);
    c = length(unique(y));

    if ~warmDone(d)
        seed = SEED_START;
        L = erdt_labeled_indices(y,seed,LABEL_FRACTION);
        [U0,V0] = erdt_initial_factors(m,n,c,seed,EPSILON);
        fprintf('Warm-up: %s\n',D(d).name);
        erdt_train(X,y,L,D(d).alpha,U0,V0,ITERATIONS,EPSILON);
        warmDone(d) = true;
        save(checkpointPath,'rows','completed','warmDone','D','runs');
    end

    for r = 1:runs
        if completed(d,r)
            fprintf('%s seed %d already complete; skipping.\n', ...
                D(d).name,SEED_START+r-1);
            continue;
        end

        seed = SEED_START+r-1;
        L = erdt_labeled_indices(y,seed,LABEL_FRACTION);
        [U0,V0] = erdt_initial_factors(m,n,c,seed,EPSILON);

        methodClock = tic;
        [V,obj0,obj1,nonincrease] = erdt_train( ...
            X,y,L,D(d).alpha,U0,V0,ITERATIONS,EPSILON);
        elapsed = toc(methodClock);

        metrics = erdt_evaluate(V,y,L,EPSILON);

        index = (d-1)*runs+r;
        row = template;
        row.dataset = D(d).name;
        row.seed = seed;
        row.alpha = D(d).alpha;
        row.parameter_source = D(d).source;
        row.factorization_seconds = elapsed;
        row.total_seconds = elapsed;
        row.unlabeled_ACC = metrics.unlabeledACC;
        row.unlabeled_NMI = metrics.unlabeledNMI;
        row.objective_initial = obj0;
        row.objective_final = obj1;
        row.objective_nonincrease = nonincrease;
        rows(index) = row;
        completed(d,r) = true;

        save(checkpointPath,'rows','completed','warmDone','D','runs');
        erdt_write_raw(rawPath,rows,completed,runs);

        fprintf(['%s seed=%d total=%.4fs ACC=%.6f NMI=%.6f ' ...
            'monotone=%d\n'], ...
            D(d).name,seed,elapsed,metrics.unlabeledACC, ...
            metrics.unlabeledNMI,nonincrease);
    end
end

if ~all(completed(:))
    error('Timing experiment is incomplete.');
end

summary = erdt_summarize(rows,D,runs);
erdt_write_raw(rawPath,rows,completed,runs);
erdt_write_summary(summaryPath,summary);
erdt_write_latex(latexPath,summary);
erdt_write_decision(decisionPath,runs,rows,summary);

save(checkpointPath,'rows','completed','warmDone','summary','D','runs');

fprintf('\nFinished.\n');
fprintf('RAW=%s\n',rawPath);
fprintf('SUMMARY=%s\n',summaryPath);
fprintf('LATEX=%s\n',latexPath);
fprintf('DECISION=%s\n',decisionPath);

if runs == 1
    fprintf('ONE_SEED_ERDNMF_TIMING_PASS=1 (confirm in decision file)\n');
elseif runs == 3
    fprintf('THREE_SEED_ERDNMF_TIMING_PASS=1 (confirm in decision file)\n');
else
    fprintf('TWENTY_SEED_ERDNMF_TIMING_COMPLETE=1 (confirm in decision file)\n');
end
end


function path = erdt_find_file_recursive(folder,candidates)
path = '';
for i = 1:length(candidates)
    direct = fullfile(folder,candidates{i});
    if exist(direct,'file') == 2
        path = direct;
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
            path = erdt_find_file_recursive(child,candidates);
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


function [X,y] = erdt_load_dataset(path,expectedSamples,expectedClasses)
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


function erdt_set_seed(seed)
rand('twister',double(seed)); %#ok<RAND>
end


function L = erdt_labeled_indices(y,seed,fraction)
erdt_set_seed(seed);
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


function [U0,V0] = erdt_initial_factors(m,n,c,seed,epsilon)
erdt_set_seed(seed+500001);
U0 = max(rand(m,c),epsilon);
V0 = max(rand(n,c),epsilon);
end


function [V,obj0,obj1,nonincrease] = erdt_train( ...
    X,y,L,alpha,U,V,iterations,epsilon)

U = max(U,epsilon);
V = max(V,epsilon);
[n,c] = size(V);

Y = zeros(n,c);
for t = 1:length(L)
    Y(L(t),y(L(t))) = 1;
end

H = zeros(n,c);
H(L,:) = 1-Y(L,:);

obj0 = erdt_objective(X,U,V,Y,H,alpha,epsilon);
previous = obj0;
nonincrease = 1;

for iter = 1:iterations
    numeratorU = X*V;
    denominatorU = U*(V'*V);
    U = U.*(numeratorU./max(denominatorU,epsilon));
    U = max(U,epsilon);

    sumCorrect = max(sum(sum(Y.*V)),epsilon);
    sumWrong = sum(sum(H.*V));

    numeratorV = X'*U + ...
        alpha*(sumWrong/(sumCorrect^2))*Y;
    denominatorV = V*(U'*U) + ...
        alpha*(1/sumCorrect)*H;

    V = V.*(numeratorV./max(denominatorV,epsilon));
    V = max(V,epsilon);

    current = erdt_objective(X,U,V,Y,H,alpha,epsilon);
    if current > previous + 1e-10*max(1,abs(previous))
        nonincrease = 0;
    end
    previous = current;
end

obj1 = previous;
end


function value = erdt_objective(X,U,V,Y,H,alpha,epsilon)
residual = X-U*V';
ratio = sum(sum(H.*V))/max(sum(sum(Y.*V)),epsilon);
value = sum(sum(residual.^2)) + 2*alpha*ratio;
end


function M = erdt_evaluate(V,y,L,epsilon)
[dummy,prediction] = max(V,[],2); %#ok<ASGLU>
mask = true(length(y),1);
mask(L) = false;
U = find(mask);

M.unlabeledACC = mean(prediction(U) == y(U));
M.unlabeledNMI = erdt_nmi(y(U),prediction(U),epsilon);
end


function value = erdt_nmi(trueLabel,predictedLabel,epsilon)
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
            mi = mi + P(i,j)*log( ...
                P(i,j)/max(pi(i)*pj(j),epsilon));
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


function erdt_write_raw(path,rows,completed,runs)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,['dataset,seed,method,alpha,parameter_source,' ...
    'graph_seconds,calibration_seconds,factorization_seconds,' ...
    'total_seconds,unlabeled_ACC,unlabeled_NMI,' ...
    'objective_initial,objective_final,objective_nonincrease\n']);

for d = 1:size(completed,1)
    for r = 1:runs
        if completed(d,r)
            index = (d-1)*runs+r;
            S = rows(index);
            fprintf(fid,['%s,%d,ERDNMF,%.15g,%s,0,0,' ...
                '%.9f,%.9f,%.15g,%.15g,%.15g,%.15g,%d\n'], ...
                S.dataset,S.seed,S.alpha,S.parameter_source, ...
                S.factorization_seconds,S.total_seconds, ...
                S.unlabeled_ACC,S.unlabeled_NMI, ...
                S.objective_initial,S.objective_final, ...
                S.objective_nonincrease);
        end
    end
end

fclose(fid);
end


function summary = erdt_summarize(rows,D,runs)
template = struct( ...
    'dataset','', ...
    'method','ERDNMF', ...
    'runs',runs, ...
    'alpha',0, ...
    'total_mean',0, ...
    'total_sd',0, ...
    'ACC_mean',0, ...
    'NMI_mean',0, ...
    'objective_pass',0);

summary = repmat(template,1,length(D));

for d = 1:length(D)
    ids = (d-1)*runs+(1:runs);
    T = [rows(ids).total_seconds];
    A = [rows(ids).unlabeled_ACC];
    N = [rows(ids).unlabeled_NMI];
    O = [rows(ids).objective_nonincrease];

    summary(d).dataset = D(d).name;
    summary(d).alpha = D(d).alpha;
    summary(d).total_mean = mean(T);
    summary(d).total_sd = erdt_std(T);
    summary(d).ACC_mean = mean(A);
    summary(d).NMI_mean = mean(N);
    summary(d).objective_pass = all(O == 1);
end
end


function value = erdt_std(x)
if length(x) > 1
    value = std(x,0);
else
    value = 0;
end
end


function erdt_write_summary(path,summary)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,['dataset,method,runs,alpha,total_mean,total_sd,' ...
    'ACC_mean,NMI_mean,objective_pass\n']);

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,'%s,ERDNMF,%d,%.15g,%.9f,%.9f,%.15g,%.15g,%d\n', ...
        S.dataset,S.runs,S.alpha,S.total_mean,S.total_sd, ...
        S.ACC_mean,S.NMI_mean,S.objective_pass);
end

fclose(fid);
end


function erdt_write_latex(path,summary)
fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'%% ERDNMF column for the revised Table 10.\n');
fprintf(fid,'%% Merge these values with the existing three-method timing table.\n');
fprintf(fid,'\\begin{tabular}{lc}\n');
fprintf(fid,'\\hline\nDataset & ERDNMF \\\\\n\\hline\n');

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,'%s & $%.3f\\pm%.3f$ \\\\\n', ...
        S.dataset,S.total_mean,S.total_sd);
end

fprintf(fid,'\\hline\n\\end{tabular}\n');
fclose(fid);
end


function erdt_write_decision(path,runs,rows,summary)
numeric = [];
for i = 1:length(rows)
    numeric = [numeric; ...
        rows(i).factorization_seconds; ...
        rows(i).total_seconds; ...
        rows(i).unlabeled_ACC; ...
        rows(i).unlabeled_NMI; ...
        rows(i).objective_initial; ...
        rows(i).objective_final]; %#ok<AGROW>
end

finitePass = all(isfinite(numeric));
positivePass = all([rows.total_seconds] > 0);
objectivePass = all([summary.objective_pass] == 1);
integrityPass = finitePass && positivePass && objectivePass;

fid = fopen(path,'w');
if fid < 0, error('Cannot create %s.',path); end

fprintf(fid,'ERDNMF dedicated Table-10 timing audit\n');
fprintf(fid,'Runs per dataset: %d\n',runs);
fprintf(fid,'Seeds: 20260617-%d\n',20260617+runs-1);
fprintf(fid,'Included: ERDNMF factorization with 200 iterations.\n');
fprintf(fid,['Excluded: data loading, common column normalization, ' ...
    'and post-training evaluation.\n']);
fprintf(fid,'Graph time=0 and calibration time=0 by model design.\n');
fprintf(fid,'One untimed warm-up was performed on every dataset.\n\n');

for i = 1:length(summary)
    S = summary(i);
    fprintf(fid,'%s ERDNMF: total=%.4f+-%.4f s, alpha=%g\n', ...
        S.dataset,S.total_mean,S.total_sd,S.alpha);
end

fprintf(fid,'\nEnvironment\n');
fprintf(fid,'MATLAB_VERSION=%s\n',version);
fprintf(fid,'COMPUTER=%s\n',computer);
fprintf(fid,'OS=%s\n',getenv('OS'));
fprintf(fid,'PROCESSOR_IDENTIFIER=%s\n',getenv('PROCESSOR_IDENTIFIER'));
fprintf(fid,'NUMBER_OF_PROCESSORS=%s\n',getenv('NUMBER_OF_PROCESSORS'));

fprintf(fid,'\nAudit\n');
fprintf(fid,'FINITE_TIMING_PASS=%d\n',finitePass);
fprintf(fid,'POSITIVE_TIMING_PASS=%d\n',positivePass);
fprintf(fid,'OBJECTIVE_NONINCREASE_PASS=%d\n',objectivePass);
fprintf(fid,'ERDNMF_TIMING_INTEGRITY_PASS=%d\n',integrityPass);

if runs == 1
    fprintf(fid,'ONE_SEED_ERDNMF_TIMING_PASS=%d\n',integrityPass);
elseif runs == 3
    fprintf(fid,'THREE_SEED_ERDNMF_TIMING_PASS=%d\n',integrityPass);
else
    fprintf(fid,'TWENTY_SEED_ERDNMF_TIMING_COMPLETE=%d\n',integrityPass);
end

fclose(fid);
end
