function run_three_methods_end_to_end_timing_R2009a(runs,dataDir,outputDir)
%RUN_THREE_METHODS_END_TO_END_TIMING_R2009A
% Dedicated end-to-end timing for GNMFLD, GOCNMF, and CGC-GOCNMF.
%
% Timed cost:
%   T_total = T_graph + T_calibration + T_factorization
%
% Data loading and common column normalization are excluded because they are
% identical preprocessing inputs for all methods. Every method independently
% constructs its own graph inside its timed block. CGC-GOCNMF additionally
% includes the complete five-fold out-of-fold Brier calibration. Even when
% theta=0, CGC-GOCNMF performs an independent factorization and never reuses
% the GOCNMF output.
%
% Usage:
%   run_three_methods_end_to_end_timing_R2009a(1,pwd,pwd);
%   run_three_methods_end_to_end_timing_R2009a(3,pwd,pwd);
%   run_three_methods_end_to_end_timing_R2009a(20,pwd,pwd);
%
% MATLAB R2009a compatible. Checkpointing occurs after every method run.

if nargin < 1 || isempty(runs), runs = 1; end
if nargin < 2 || isempty(dataDir), dataDir = pwd; end
if nargin < 3 || isempty(outputDir), outputDir = pwd; end
if ~(runs==1 || runs==3 || runs==20)
    error('runs must be 1, 3, or 20.');
end
if exist(outputDir,'dir')~=7, mkdir(outputDir); end

EPSILON = 1e-12;
TOL = 1e-12;
LABEL_FRACTION = 0.10;
GOC_ITERATIONS = 50;
GNMFLD_ITERATIONS = 200;
FOLDS = 5;
SEED_START = 20260617;
BLOCK_SIZE = 256;
GNMFLD_P = 5;

D(1).name='PIE';
D(1).files={'CMU_PIE_fac.mat','CMU_PIE.mat','PIE.mat','PIE_fac.mat'};
D(1).expectedSamples=2856; D(1).expectedClasses=68;
D(1).gocP=3; D(1).gocAlpha=1000;
D(1).gnAlpha=1e5; D(1).gnBeta=10;

D(2).name='YaleB';
D(2).files={'YaleB.mat','YaleB(1).mat','YaleB_32x32.mat'};
D(2).expectedSamples=2414; D(2).expectedClasses=38;
D(2).gocP=2; D(2).gocAlpha=1000;
D(2).gnAlpha=1e5; D(2).gnBeta=100;

D(3).name='COIL20';
D(3).files={'COIL20_Obj.mat','COIL20.mat','COIL20_Obj(1).mat'};
D(3).expectedSamples=1440; D(3).expectedClasses=20;
D(3).gocP=3; D(3).gocAlpha=10;
D(3).gnAlpha=1e4; D(3).gnBeta=10;

D(4).name='COIL100';
D(4).files={'COIL100_Obj.mat','COIL100.mat','COIL100_Obj(1).mat'};
D(4).expectedSamples=7200; D(4).expectedClasses=100;
D(4).gocP=3; D(4).gocAlpha=10;
D(4).gnAlpha=1e4; D(4).gnBeta=10;

D(5).name='Optdigits';
D(5).files={'Optdigits_Han.mat','Optdigits.mat','optdigits.mat'};
D(5).expectedSamples=5620; D(5).expectedClasses=10;
D(5).gocP=4; D(5).gocAlpha=10;
D(5).gnAlpha=1e5; D(5).gnBeta=10;

D(6).name='MNIST';
D(6).files={'MNIST_Han.mat','MNIST.mat','mnist.mat'};
D(6).expectedSamples=6996; D(6).expectedClasses=10;
D(6).gocP=4; D(6).gocAlpha=10;
D(6).gnAlpha=1e5; D(6).gnBeta=10;

methodNames={'GNMFLD','GOCNMF','CGC-GOCNMF'};
numberOfDatasets=length(D);
numberOfMethods=length(methodNames);

prefix=sprintf('three_methods_end_to_end_timing_%dseed',runs);
rawPath=fullfile(outputDir,[prefix '_raw.csv']);
summaryPath=fullfile(outputDir,[prefix '_summary.csv']);
decisionPath=fullfile(outputDir,[prefix '_decision.txt']);
latexPath=fullfile(outputDir,[prefix '_table.tex']);
checkpointPath=fullfile(outputDir,[prefix '_checkpoint.mat']);

rowTemplate=struct( ...
    'dataset','', 'seed',0, 'method','', ...
    'graph_seconds',0, 'calibration_seconds',0, ...
    'factorization_seconds',0, 'total_seconds',0, ...
    'wall_seconds',0, 'component_gap',0, ...
    'theta',NaN, 'delta_cv',NaN, 'se_cv',NaN, ...
    'unlabeled_ACC',0, 'unlabeled_NMI',0, ...
    'objective_nonincrease',0, 'execution_order',0);

totalRows=numberOfDatasets*runs*numberOfMethods;
rows=repmat(rowTemplate,1,totalRows);
completed=false(numberOfDatasets,runs,numberOfMethods);
warmDone=false(numberOfDatasets,numberOfMethods);

if exist(checkpointPath,'file')==2
    P=load(checkpointPath);
    if isfield(P,'rows') && isfield(P,'completed') && ...
            length(P.rows)==totalRows && ...
            all(size(P.completed)==[numberOfDatasets runs numberOfMethods])
        rows=P.rows;
        completed=P.completed;
        if isfield(P,'warmDone'), warmDone=P.warmDone; end
        fprintf('Resuming checkpoint: %s\n',checkpointPath);
    else
        error('Existing checkpoint is incompatible with this run.');
    end
end

fprintf('\nDEDICATED END-TO-END TIMING -- 2026-07-31\n');
fprintf('Runs per dataset: %d\n',runs);
fprintf('Data loading/common normalization excluded.\n');
fprintf('Every method builds its own graph; no graph or factor output is shared.\n\n');

for d=1:numberOfDatasets
    dataPath=tdm_find_file_recursive(dataDir,D(d).files);
    [X,y]=tdm_load_dataset(dataPath,D(d).expectedSamples,D(d).expectedClasses);
    [m,n]=size(X);
    c=length(unique(y));

    % Full untimed warm-up once per method and dataset.
    for k=1:numberOfMethods
        if ~warmDone(d,k)
            seed=SEED_START;
            L=tdm_labeled_indices(y,seed,LABEL_FRACTION);
            [U0,V0]=tdm_initial_factors(m,n,c,seed,EPSILON);
            fprintf('Warm-up: %s %s\n',D(d).name,methodNames{k});
            tdm_time_one_method(k,X,y,L,U0,V0,D(d),GNMFLD_P, ...
                GNMFLD_ITERATIONS,GOC_ITERATIONS,FOLDS,seed, ...
                BLOCK_SIZE,EPSILON);
            warmDone(d,k)=true;
            save(checkpointPath,'rows','completed','warmDone','D','methodNames');
        end
    end

    for r=1:runs
        seed=SEED_START+r-1;
        L=tdm_labeled_indices(y,seed,LABEL_FRACTION);
        [U0,V0]=tdm_initial_factors(m,n,c,seed,EPSILON);

        % Reproducible random method order controls fixed-order/cache bias.
        tdm_set_seed(seed+900001+1000*d);
        order=randperm(numberOfMethods);

        fprintf('\n%s seed %d, order: %s -> %s -> %s\n', ...
            D(d).name,seed,methodNames{order(1)}, ...
            methodNames{order(2)},methodNames{order(3)});

        for position=1:numberOfMethods
            k=order(position);
            if completed(d,r,k)
                fprintf('  %s already complete; skipping.\n',methodNames{k});
                continue;
            end

            R=tdm_time_one_method(k,X,y,L,U0,V0,D(d),GNMFLD_P, ...
                GNMFLD_ITERATIONS,GOC_ITERATIONS,FOLDS,seed, ...
                BLOCK_SIZE,EPSILON);

            index=((d-1)*runs+(r-1))*numberOfMethods+k;
            row=rowTemplate;
            row.dataset=D(d).name;
            row.seed=seed;
            row.method=methodNames{k};
            row.graph_seconds=R.graphSeconds;
            row.calibration_seconds=R.calibrationSeconds;
            row.factorization_seconds=R.factorizationSeconds;
            row.total_seconds=R.totalSeconds;
            row.wall_seconds=R.wallSeconds;
            row.component_gap=R.componentGap;
            row.theta=R.theta;
            row.delta_cv=R.deltaCV;
            row.se_cv=R.seCV;
            row.unlabeled_ACC=R.metrics.unlabeledACC;
            row.unlabeled_NMI=R.metrics.unlabeledNMI;
            row.objective_nonincrease=R.objectiveNonincrease;
            row.execution_order=position;
            rows(index)=row;
            completed(d,r,k)=true;

            save(checkpointPath,'rows','completed','warmDone','D','methodNames');
            tdm_time_write_raw(rawPath,rows,completed,runs,numberOfMethods);

            fprintf(['  %s: graph=%.3f, calibration=%.3f, ' ...
                'factorization=%.3f, total=%.3f s, ACC=%.4f, NMI=%.4f\n'], ...
                methodNames{k},R.graphSeconds,R.calibrationSeconds, ...
                R.factorizationSeconds,R.totalSeconds, ...
                R.metrics.unlabeledACC,R.metrics.unlabeledNMI);
        end
    end
end

if ~all(completed(:))
    error('Timing experiment incomplete.');
end

summaryRows=tdm_time_summarize(rows,D,methodNames,runs);
tdm_time_write_summary(summaryPath,summaryRows);
tdm_time_write_latex(latexPath,summaryRows,D,methodNames);
tdm_time_write_decision(decisionPath,runs,rows,summaryRows,TOL);
save(checkpointPath,'rows','completed','warmDone','summaryRows','D','methodNames');

fprintf('\nFinished.\n');
fprintf('Raw: %s\n',rawPath);
fprintf('Summary: %s\n',summaryPath);
fprintf('LaTeX table: %s\n',latexPath);
fprintf('Decision: %s\n',decisionPath);
if runs==1
    fprintf('ONE_SEED_TIMING_SMOKE_PASS=1 (confirm in decision file)\n');
elseif runs==3
    fprintf('THREE_SEED_TIMING_INTEGRITY_PASS=1 (confirm in decision file)\n');
else
    fprintf('TWENTY_SEED_END_TO_END_TIMING_COMPLETE=1 (confirm in decision file)\n');
end
end


function R=tdm_time_one_method(k,X,y,L,U0,V0,D,gnP, ...
    gnIterations,gocIterations,folds,seed,blockSize,epsilon)

n=size(X,2);
c=length(unique(y));
R.theta=NaN; R.deltaCV=NaN; R.seCV=NaN;
R.calibrationSeconds=0;

wallClock=tic;

graphClock=tic;
if k==1
    [neighborIndex,distance]=tdm_time_knn_quiet(X,gnP,blockSize);
    W=tdm_binary_graph(neighborIndex,n);
elseif k==2
    [neighborIndex,distance]=tdm_time_knn_quiet(X,D.gocP,blockSize);
    W=tdm_binary_graph(neighborIndex,n);
else
    [neighborIndex,distance]=tdm_time_knn_quiet(X,D.gocP,blockSize);
    W0=tdm_binary_graph(neighborIndex,n);
    W1=tdm_local_graph(neighborIndex,distance,n,epsilon);
end
R.graphSeconds=toc(graphClock);

if k==3
    calibrationClock=tic;
    [R.theta,R.deltaCV,R.seCV]=tdm_shrinkage_weight( ...
        W0,W1,L,y,c,seed,folds,epsilon);
    W=(1-R.theta)*W0+R.theta*W1;
    R.calibrationSeconds=toc(calibrationClock);
end

factorClock=tic;
if k==1
    [V,obj0,obj1,monotone]=tdm_train_gnmfld( ...
        X,y,L,W,D.gnAlpha,D.gnBeta,U0,V0,gnIterations,epsilon);
else
    [V,obj0,obj1,monotone]=tdm_train_goc( ...
        X,y,L,W,D.gocAlpha,U0,V0,gocIterations,epsilon);
end
R.factorizationSeconds=toc(factorClock);

R.wallSeconds=toc(wallClock);
R.totalSeconds=R.graphSeconds+R.calibrationSeconds+ ...
    R.factorizationSeconds;
R.componentGap=R.wallSeconds-R.totalSeconds;
R.objectiveNonincrease=monotone;
R.metrics=tdm_evaluate(V,y,L,epsilon);

if any(~isfinite([R.graphSeconds R.calibrationSeconds ...
        R.factorizationSeconds R.totalSeconds R.wallSeconds ...
        R.metrics.unlabeledACC R.metrics.unlabeledNMI obj0 obj1]))
    error('Non-finite timing or metric result.');
end
end


function tdm_time_write_raw(path,rows,completed,runs,numberOfMethods)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,seed,method,graph_seconds,calibration_seconds,' ...
    'factorization_seconds,total_seconds,wall_seconds,component_gap,' ...
    'theta,delta_cv,se_cv,unlabeled_ACC,unlabeled_NMI,' ...
    'objective_nonincrease,execution_order\n']);
numberOfDatasets=size(completed,1);
for d=1:numberOfDatasets
    for r=1:runs
        for k=1:numberOfMethods
            if completed(d,r,k)
                index=((d-1)*runs+(r-1))*numberOfMethods+k;
                S=rows(index);
                fprintf(fid,['%s,%d,%s,%.9f,%.9f,%.9f,%.9f,%.9f,' ...
                    '%.9f,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%d\n'], ...
                    S.dataset,S.seed,S.method,S.graph_seconds, ...
                    S.calibration_seconds,S.factorization_seconds, ...
                    S.total_seconds,S.wall_seconds,S.component_gap, ...
                    S.theta,S.delta_cv,S.se_cv,S.unlabeled_ACC, ...
                    S.unlabeled_NMI,S.objective_nonincrease, ...
                    S.execution_order);
            end
        end
    end
end
fclose(fid);
end


function summary=tdm_time_summarize(rows,D,methodNames,runs)
template=struct('dataset','', 'method','', 'runs',0, ...
    'graph_mean',0,'graph_sd',0, ...
    'calibration_mean',0,'calibration_sd',0, ...
    'factorization_mean',0,'factorization_sd',0, ...
    'total_mean',0,'total_sd',0, ...
    'wall_mean',0,'wall_sd',0, ...
    'ACC_mean',0,'NMI_mean',0, ...
    'theta_mean',NaN,'objective_pass',0);
summary=repmat(template,1,length(D)*length(methodNames));
counter=0;
for d=1:length(D)
    for k=1:length(methodNames)
        counter=counter+1;
        graph=zeros(runs,1); calibration=zeros(runs,1);
        factorization=zeros(runs,1); total=zeros(runs,1);
        wall=zeros(runs,1); acc=zeros(runs,1); nmi=zeros(runs,1);
        theta=zeros(runs,1); objective=zeros(runs,1);
        for r=1:runs
            index=((d-1)*runs+(r-1))*length(methodNames)+k;
            S=rows(index);
            graph(r)=S.graph_seconds;
            calibration(r)=S.calibration_seconds;
            factorization(r)=S.factorization_seconds;
            total(r)=S.total_seconds;
            wall(r)=S.wall_seconds;
            acc(r)=S.unlabeled_ACC;
            nmi(r)=S.unlabeled_NMI;
            theta(r)=S.theta;
            objective(r)=S.objective_nonincrease;
        end
        Q=template;
        Q.dataset=D(d).name; Q.method=methodNames{k}; Q.runs=runs;
        Q.graph_mean=mean(graph); Q.graph_sd=tdm_time_std(graph);
        Q.calibration_mean=mean(calibration);
        Q.calibration_sd=tdm_time_std(calibration);
        Q.factorization_mean=mean(factorization);
        Q.factorization_sd=tdm_time_std(factorization);
        Q.total_mean=mean(total); Q.total_sd=tdm_time_std(total);
        Q.wall_mean=mean(wall); Q.wall_sd=tdm_time_std(wall);
        Q.ACC_mean=mean(acc); Q.NMI_mean=mean(nmi);
        if k==3, Q.theta_mean=mean(theta); end
        Q.objective_pass=all(objective==1);
        summary(counter)=Q;
    end
end
end


function value=tdm_time_std(x)
if length(x)>1, value=std(x,0); else value=0; end
end


function tdm_time_write_summary(path,rows)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,method,runs,graph_mean,graph_sd,' ...
    'calibration_mean,calibration_sd,factorization_mean,' ...
    'factorization_sd,total_mean,total_sd,wall_mean,wall_sd,' ...
    'ACC_mean,NMI_mean,theta_mean,objective_pass\n']);
for i=1:length(rows)
    S=rows(i);
    fprintf(fid,['%s,%s,%d,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,' ...
        '%.9f,%.9f,%.9f,%.9f,%.15g,%.15g,%.15g,%d\n'], ...
        S.dataset,S.method,S.runs,S.graph_mean,S.graph_sd, ...
        S.calibration_mean,S.calibration_sd,S.factorization_mean, ...
        S.factorization_sd,S.total_mean,S.total_sd,S.wall_mean, ...
        S.wall_sd,S.ACC_mean,S.NMI_mean,S.theta_mean,S.objective_pass);
end
fclose(fid);
end


function tdm_time_write_latex(path,rows,D,methodNames)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'%% Requires booktabs.\n');
fprintf(fid,'\\begin{table*}[t]\n\\centering\n');
fprintf(fid,['\\caption{End-to-end running time in seconds. Each value is ' ...
    'the mean $\\pm$ standard deviation over paired trials and includes ' ...
    'method-specific graph construction, calibration when applicable, ' ...
    'and matrix factorization. Common data loading and normalization are ' ...
    'excluded.}\\label{tab:end_to_end_time}\n']);
fprintf(fid,'\\setlength{\\tabcolsep}{6pt}\n');
fprintf(fid,'\\begin{tabular}{lccc}\n\\toprule\n');
fprintf(fid,'Dataset & GNMFLD & GOCNMF & CGC-GOCNMF \\\\\n\\midrule\n');
for d=1:length(D)
    values=zeros(3,2);
    for k=1:length(methodNames)
        for i=1:length(rows)
            if strcmp(rows(i).dataset,D(d).name) && ...
                    strcmp(rows(i).method,methodNames{k})
                values(k,1)=rows(i).total_mean;
                values(k,2)=rows(i).total_sd;
            end
        end
    end
    fprintf(fid,'%s & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ \\\\\n', ...
        D(d).name,values(1,1),values(1,2),values(2,1),values(2,2), ...
        values(3,1),values(3,2));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table*}\n');
fclose(fid);
end


function tdm_time_write_decision(path,runs,rows,summary,tol)
numeric=[];
gapRatio=zeros(length(rows),1);
for i=1:length(rows)
    numeric=[numeric;rows(i).graph_seconds; ...
        rows(i).calibration_seconds;rows(i).factorization_seconds; ...
        rows(i).total_seconds;rows(i).wall_seconds; ...
        rows(i).unlabeled_ACC;rows(i).unlabeled_NMI]; %#ok<AGROW>
    gapRatio(i)=abs(rows(i).component_gap)/max(rows(i).wall_seconds,tol);
end
finitePass=all(isfinite(numeric));
positivePass=all(numeric(1:7:end)>=0) && ...
    all(numeric(2:7:end)>=0) && all(numeric(3:7:end)>=0) && ...
    all(numeric(4:7:end)>0);
componentPass=max(gapRatio)<0.02;
objectivePass=all([summary.objective_pass]==1);
cgcCalibrationPass=1;
for i=1:length(rows)
    if strcmp(rows(i).method,'CGC-GOCNMF')
        if rows(i).calibration_seconds<=0 || ~isfinite(rows(i).theta)
            cgcCalibrationPass=0;
        end
    elseif abs(rows(i).calibration_seconds)>tol
        cgcCalibrationPass=0;
    end
end
integrityPass=finitePass && positivePass && componentPass && ...
    objectivePass && cgcCalibrationPass;

fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'Dedicated end-to-end timing audit\n');
fprintf(fid,'Runs per dataset: %d\n',runs);
fprintf(fid,['Included: independent graph construction, complete CGC ' ...
    'calibration, independent factorization.\n']);
fprintf(fid,['Excluded: common data loading, common column normalization, ' ...
    'and post-training evaluation.\n']);
fprintf(fid,'Method execution order was randomized reproducibly per seed.\n');
fprintf(fid,'CGC factorization was never reused from GOCNMF.\n\n');

for i=1:length(summary)
    S=summary(i);
    fprintf(fid,['%s %s: graph=%.4f+-%.4f, calibration=%.4f+-%.4f, ' ...
        'factorization=%.4f+-%.4f, total=%.4f+-%.4f s\n'], ...
        S.dataset,S.method,S.graph_mean,S.graph_sd, ...
        S.calibration_mean,S.calibration_sd,S.factorization_mean, ...
        S.factorization_sd,S.total_mean,S.total_sd);
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
fprintf(fid,'COMPONENT_WALL_AGREEMENT_PASS=%d\n',componentPass);
fprintf(fid,'MAX_COMPONENT_WALL_RELATIVE_GAP=%.9g\n',max(gapRatio));
fprintf(fid,'OBJECTIVE_NONINCREASE_PASS=%d\n',objectivePass);
fprintf(fid,'CGC_CALIBRATION_INCLUDED_PASS=%d\n',cgcCalibrationPass);
fprintf(fid,'END_TO_END_TIMING_INTEGRITY_PASS=%d\n',integrityPass);
if runs==1
    fprintf(fid,'ONE_SEED_TIMING_SMOKE_PASS=%d\n',integrityPass);
elseif runs==3
    fprintf(fid,'THREE_SEED_TIMING_INTEGRITY_PASS=%d\n',integrityPass);
else
    fprintf(fid,'TWENTY_SEED_END_TO_END_TIMING_COMPLETE=%d\n',integrityPass);
end
fclose(fid);
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

function [neighborIndex,neighborDistance] = tdm_time_knn_quiet(X,p,blockSize)
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
    clear similarities sortedSimilarity sortedIndex;
end
end

