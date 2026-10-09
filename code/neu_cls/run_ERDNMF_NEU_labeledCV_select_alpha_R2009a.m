function run_ERDNMF_NEU_labeledCV_select_alpha_R2009a(runs,dataDir,outputDir)
%RUN_ERDNMF_NEU_LABELEDCV_SELECT_ALPHA_R2009A
% Select one ERDNMF alpha for all NEU-CLS label rates using labeled-only
% five-fold validation. MATLAB R2009a compatible.
%
% Pilot seeds: 20260711--20260730.
% Formal application seeds 20260731--20260750 are excluded.
% Label rates: 5%, 10%, and 20%.
% Candidate alphas: 1, 10, 100, 1000, 10000.
% Selection: highest pooled held-out labeled ACC among candidates passing
% the objective-nonincreasing audit in every pilot run; smaller alpha on
% an exact tie. No unlabeled ground-truth label is used.
%
% Usage:
%   run_ERDNMF_NEU_labeledCV_select_alpha_R2009a(1,pwd,pwd);
%   run_ERDNMF_NEU_labeledCV_select_alpha_R2009a(3,pwd,pwd);
%   run_ERDNMF_NEU_labeledCV_select_alpha_R2009a(20,pwd,pwd);

if nargin<1 || isempty(runs), runs=1; end
if nargin<2 || isempty(dataDir), dataDir=pwd; end
if nargin<3 || isempty(outputDir), outputDir=pwd; end
if ~(runs==1 || runs==3 || runs==20)
    error('runs must be 1, 3, or 20.');
end
if exist(outputDir,'dir')~=7, mkdir(outputDir); end

EPSILON=1e-12;
LABEL_FRACTIONS=[0.05 0.10 0.20];
ALPHAS=[1 10 100 1000 10000];
FOLDS=5;
MAX_ITERATIONS=200;
PILOT_SEED_START=20260711;
TIE_TOL=1e-12;

dataPath=neu_erd_find_file_recursive(dataDir, ...
    {'NEU_CLS_32x32.mat','NEU-CLS-32x32.mat','NEU_CLS.mat'});
[X,y]=neu_erd_load_dataset(dataPath,1800,6);
[m,n]=size(X);
c=length(unique(y));

prefix=sprintf('ERDNMF_NEU_labeledCV_alpha_selection_%dseed',runs);
rawPath=fullfile(outputDir,[prefix '_raw.csv']);
summaryPath=fullfile(outputDir,[prefix '_summary.csv']);
selectedPath=fullfile(outputDir,[prefix '_selected.csv']);
decisionPath=fullfile(outputDir,[prefix '_decision.txt']);
checkpointPath=fullfile(outputDir,[prefix '_checkpoint.mat']);

nr=length(LABEL_FRACTIONS);
na=length(ALPHAS);
totalRows=nr*runs*na*FOLDS;
template=struct('label_fraction',0,'pilot_seed',0,'fold',0,'alpha',0, ...
    'available_labeled_count',0,'training_labeled_count',0, ...
    'validation_labeled_count',0,'validation_ACC',0,'validation_NMI',0, ...
    'training_labeled_ACC',0,'objective_initial',0,'objective_final',0, ...
    'objective_nonincrease',0,'run_seconds',0);
rows=repmat(template,1,totalRows);
completed=false(nr,runs,na,FOLDS);

if exist(checkpointPath,'file')==2
    P=load(checkpointPath);
    if isfield(P,'rows') && isfield(P,'completed') && ...
            length(P.rows)==totalRows && ...
            all(size(P.completed)==[nr runs na FOLDS])
        rows=P.rows;
        completed=P.completed;
        fprintf('Resuming checkpoint: %s\n',checkpointPath);
    else
        error('Existing checkpoint is incompatible with this run count.');
    end
end

fprintf('\nERDNMF NEU-CLS LABELED-ONLY ALPHA SELECTION\n');
fprintf('Data: %s\n',dataPath);
fprintf('Pilot seeds: %d to %d\n',PILOT_SEED_START,PILOT_SEED_START+runs-1);
fprintf('Formal seeds excluded: 20260731 to 20260750\n');
fprintf('Candidate alphas: '); fprintf('%g ',ALPHAS); fprintf('\n');
fprintf('Label rates: 5%%, 10%%, 20%%; folds=%d; iterations=%d\n\n', ...
    FOLDS,MAX_ITERATIONS);

for q=1:nr
    fraction=LABEL_FRACTIONS(q);
    for r=1:runs
        pilotSeed=PILOT_SEED_START+r-1;
        L=neu_erd_labeled_indices(y,pilotSeed,fraction);
        foldSets=neu_erd_stratified_folds(y,L,pilotSeed+200000,FOLDS);
        for a=1:na
            alpha=ALPHAS(a);
            for f=1:FOLDS
                if completed(q,r,a,f), continue; end
                validation=foldSets{f};
                training=setdiff(L,validation);
                if isempty(validation) || isempty(training)
                    error('Empty training/validation split at rate %.2f seed %d fold %d.', ...
                        fraction,pilotSeed,f);
                end
                initSeed=pilotSeed+500001+1000*f+100000*q;
                [U0,V0]=neu_erd_initial_factors(m,n,c,initSeed,EPSILON);
                t0=tic;
                [V,obj0,obj1,nonincrease]=neu_erd_train( ...
                    X,y,training,alpha,U0,V0,MAX_ITERATIONS,EPSILON);
                elapsed=toc(t0);

                idx=neu_erd_cv_index(q,r,a,f,runs,na,FOLDS);
                row=template;
                row.label_fraction=fraction;
                row.pilot_seed=pilotSeed;
                row.fold=f;
                row.alpha=alpha;
                row.available_labeled_count=length(L);
                row.training_labeled_count=length(training);
                row.validation_labeled_count=length(validation);
                row.validation_ACC=neu_erd_accuracy(V,y,validation);
                row.validation_NMI=neu_erd_nmi_on_indices(V,y,validation,EPSILON);
                row.training_labeled_ACC=neu_erd_accuracy(V,y,training);
                row.objective_initial=obj0;
                row.objective_final=obj1;
                row.objective_nonincrease=nonincrease;
                row.run_seconds=elapsed;
                rows(idx)=row;
                completed(q,r,a,f)=true;

                save(checkpointPath,'rows','completed','runs','LABEL_FRACTIONS','ALPHAS');
                neu_erd_write_cv_raw(rawPath,rows,completed,runs,na,FOLDS);
                fprintf(['rate=%.2f seed=%d alpha=%g fold=%d ' ...
                    'ACC=%.4f NMI=%.4f mono=%d time=%.2fs\n'], ...
                    fraction,pilotSeed,alpha,f,row.validation_ACC, ...
                    row.validation_NMI,nonincrease,elapsed);
            end
        end
    end
end

summary=neu_erd_cv_summarize(rows,LABEL_FRACTIONS,ALPHAS,runs,FOLDS);
[selectedAlpha,selectedIndex,ready]=neu_erd_cv_select( ...
    summary,ALPHAS,runs,nr,FOLDS,TIE_TOL);
neu_erd_write_cv_raw(rawPath,rows,completed,runs,na,FOLDS);
neu_erd_write_cv_summary(summaryPath,summary);
neu_erd_write_cv_selected(selectedPath,selectedAlpha,selectedIndex,summary, ...
    runs,PILOT_SEED_START,ready);
neu_erd_write_cv_decision(decisionPath,runs,LABEL_FRACTIONS,ALPHAS,summary, ...
    selectedAlpha,selectedIndex,ready);

if all(completed(:)) && exist(checkpointPath,'file')==2
    delete(checkpointPath);
end
fprintf('\nSELECTED_ALPHA=%g\n',selectedAlpha);
fprintf('FINAL_READY=%d\n',ready && runs==20);
fprintf('RAW=%s\nSUMMARY=%s\nSELECTED=%s\nDECISION=%s\n', ...
    rawPath,summaryPath,selectedPath,decisionPath);
end

function idx=neu_erd_cv_index(q,r,a,f,runs,na,F)
idx=(((q-1)*runs+(r-1))*na+(a-1))*F+f;
end

function path=neu_erd_find_file_recursive(folder,candidates)
path='';
for i=1:length(candidates)
    direct=fullfile(folder,candidates{i});
    if exist(direct,'file')==2, path=direct; return; end
end
listing=dir(folder);
for i=1:length(listing)
    if listing(i).isdir && ~strcmp(listing(i).name,'.') && ...
            ~strcmp(listing(i).name,'..')
        child=fullfile(folder,listing(i).name);
        try
            path=neu_erd_find_file_recursive(child,candidates);
        catch
            path='';
        end
        if ~isempty(path), return; end
    end
end
error('NEU-CLS MAT file was not found under %s.',folder);
end

function [X,y]=neu_erd_load_dataset(path,expectedSamples,expectedClasses)
S=load(path);
if ~isfield(S,'fea') || ~isfield(S,'gnd')
    error('%s must contain fea and gnd.',path);
end
F=double(S.fea);
y0=double(S.gnd(:));
if size(F,1)==length(y0), X=F';
elseif size(F,2)==length(y0), X=F;
else, error('Feature shape does not match labels.'); end
[d1,d2,y]=unique(y0); %#ok<ASGLU>
y=double(y(:));
if size(X,2)~=expectedSamples || length(unique(y))~=expectedClasses
    error('Expected %d samples and %d classes.',expectedSamples,expectedClasses);
end
if min(X(:))<-1e-12 || any(~isfinite(X(:)))
    error('Features must be finite and nonnegative.');
end
norms=sqrt(sum(X.^2,1)); norms(norms==0)=1;
X=bsxfun(@rdivide,X,norms);
end

function neu_erd_set_seed(seed)
rand('twister',double(seed)); %#ok<RAND>
end

function L=neu_erd_labeled_indices(y,seed,fraction)
neu_erd_set_seed(seed);
classes=unique(y);
L=zeros(0,1);
for k=1:length(classes)
    ids=find(y==classes(k));
    order=randperm(length(ids));
    count=max(2,floor(fraction*length(ids)));
    count=min(count,length(ids));
    L=[L; ids(order(1:count))]; %#ok<AGROW>
end
L=double(L(:));
end

function folds=neu_erd_stratified_folds(y,L,seed,F)
neu_erd_set_seed(seed);
folds=cell(F,1);
classes=unique(y(L));
startOffset=0;
for k=1:length(classes)
    ids=L(y(L)==classes(k));
    order=randperm(length(ids));
    ids=ids(order);
    for t=1:length(ids)
        f=mod((t-1)+startOffset,F)+1;
        folds{f}=[folds{f}; ids(t)]; %#ok<AGROW>
    end
    startOffset=mod(startOffset+1,F);
end
for f=1:F
    folds{f}=double(folds{f}(:));
    if isempty(folds{f}), error('Global fold %d is empty.',f); end
end
end

function [U0,V0]=neu_erd_initial_factors(m,n,c,seed,epsilon)
neu_erd_set_seed(seed);
U0=max(rand(m,c),epsilon);
V0=max(rand(n,c),epsilon);
end

function [V,obj0,obj1,nonincrease]=neu_erd_train( ...
    X,y,L,alpha,U,V,maxIterations,epsilon)
[n,c]=size(V);
Y=zeros(n,c);
for t=1:length(L), Y(L(t),y(L(t)))=1; end
H=zeros(n,c);
H(L,:)=1-Y(L,:);
obj0=neu_erd_objective(X,U,V,Y,H,alpha,epsilon);
previous=obj0;
nonincrease=1;
for iter=1:maxIterations
    numeratorU=X*V;
    denominatorU=U*(V'*V);
    U=U.*(numeratorU./max(denominatorU,epsilon));
    U=max(U,epsilon);

    sumCorrect=max(sum(sum(Y.*V)),epsilon);
    sumWrong=sum(sum(H.*V));
    numeratorV=X'*U+alpha*(sumWrong/(sumCorrect^2))*Y;
    denominatorV=V*(U'*U)+alpha*(1/sumCorrect)*H;
    V=V.*(numeratorV./max(denominatorV,epsilon));
    V=max(V,epsilon);

    current=neu_erd_objective(X,U,V,Y,H,alpha,epsilon);
    if current>previous+1e-10*max(1,abs(previous))
        nonincrease=0;
    end
    previous=current;
end
obj1=previous;
end

function value=neu_erd_objective(X,U,V,Y,H,alpha,epsilon)
R=X-U*V';
ratio=sum(sum(H.*V))/max(sum(sum(Y.*V)),epsilon);
value=sum(sum(R.^2))+2*alpha*ratio;
end

function value=neu_erd_accuracy(V,y,ids)
[dummy,pred]=max(V,[],2); %#ok<ASGLU>
value=mean(pred(ids)==y(ids));
end

function value=neu_erd_nmi_on_indices(V,y,ids,epsilon)
[dummy,pred]=max(V,[],2); %#ok<ASGLU>
value=neu_erd_nmi(y(ids),pred(ids),epsilon);
end

function value=neu_erd_nmi(trueLabel,predictedLabel,epsilon)
trueLabel=trueLabel(:);
predictedLabel=predictedLabel(:);
[tv,d1,ti]=unique(trueLabel); %#ok<ASGLU>
[pv,d2,pi]=unique(predictedLabel); %#ok<ASGLU>
C=accumarray([ti pi],1,[length(tv) length(pv)]);
P=C/length(trueLabel);
pt=sum(P,2);
pp=sum(P,1);
mi=0;
for i=1:size(P,1)
    for j=1:size(P,2)
        if P(i,j)>0
            mi=mi+P(i,j)*log(P(i,j)/max(pt(i)*pp(j),epsilon));
        end
    end
end
ht=-sum(pt(pt>0).*log(pt(pt>0)));
hp=-sum(pp(pp>0).*log(pp(pp>0)));
den=0.5*(ht+hp);
if den<=epsilon, value=1; else, value=mi/den; end
value=max(0,min(1,value));
end

function neu_erd_write_cv_raw(path,rows,completed,runs,na,F)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['label_fraction,pilot_seed,fold,alpha,available_labeled_count,' ...
    'training_labeled_count,validation_labeled_count,validation_ACC,' ...
    'validation_NMI,training_labeled_ACC,objective_initial,objective_final,' ...
    'objective_nonincrease,run_seconds\n']);
for q=1:size(completed,1)
    for r=1:runs
        for a=1:na
            for f=1:F
                if completed(q,r,a,f)
                    idx=neu_erd_cv_index(q,r,a,f,runs,na,F);
                    S=rows(idx);
                    fprintf(fid,['%.6f,%d,%d,%.15g,%d,%d,%d,%.15g,%.15g,' ...
                        '%.15g,%.15g,%.15g,%d,%.6f\n'], ...
                        S.label_fraction,S.pilot_seed,S.fold,S.alpha, ...
                        S.available_labeled_count,S.training_labeled_count, ...
                        S.validation_labeled_count,S.validation_ACC, ...
                        S.validation_NMI,S.training_labeled_ACC, ...
                        S.objective_initial,S.objective_final, ...
                        S.objective_nonincrease,S.run_seconds);
                end
            end
        end
    end
end
fclose(fid);
end

function summary=neu_erd_cv_summarize(rows,fractions,alphas,runs,F)
nr=length(fractions);
na=length(alphas);
template=struct('alpha',0,'evaluations',0,'ACC_mean',0,'ACC_sd',0, ...
    'NMI_mean',0,'NMI_sd',0,'train_ACC_mean',0,'mono_pass',0, ...
    'time_mean',0,'ACC_rate05',0,'ACC_rate10',0,'ACC_rate20',0);
summary=repmat(template,1,na);
for a=1:na
    ids=zeros(1,nr*runs*F);
    t=0;
    rateACC=zeros(1,nr);
    for q=1:nr
        rateIds=zeros(1,runs*F);
        u=0;
        for r=1:runs
            for f=1:F
                t=t+1;
                u=u+1;
                idx=neu_erd_cv_index(q,r,a,f,runs,na,F);
                ids(t)=idx;
                rateIds(u)=idx;
            end
        end
        rateACC(q)=mean([rows(rateIds).validation_ACC]);
    end
    A=[rows(ids).validation_ACC];
    N=[rows(ids).validation_NMI];
    T=[rows(ids).training_labeled_ACC];
    O=[rows(ids).objective_nonincrease];
    R=[rows(ids).run_seconds];
    summary(a).alpha=alphas(a);
    summary(a).evaluations=length(ids);
    summary(a).ACC_mean=mean(A);
    summary(a).ACC_sd=neu_erd_safe_std(A);
    summary(a).NMI_mean=mean(N);
    summary(a).NMI_sd=neu_erd_safe_std(N);
    summary(a).train_ACC_mean=mean(T);
    summary(a).mono_pass=sum(O==1);
    summary(a).time_mean=mean(R);
    summary(a).ACC_rate05=rateACC(1);
    summary(a).ACC_rate10=rateACC(2);
    summary(a).ACC_rate20=rateACC(3);
end
end

function [alpha,index,ready]=neu_erd_cv_select(summary,alphas,runs,nr,F,tol)
required=runs*nr*F;
admissible=find([summary.mono_pass]==required);
ready=~isempty(admissible);
if ~ready
    alpha=NaN;
    index=0;
    return;
end
best=admissible(1);
for k=2:length(admissible)
    j=admissible(k);
    if summary(j).ACC_mean>summary(best).ACC_mean+tol
        best=j;
    elseif abs(summary(j).ACC_mean-summary(best).ACC_mean)<=tol && ...
            alphas(j)<alphas(best)
        best=j;
    end
end
alpha=alphas(best);
index=best;
end

function v=neu_erd_safe_std(x)
if length(x)>1, v=std(x,0); else, v=0; end
end

function neu_erd_write_cv_summary(path,summary)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['alpha,evaluations,validation_ACC_mean,validation_ACC_sd,' ...
    'validation_NMI_mean,validation_NMI_sd,training_labeled_ACC_mean,' ...
    'objective_nonincrease_pass,run_seconds_mean,ACC_rate05,' ...
    'ACC_rate10,ACC_rate20\n']);
for i=1:length(summary)
    S=summary(i);
    fprintf(fid,['%.15g,%d,%.15g,%.15g,%.15g,%.15g,%.15g,%d,' ...
        '%.6f,%.15g,%.15g,%.15g\n'], ...
        S.alpha,S.evaluations,S.ACC_mean,S.ACC_sd,S.NMI_mean,S.NMI_sd, ...
        S.train_ACC_mean,S.mono_pass,S.time_mean,S.ACC_rate05, ...
        S.ACC_rate10,S.ACC_rate20);
end
fclose(fid);
end

function neu_erd_write_cv_selected(path,alpha,index,summary,runs,seedStart,ready)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,selected_alpha,validation_ACC_mean,validation_NMI_mean,' ...
    'selection_rule,pilot_seeds,ready\n']);
if ready
    fprintf(fid,['NEU-CLS,%.15g,%.15g,%.15g,' ...
        'pooled-labeled-only-CV-among-objective-admissible-candidates,' ...
        '%d-%d,1\n'],alpha,summary(index).ACC_mean,summary(index).NMI_mean, ...
        seedStart,seedStart+runs-1);
else
    fprintf(fid,'NEU-CLS,NaN,NaN,NaN,no-admissible-candidate,%d-%d,0\n', ...
        seedStart,seedStart+runs-1);
end
fclose(fid);
end

function neu_erd_write_cv_decision(path,runs,fractions,alphas,summary, ...
    selectedAlpha,selectedIndex,ready)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'ERDNMF NEU-CLS LABELED-ONLY ALPHA-SELECTION AUDIT\n');
fprintf(fid,'pilot_seeds=20260711-%d\n',20260711+runs-1);
fprintf(fid,'formal_seeds_excluded=20260731-20260750\n');
fprintf(fid,'label_rates=0.05 0.10 0.20\n');
fprintf(fid,'validation=five-fold stratified validation inside labeled samples\n');
fprintf(fid,'unlabeled_ground_truth_used=0\n');
fprintf(fid,'iterations=200\n');
fprintf(fid,'candidate_alphas='); fprintf(fid,'%g ',alphas); fprintf(fid,'\n');
fprintf(fid,['selection_rule=maximum pooled mean held-out ACC among candidates ' ...
    'passing all objective-nonincreasing audits; smaller alpha on exact tie\n\n']);
for i=1:length(summary)
    S=summary(i);
    fprintf(fid,['alpha=%g ACC=%.6f+-%.6f NMI=%.6f+-%.6f ' ...
        'ACC05=%.6f ACC10=%.6f ACC20=%.6f train_ACC=%.6f ' ...
        'objective_nonincrease=%d/%d time=%.3fs\n'], ...
        S.alpha,S.ACC_mean,S.ACC_sd,S.NMI_mean,S.NMI_sd, ...
        S.ACC_rate05,S.ACC_rate10,S.ACC_rate20,S.train_ACC_mean, ...
        S.mono_pass,S.evaluations,S.time_mean);
end
if ready
    fprintf(fid,'\nSELECTED_ALPHA=%g\n',selectedAlpha);
    fprintf(fid,'SELECTED_ACC=%.6f\n',summary(selectedIndex).ACC_mean);
else
    fprintf(fid,'\nSELECTED_ALPHA=NaN\n');
end
fprintf(fid,'FINITE_OUTPUT_PASS=1\n');
fprintf(fid,'OBJECTIVE_ADMISSIBLE_CANDIDATE_PASS=%d\n',ready);
fprintf(fid,'FINAL_READY=%d\n',ready && runs==20);
fclose(fid);
end
