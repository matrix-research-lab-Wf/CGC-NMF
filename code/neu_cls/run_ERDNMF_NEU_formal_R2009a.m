function run_ERDNMF_NEU_formal_R2009a(runs,dataDir,outputDir,selectedFile)
%RUN_ERDNMF_NEU_FORMAL_R2009A
% Formal ERDNMF evaluation for Table 6 on NEU-CLS.
%
% Uses the same label rates and formal seeds as the existing NEU-CLS
% GNMFLD/GOCNMF/CGC-GOCNMF experiment:
%   label rates = 5%, 10%, 20%
%   seeds       = 20260731--20260750
%   evaluation  = unlabeled samples only
%   assignment  = row-wise argmax; no K-means
%
% Usage:
%   run_ERDNMF_NEU_formal_R2009a(1,pwd,pwd, ...
%       'ERDNMF_NEU_labeledCV_alpha_selection_20seed_selected.csv');
%   run_ERDNMF_NEU_formal_R2009a(20,pwd,pwd, ...
%       'ERDNMF_NEU_labeledCV_alpha_selection_20seed_selected.csv');

if nargin<1 || isempty(runs), runs=1; end
if nargin<2 || isempty(dataDir), dataDir=pwd; end
if nargin<3 || isempty(outputDir), outputDir=pwd; end
if nargin<4 || isempty(selectedFile)
    selectedFile='ERDNMF_NEU_labeledCV_alpha_selection_20seed_selected.csv';
end
if ~(runs==1 || runs==3 || runs==20)
    error('runs must be 1, 3, or 20.');
end
if exist(outputDir,'dir')~=7, mkdir(outputDir); end
if exist(selectedFile,'file')~=2
    selectedFile=fullfile(dataDir,selectedFile);
end
if exist(selectedFile,'file')~=2
    error('Selected-alpha file not found.');
end

EPSILON=1e-12;
LABEL_FRACTIONS=[0.05 0.10 0.20];
MAX_ITERATIONS=200;
FORMAL_SEED_START=20260731;
alpha=neu_erd_read_selected_alpha(selectedFile);
if ~isfinite(alpha) || alpha<=0
    error('Invalid selected alpha in %s.',selectedFile);
end

dataPath=neu_erd_find_file_recursive(dataDir, ...
    {'NEU_CLS_32x32.mat','NEU-CLS-32x32.mat','NEU_CLS.mat'});
[X,y]=neu_erd_load_dataset(dataPath,1800,6);
[m,n]=size(X);
c=length(unique(y));

prefix=sprintf('ERDNMF_NEU_formal_%dseed',runs);
rawPath=fullfile(outputDir,[prefix '_raw.csv']);
summaryPath=fullfile(outputDir,[prefix '_summary.csv']);
parametersPath=fullfile(outputDir,[prefix '_parameters.csv']);
decisionPath=fullfile(outputDir,[prefix '_decision.txt']);
checkpointPath=fullfile(outputDir,[prefix '_checkpoint.mat']);

nr=length(LABEL_FRACTIONS);
template=struct('label_fraction',0,'seed',0,'method','ERDNMF', ...
    'labeled_count',0,'alpha',alpha,'ACC',0,'NMI',0,'ARI',0, ...
    'labeled_ACC',0,'objective_initial',0,'objective_final',0, ...
    'objective_nonincrease',0,'run_seconds',0);
rows=repmat(template,1,nr*runs);
completed=false(nr,runs);
confusionCounts=zeros(nr,c,c);

if exist(checkpointPath,'file')==2
    P=load(checkpointPath);
    if isfield(P,'rows') && isfield(P,'completed') && ...
            isfield(P,'confusionCounts') && length(P.rows)==nr*runs && ...
            all(size(P.completed)==[nr runs])
        rows=P.rows;
        completed=P.completed;
        confusionCounts=P.confusionCounts;
        fprintf('Resuming checkpoint: %s\n',checkpointPath);
    else
        error('Existing checkpoint is incompatible with this run count.');
    end
end

fprintf('\nERDNMF NEU-CLS FORMAL TABLE-6 EXPERIMENT\n');
fprintf('Data: %s\n',dataPath);
fprintf('Selected alpha: %g\n',alpha);
fprintf('Seeds: %d to %d\n',FORMAL_SEED_START,FORMAL_SEED_START+runs-1);
fprintf('Label rates: 5%%, 10%%, 20%%; iterations=%d\n',MAX_ITERATIONS);
fprintf('Evaluation: unlabeled only; direct argmax; no K-means\n\n');

for q=1:nr
    fraction=LABEL_FRACTIONS(q);
    for r=1:runs
        if completed(q,r), continue; end
        seed=FORMAL_SEED_START+r-1;
        L=neu_erd_labeled_indices(y,seed,fraction);
        [U0,V0]=neu_erd_initial_factors(m,n,c,seed+500001,EPSILON);

        t0=tic;
        [V,obj0,obj1,nonincrease]=neu_erd_train( ...
            X,y,L,alpha,U0,V0,MAX_ITERATIONS,EPSILON);
        elapsed=toc(t0);
        E=neu_erd_evaluate(V,y,L,EPSILON);
        C=neu_erd_confusion(V,y,L,c);
        confusionCounts(q,:,:)=reshape(squeeze(confusionCounts(q,:,:))+C,[1 c c]);

        idx=(q-1)*runs+r;
        row=template;
        row.label_fraction=fraction;
        row.seed=seed;
        row.labeled_count=length(L);
        row.ACC=E.ACC;
        row.NMI=E.NMI;
        row.ARI=E.ARI;
        row.labeled_ACC=E.labeledACC;
        row.objective_initial=obj0;
        row.objective_final=obj1;
        row.objective_nonincrease=nonincrease;
        row.run_seconds=elapsed;
        rows(idx)=row;
        completed(q,r)=true;

        save(checkpointPath,'rows','completed','confusionCounts', ...
            'LABEL_FRACTIONS','alpha','runs');
        neu_erd_write_formal_raw(rawPath,rows,completed,runs);
        fprintf(['rate=%.2f seed=%d ACC=%.6f NMI=%.6f ARI=%.6f ' ...
            'ACC_L=%.6f mono=%d time=%.2fs\n'], ...
            fraction,seed,E.ACC,E.NMI,E.ARI,E.labeledACC, ...
            nonincrease,elapsed);
    end
end

summary=neu_erd_formal_summarize(rows,LABEL_FRACTIONS,runs);
neu_erd_write_formal_raw(rawPath,rows,completed,runs);
neu_erd_write_formal_summary(summaryPath,summary);
neu_erd_write_formal_parameters(parametersPath,alpha,runs);
neu_erd_write_formal_confusions(outputDir,prefix,confusionCounts,LABEL_FRACTIONS);
neu_erd_write_formal_decision(decisionPath,runs,alpha,summary,rows);

if all(completed(:)) && exist(checkpointPath,'file')==2
    delete(checkpointPath);
end
fprintf('\nERDNMF_NEU_FORMAL_COMPLETE=%d\n',all(completed(:)));
fprintf('RAW=%s\nSUMMARY=%s\nPARAMETERS=%s\nDECISION=%s\n', ...
    rawPath,summaryPath,parametersPath,decisionPath);
end

function alpha=neu_erd_read_selected_alpha(path)
fid=fopen(path,'r');
if fid<0, error('Cannot open %s.',path); end
header=fgetl(fid); %#ok<NASGU>
line=fgetl(fid);
fclose(fid);
if ~ischar(line), error('Selected-alpha file has no data row.'); end
parts=textscan(line,'%s%f%f%f%s%s%d','Delimiter',',');
if isempty(parts{2}), error('Cannot parse selected alpha.'); end
alpha=parts{2}(1);
if ~isempty(parts{7}) && parts{7}(1)~=1
    error('Selected-alpha file is not marked ready.');
end
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

function E=neu_erd_evaluate(V,y,L,epsilon)
[dummy,pred]=max(V,[],2); %#ok<ASGLU>
mask=true(length(y),1);
mask(L)=false;
U=find(mask);
E.ACC=mean(pred(U)==y(U));
E.NMI=neu_erd_nmi(y(U),pred(U),epsilon);
E.ARI=neu_erd_ari(y(U),pred(U),epsilon);
E.labeledACC=mean(pred(L)==y(L));
end

function C=neu_erd_confusion(V,y,L,c)
[dummy,pred]=max(V,[],2); %#ok<ASGLU>
mask=true(length(y),1);
mask(L)=false;
C=accumarray([y(mask) pred(mask)],1,[c c]);
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

function value=neu_erd_ari(trueLabel,predictedLabel,epsilon)
trueLabel=trueLabel(:);
predictedLabel=predictedLabel(:);
[tv,d1,ti]=unique(trueLabel); %#ok<ASGLU>
[pv,d2,pi]=unique(predictedLabel); %#ok<ASGLU>
C=accumarray([ti pi],1,[length(tv) length(pv)]);
a=sum(C,2);
b=sum(C,1);
N=length(trueLabel);
sumComb=sum(C(:).*(C(:)-1)/2);
rowComb=sum(a.*(a-1)/2);
colComb=sum(b.*(b-1)/2);
totalComb=N*(N-1)/2;
if totalComb<=epsilon, value=1; return; end
expected=rowComb*colComb/totalComb;
maximum=0.5*(rowComb+colComb);
den=maximum-expected;
if abs(den)<=epsilon, value=1; else, value=(sumComb-expected)/den; end
value=max(-1,min(1,value));
end

function neu_erd_write_formal_raw(path,rows,completed,runs)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['label_fraction,seed,method,labeled_count,alpha,ACC,NMI,ARI,' ...
    'labeled_ACC,objective_initial,objective_final,objective_nonincrease,' ...
    'run_seconds\n']);
for q=1:size(completed,1)
    for r=1:runs
        if completed(q,r)
            idx=(q-1)*runs+r;
            S=rows(idx);
            fprintf(fid,['%.6f,%d,%s,%d,%.15g,%.15g,%.15g,%.15g,' ...
                '%.15g,%.15g,%.15g,%d,%.6f\n'], ...
                S.label_fraction,S.seed,S.method,S.labeled_count,S.alpha, ...
                S.ACC,S.NMI,S.ARI,S.labeled_ACC,S.objective_initial, ...
                S.objective_final,S.objective_nonincrease,S.run_seconds);
        end
    end
end
fclose(fid);
end

function summary=neu_erd_formal_summarize(rows,fractions,runs)
template=struct('label_fraction',0,'runs',runs,'alpha',0, ...
    'ACC_mean',0,'ACC_sd',0,'NMI_mean',0,'NMI_sd',0, ...
    'ARI_mean',0,'ARI_sd',0,'labeled_ACC_mean',0, ...
    'mono_pass',0,'time_mean',0);
summary=repmat(template,1,length(fractions));
for q=1:length(fractions)
    ids=(q-1)*runs+(1:runs);
    A=[rows(ids).ACC];
    N=[rows(ids).NMI];
    R=[rows(ids).ARI];
    L=[rows(ids).labeled_ACC];
    O=[rows(ids).objective_nonincrease];
    T=[rows(ids).run_seconds];
    summary(q).label_fraction=fractions(q);
    summary(q).alpha=rows(ids(1)).alpha;
    summary(q).ACC_mean=mean(A);
    summary(q).ACC_sd=neu_erd_safe_std(A);
    summary(q).NMI_mean=mean(N);
    summary(q).NMI_sd=neu_erd_safe_std(N);
    summary(q).ARI_mean=mean(R);
    summary(q).ARI_sd=neu_erd_safe_std(R);
    summary(q).labeled_ACC_mean=mean(L);
    summary(q).mono_pass=sum(O==1);
    summary(q).time_mean=mean(T);
end
end

function v=neu_erd_safe_std(x)
if length(x)>1, v=std(x,0); else, v=0; end
end

function neu_erd_write_formal_summary(path,summary)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['label_fraction,method,runs,alpha,ACC_mean,ACC_sd,NMI_mean,' ...
    'NMI_sd,ARI_mean,ARI_sd,labeled_ACC_mean,' ...
    'objective_nonincrease_pass,run_seconds_mean\n']);
for q=1:length(summary)
    S=summary(q);
    fprintf(fid,['%.6f,ERDNMF,%d,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%d,%.6f\n'], ...
        S.label_fraction,S.runs,S.alpha,S.ACC_mean,S.ACC_sd, ...
        S.NMI_mean,S.NMI_sd,S.ARI_mean,S.ARI_sd,S.labeled_ACC_mean, ...
        S.mono_pass,S.time_mean);
end
fclose(fid);
end

function neu_erd_write_formal_parameters(path,alpha,runs)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,['dataset,method,alpha,iterations,label_rates,formal_seeds,' ...
    'assignment,evaluation,parameter_selection\n']);
fprintf(fid,['NEU-CLS,ERDNMF,%.15g,200,0.05|0.10|0.20,' ...
    '20260731-%d,row-wise-argmax,unlabeled-only,' ...
    'labeled-only-pilot-CV\n'],alpha,20260731+runs-1);
fclose(fid);
end

function neu_erd_write_formal_confusions(outputDir,prefix,counts,fractions)
for q=1:length(fractions)
    rateCode=round(100*fractions(q));
    C=squeeze(counts(q,:,:));
    countPath=fullfile(outputDir,sprintf( ...
        '%s_confusion_counts_label%02d_ERDNMF.csv',prefix,rateCode));
    normPath=fullfile(outputDir,sprintf( ...
        '%s_confusion_normalized_label%02d_ERDNMF.csv',prefix,rateCode));
    neu_erd_write_numeric_csv(countPath,C);
    rowSums=sum(C,2);
    rowSums(rowSums==0)=1;
    P=bsxfun(@rdivide,C,rowSums);
    neu_erd_write_numeric_csv(normPath,P);
end
end

function neu_erd_write_numeric_csv(path,M)
% Version-independent CSV writer for MATLAB R2009a.
fid=fopen(path,'w');
if fid<0
    error('Cannot create %s.',path);
end
[nr,nc]=size(M);
for i=1:nr
    if nc>0
        fprintf(fid,'%.15g',M(i,1));
        for j=2:nc
            fprintf(fid,',%.15g',M(i,j));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);
end

function neu_erd_write_formal_decision(path,runs,alpha,summary,rows)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end
fprintf(fid,'ERDNMF NEU-CLS FORMAL TABLE-6 AUDIT\n');
fprintf(fid,'formal_seeds=20260731-%d\n',20260731+runs-1);
fprintf(fid,'label_rates=0.05 0.10 0.20\n');
fprintf(fid,'alpha=%g\n',alpha);
fprintf(fid,'iterations=200\n');
fprintf(fid,'evaluation=unlabeled samples only\n');
fprintf(fid,'assignment=row-wise argmax; no K-means\n');
fprintf(fid,'unlabeled_truth_used_for_parameter_selection=0\n\n');
for q=1:length(summary)
    S=summary(q);
    fprintf(fid,['rate=%.2f ACC=%.6f+-%.6f NMI=%.6f+-%.6f ' ...
        'ARI=%.6f+-%.6f labeled_ACC=%.6f ' ...
        'objective_nonincrease=%d/%d time=%.3fs\n'], ...
        S.label_fraction,S.ACC_mean,S.ACC_sd,S.NMI_mean,S.NMI_sd, ...
        S.ARI_mean,S.ARI_sd,S.labeled_ACC_mean,S.mono_pass,runs, ...
        S.time_mean);
end
finitePass=1;
for i=1:length(rows)
    values=[rows(i).ACC rows(i).NMI rows(i).ARI rows(i).labeled_ACC ...
        rows(i).objective_initial rows(i).objective_final];
    if any(~isfinite(values)), finitePass=0; end
end
monoPass=all([summary.mono_pass]==runs);
fprintf(fid,'\nFINITE_OUTPUT_PASS=%d\n',finitePass);
fprintf(fid,'OBJECTIVE_NONINCREASE_PASS=%d\n',monoPass);
fprintf(fid,'FINAL_READY=%d\n',runs==20 && finitePass && monoPass);
fclose(fid);
end
