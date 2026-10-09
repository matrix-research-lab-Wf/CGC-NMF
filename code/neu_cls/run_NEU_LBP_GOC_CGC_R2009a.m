function run_NEU_LBP_GOC_CGC_R2009a(runs,dataDir,outputDir)
%RUN_NEU_LBP_GOC_CGC_R2009A
% Formal paired comparison of GOCNMF and CGC-GOCNMF on NEU-CLS using
% nonnegative spatial uniform-LBP features.
%
% Fixed protocol:
%   feature       = LBP_{8,1}^{u2}, 4-by-4 cells, 944 dimensions
%   label rates   = 5%, 10%, 20% per class
%   seeds         = 20260731 onward
%   evaluation    = unlabeled samples only
%   assignment    = direct row-wise argmax; no K-means
%   GOC/CGC       = p=5, alpha=10, 50 iterations
%   calibration   = five stratified folds, kappa=1
%
% Usage:
%   run_NEU_LBP_GOC_CGC_R2009a(3,pwd,pwd);
%   run_NEU_LBP_GOC_CGC_R2009a(20,pwd,pwd);
%
% MATLAB R2009a compatible.

if nargin<1 || isempty(runs), runs=3; end
if nargin<2 || isempty(dataDir), dataDir=pwd; end
if nargin<3 || isempty(outputDir), outputDir=pwd; end

if ~(runs==3 || runs==20)
    error('runs must be 3 or 20.');
end
if exist(outputDir,'dir')~=7
    mkdir(outputDir);
end

EPSILON=1e-12;
TOL=1e-12;
LABEL_FRACTIONS=[0.05 0.10 0.20];
P=5;
ALPHA=10;
ITERATIONS=50;
FOLDS=5;
SEED_START=20260731;
BLOCK_SIZE=256;
methodNames={'GOCNMF','CGC-GOCNMF'};

dataPath=neu_lbp_find_file_recursive(dataDir, ...
    {'NEU_CLS_LBP59_4x4.mat'});
S=load(dataPath);
if ~isfield(S,'fea') || ~isfield(S,'gnd')
    error('%s must contain fea and gnd.',dataPath);
end

F=double(S.fea);
y0=double(S.gnd(:));
if size(F,1)==length(y0)
    X=F';
elseif size(F,2)==length(y0)
    X=F;
else
    error('Feature shape does not match labels.');
end

[dummy1,dummy2,y]=unique(y0); %#ok<ASGLU>
y=double(y(:));

if size(X,2)~=1800 || length(unique(y))~=6
    error('Expected 1800 samples and 6 classes.');
end
if size(X,1)~=944
    error('Expected a 944-dimensional LBP feature, obtained %d.',size(X,1));
end
if min(X(:))<-TOL || any(~isfinite(X(:)))
    error('Features must be finite and nonnegative.');
end

norms=sqrt(sum(X.^2,1));
norms(norms==0)=1;
X=bsxfun(@rdivide,X,norms);

[m,n]=size(X);
c=length(unique(y));

fprintf('\nNEU-CLS LBP: GOCNMF VERSUS CGC-GOCNMF\n');
fprintf('Data: %s\n',dataPath);
fprintf('Feature dimension: %d\n',m);
fprintf('Runs per label rate: %d\n',runs);
fprintf('Rates: 5%%, 10%%, 20%%; p=%d; alpha=%g; iterations=%d\n', ...
    P,ALPHA,ITERATIONS);
fprintf('Direct argmax; no K-means; unlabeled evaluation only.\n\n');

clockGraph=tic;
[neighborIndex,neighborDistance]=neu_lbp_knn(X,P,BLOCK_SIZE);
W0=neu_lbp_binary_graph(neighborIndex,n);
W1=neu_lbp_local_graph( ...
    neighborIndex,neighborDistance,n,EPSILON);
graphSeconds=toc(clockGraph);
fprintf('Candidate graphs constructed in %.3f seconds.\n',graphSeconds);

numberOfRates=length(LABEL_FRACTIONS);
numberOfMethods=length(methodNames);
totalRows=numberOfRates*runs*numberOfMethods;

prefix=sprintf('NEU_LBP_GOC_CGC_%dseed',runs);
rawPath=fullfile(outputDir,[prefix '_raw.csv']);
summaryPath=fullfile(outputDir,[prefix '_summary.csv']);
pairedPath=fullfile(outputDir,[prefix '_paired_tests.csv']);
decisionPath=fullfile(outputDir,[prefix '_decision.txt']);
latexPath=fullfile(outputDir,[prefix '_table.tex']);
checkpointPath=fullfile(outputDir,[prefix '_checkpoint.mat']);

template=struct( ...
    'label_fraction',0,'seed',0,'method','', ...
    'labeled_count',0,'theta',NaN,'delta_cv',NaN,'se_cv',NaN, ...
    'fallback',0,'ACC',0,'NMI',0,'ARI',0,'labeled_ACC',0, ...
    'objective_initial',0,'objective_final',0, ...
    'objective_nonincrease',0,'run_seconds',0);

rows=repmat(template,1,totalRows);
completed=false(numberOfRates,runs);
confusionCounts=zeros(numberOfRates,numberOfMethods,c,c);

if exist(checkpointPath,'file')==2
    Pcheckpoint=load(checkpointPath);
    if isfield(Pcheckpoint,'rows') && ...
            isfield(Pcheckpoint,'completed') && ...
            isfield(Pcheckpoint,'confusionCounts') && ...
            length(Pcheckpoint.rows)==totalRows && ...
            all(size(Pcheckpoint.completed)==[numberOfRates runs])
        rows=Pcheckpoint.rows;
        completed=Pcheckpoint.completed;
        confusionCounts=Pcheckpoint.confusionCounts;
        fprintf('Resuming checkpoint: %s\n',checkpointPath);
    else
        error('Existing checkpoint is incompatible with this run count.');
    end
end

for q=1:numberOfRates
    fraction=LABEL_FRACTIONS(q);

    for r=1:runs
        if completed(q,r)
            continue;
        end

        seed=SEED_START+r-1;
        fprintf('\nRate %.2f, seed %d (%d/%d)\n', ...
            fraction,seed,r,runs);

        L=neu_lbp_labeled_indices(y,seed,fraction);
        [U0,V0]=neu_lbp_initial_factors(m,n,c,seed,EPSILON);

        [theta,deltaCV,seCV]=neu_lbp_shrinkage_weight( ...
            W0,W1,L,y,c,seed,FOLDS,EPSILON);
        Wtheta=(1-theta)*W0+theta*W1;

        timer=tic;
        [Vgoc,gocInitial,gocFinal,gocMonotone]=neu_lbp_train_goc( ...
            X,y,L,W0,ALPHA,U0,V0,ITERATIONS,EPSILON);
        timeGoc=toc(timer);
        eGoc=neu_lbp_evaluate(Vgoc,y,L,EPSILON);

        if theta<=TOL
            Vcgc=Vgoc;
            cgcInitial=gocInitial;
            cgcFinal=gocFinal;
            cgcMonotone=gocMonotone;
            timeCgc=0;
        else
            timer=tic;
            [Vcgc,cgcInitial,cgcFinal,cgcMonotone]=neu_lbp_train_goc( ...
                X,y,L,Wtheta,ALPHA,U0,V0,ITERATIONS,EPSILON);
            timeCgc=toc(timer);
        end
        eCgc=neu_lbp_evaluate(Vcgc,y,L,EPSILON);

        evaluations={eGoc,eCgc};
        factors={Vgoc,Vcgc};
        initialObjectives=[gocInitial cgcInitial];
        finalObjectives=[gocFinal cgcFinal];
        monotoneFlags=[gocMonotone cgcMonotone];
        runTimes=[timeGoc timeCgc];

        baseIndex=((q-1)*runs+(r-1))*numberOfMethods;
        for k=1:numberOfMethods
            row=template;
            row.label_fraction=fraction;
            row.seed=seed;
            row.method=methodNames{k};
            row.labeled_count=length(L);
            row.delta_cv=deltaCV;
            row.se_cv=seCV;

            if k==1
                row.theta=0;
                row.fallback=1;
            else
                row.theta=theta;
                row.fallback=theta<=TOL;
            end

            row.ACC=evaluations{k}.ACC;
            row.NMI=evaluations{k}.NMI;
            row.ARI=evaluations{k}.ARI;
            row.labeled_ACC=evaluations{k}.labeledACC;
            row.objective_initial=initialObjectives(k);
            row.objective_final=finalObjectives(k);
            row.objective_nonincrease=monotoneFlags(k);
            row.run_seconds=runTimes(k);
            rows(baseIndex+k)=row;

            updated=squeeze(confusionCounts(q,k,:,:))+ ...
                neu_lbp_confusion(factors{k},y,L,c);
            confusionCounts(q,k,:,:)=reshape(updated,[1 1 c c]);
        end

        completed(q,r)=true;
        save(checkpointPath,'rows','completed','confusionCounts', ...
            'LABEL_FRACTIONS','methodNames','graphSeconds','P','ALPHA');

        neu_lbp_write_raw( ...
            rawPath,rows,completed,runs,numberOfMethods);

        fprintf(['GOC=%.4f/%.4f/%.4f; ' ...
            'CGC=%.4f/%.4f/%.4f; theta=%.4f\n'], ...
            eGoc.ACC,eGoc.NMI,eGoc.ARI, ...
            eCgc.ACC,eCgc.NMI,eCgc.ARI,theta);
    end
end

if ~all(completed(:))
    error('The experiment is incomplete.');
end

[summary,tests]=neu_lbp_summarize( ...
    rows,LABEL_FRACTIONS,runs,methodNames,TOL);

neu_lbp_write_raw(rawPath,rows,completed,runs,numberOfMethods);
neu_lbp_write_summary(summaryPath,summary);
neu_lbp_write_paired_tests(pairedPath,tests);
neu_lbp_write_confusions( ...
    outputDir,prefix,confusionCounts,methodNames,LABEL_FRACTIONS);
neu_lbp_write_latex(latexPath,summary,tests);
neu_lbp_write_decision( ...
    decisionPath,runs,rows,summary,tests,TOL,graphSeconds);
save(checkpointPath,'rows','completed','confusionCounts', ...
    'summary','tests','LABEL_FRACTIONS','methodNames', ...
    'graphSeconds','P','ALPHA');

fprintf('\nNEU-CLS LBP comparison complete.\n');
fprintf('RAW=%s\n',rawPath);
fprintf('SUMMARY=%s\n',summaryPath);
fprintf('PAIRED_TESTS=%s\n',pairedPath);
fprintf('TABLE=%s\n',latexPath);
fprintf('DECISION=%s\n',decisionPath);

if runs==3
    fprintf('NEU_LBP_THREE_SEED_INTEGRITY_COMPLETE=1\n');
else
    fprintf('NEU_LBP_TWENTY_SEED_FORMAL_COMPLETE=1\n');
end
end


function path=neu_lbp_find_file_recursive(folder,candidates)
path='';

for i=1:length(candidates)
    candidate=fullfile(folder,candidates{i});
    if exist(candidate,'file')==2
        path=candidate;
        return;
    end
end

entries=dir(folder);
for i=1:length(entries)
    if entries(i).isdir && ...
            ~strcmp(entries(i).name,'.') && ...
            ~strcmp(entries(i).name,'..')
        child=fullfile(folder,entries(i).name);
        path=neu_lbp_find_file_recursive(child,candidates);
        if ~isempty(path)
            return;
        end
    end
end

error('NEU_CLS_LBP59_4x4.mat was not found under %s.',folder);
end


function neu_lbp_set_seed(seed)
rand('twister',double(seed)); %#ok<RAND>
end


function L=neu_lbp_labeled_indices(y,seed,fraction)
neu_lbp_set_seed(seed);
classes=unique(y);
L=zeros(0,1);

for k=1:length(classes)
    ids=find(y==classes(k));
    order=randperm(length(ids));
    count=max(2,floor(fraction*length(ids)));
    count=min(count,length(ids));
    L=[L;ids(order(1:count))]; %#ok<AGROW>
end
L=double(L(:));
end


function [U0,V0]=neu_lbp_initial_factors(m,n,c,seed,epsilon)
neu_lbp_set_seed(seed+500001);
U0=max(rand(m,c),epsilon);
V0=max(rand(n,c),epsilon);
end


function [neighborIndex,neighborDistance]=neu_lbp_knn(X,p,blockSize)
n=size(X,2);
neighborIndex=zeros(n,p);
neighborDistance=zeros(n,p);

for first=1:blockSize:n
    last=min(first+blockSize-1,n);
    ids=first:last;
    similarities=X(:,ids)'*X;

    for q=1:length(ids)
        similarities(q,ids(q))=-Inf;
    end

    [sortedSimilarity,sortedIndex]=sort(similarities,2,'descend');
    neighborIndex(ids,:)=sortedIndex(:,1:p);
    neighborDistance(ids,:)=max(0,1-sortedSimilarity(:,1:p));

    fprintf('  kNN block %d:%d of %d\n',first,last,n);
    clear similarities sortedSimilarity sortedIndex;
end
end


function W=neu_lbp_binary_graph(neighborIndex,n)
p=size(neighborIndex,2);
rows=repmat((1:n)',1,p);
A=sparse(rows(:),neighborIndex(:),1,n,n);
W=spones(A+A');
W=W-spdiags(diag(W),0,n,n);
W=sparse(W);
end


function W=neu_lbp_local_graph( ...
    neighborIndex,neighborDistance,n,epsilon)

p=size(neighborIndex,2);
rows=repmat((1:n)',1,p);
cols=neighborIndex;
sigma=max(neighborDistance(:,p),epsilon);
denominator=sqrt(sigma(rows(:)).*sigma(cols(:)))+epsilon;
weights=exp(-(neighborDistance(:)./denominator).^2);

A=sparse(rows(:),cols(:),weights,n,n);
W=max(A,A');
W=W-spdiags(diag(W),0,n,n);

values=nonzeros(W);
if isempty(values)
    error('The locally weighted graph is empty.');
end
W=W/mean(values);
W=sparse(W);

W0=neu_lbp_binary_graph(neighborIndex,n);
if nnz(W)~=nnz(W0)
    error('The two candidate graphs do not have identical support.');
end
end


function folds=neu_lbp_stratified_folds(L,y,seed,numberOfFolds)
neu_lbp_set_seed(seed);
folds=cell(numberOfFolds,1);
for f=1:numberOfFolds
    folds{f}=zeros(0,1);
end

classes=unique(y(L));
for k=1:length(classes)
    ids=L(y(L)==classes(k));
    ids=ids(randperm(length(ids)));
    classOffset=mod(k-1,numberOfFolds);

    for j=1:length(ids)
        f=mod((j-1)+classOffset,numberOfFolds)+1;
        folds{f}=[folds{f};ids(j)]; %#ok<AGROW>
    end
end

for f=1:numberOfFolds
    if isempty(folds{f})
        error('Cross-validation fold %d is empty.',f);
    end
end
end


function P=neu_lbp_harmonic_predictions(W,train,y,c,query,epsilon)
n=size(W,1);
unknownMask=true(n,1);
unknownMask(train)=false;
unknown=find(unknownMask);

position=zeros(n,1);
position(unknown)=1:length(unknown);

degree=full(sum(W,2));
A=spdiags(degree(unknown),0,length(unknown),length(unknown))- ...
    W(unknown,unknown);

C=zeros(length(train),c);
index=sub2ind(size(C),(1:length(train))',y(train));
C(index)=1;
rhs=W(unknown,train)*C;

lastwarn('');
try
    Fu=A\rhs;
catch
    Fu=(A+1e-10*speye(length(unknown)))\rhs;
end
[warningMessage,dummyWarningId]=lastwarn; %#ok<ASGLU>

if any(~isfinite(Fu(:))) || ...
        ~isempty(strfind(lower(warningMessage),'singular'))
    Fu=(A+1e-10*speye(length(unknown)))\rhs;
end

Fu=max(Fu,0);
rowSum=sum(Fu,2);
positive=rowSum>epsilon;
if any(positive)
    Fu(positive,:)=bsxfun(@rdivide,Fu(positive,:),rowSum(positive));
end

P=Fu(position(query),:);
end


function losses=neu_lbp_fold_losses(W,folds,L,y,c,epsilon)
losses=zeros(length(folds),1);

for f=1:length(folds)
    heldOut=folds{f};
    train=setdiff(L,heldOut);
    prediction=neu_lbp_harmonic_predictions( ...
        W,train,y,c,heldOut,epsilon);

    target=zeros(length(heldOut),c);
    index=sub2ind(size(target),(1:length(heldOut))',y(heldOut));
    target(index)=1;
    losses(f)=mean(sum((prediction-target).^2,2));
end
end


function [theta,delta,se]=neu_lbp_shrinkage_weight( ...
    W0,W1,L,y,c,seed,numberOfFolds,epsilon)

folds=neu_lbp_stratified_folds(L,y,seed+991,numberOfFolds);
loss0=neu_lbp_fold_losses(W0,folds,L,y,c,epsilon);
loss1=neu_lbp_fold_losses(W1,folds,L,y,c,epsilon);

if any(~isfinite(loss0)) || any(~isfinite(loss1))
    error('Nonfinite foldwise prediction loss.');
end

difference=loss0-loss1;
delta=mean(difference);

if length(difference)>1
    se=std(difference,0)/sqrt(length(difference));
else
    se=0;
end

if delta>0
    theta=max(0,(delta-se)/max(delta,epsilon));
else
    theta=0;
end
theta=min(1,max(0,theta));
end


function [V,obj0,obj1,monotone]=neu_lbp_train_goc( ...
    X,y,L,W,alpha,U0,V0,iterations,epsilon)

U=max(U0,epsilon);
V=max(V0,epsilon);
c=size(V,2);

C=zeros(length(L),c);
index=sub2ind(size(C),(1:length(L))',y(L));
C(index)=1;
V(L,:)=C;

degree=full(sum(W,2));
obj0=neu_lbp_objective(X,U,V,W,degree,alpha);
previous=obj0;
monotone=1;

for iter=1:iterations
    numeratorU=X*V;
    denominatorU=U*(V'*V);
    U=U.*(numeratorU./max(denominatorU,epsilon));

    numeratorV=X'*U+alpha*(W*V);
    denominatorV=V*(U'*U)+ ...
        alpha*bsxfun(@times,degree,V);
    V=V.*(numeratorV./max(denominatorV,epsilon));
    V(L,:)=C;

    if iter==1 || mod(iter,10)==0 || iter==iterations
        current=neu_lbp_objective(X,U,V,W,degree,alpha);
        if current>previous+1e-8*max(1,abs(previous))
            monotone=0;
        end
        previous=current;
    end
end

obj1=neu_lbp_objective(X,U,V,W,degree,alpha);
end


function value=neu_lbp_objective(X,U,V,W,degree,alpha)
residual=X-U*V';
reconstruction=sum(residual(:).^2);
graphTerm=sum(sum(bsxfun(@times,degree,V).*V))- ...
    sum(sum((W*V).*V));
value=reconstruction+alpha*graphTerm;
end


function E=neu_lbp_evaluate(V,y,L,epsilon)
[dummy,prediction]=max(V,[],2); %#ok<ASGLU>
mask=true(length(y),1);
mask(L)=false;
U=find(mask);

E.ACC=mean(prediction(U)==y(U));
E.NMI=neu_lbp_nmi(y(U),prediction(U),epsilon);
E.ARI=neu_lbp_ari(y(U),prediction(U),epsilon);
E.labeledACC=mean(prediction(L)==y(L));
end


function value=neu_lbp_nmi(trueLabel,predictedLabel,epsilon)
trueLabel=trueLabel(:);
predictedLabel=predictedLabel(:);

[trueValues,dummy1,trueIndex]=unique(trueLabel); %#ok<ASGLU>
[predValues,dummy2,predIndex]=unique(predictedLabel); %#ok<ASGLU>

contingency=accumarray( ...
    [trueIndex predIndex],1,[length(trueValues) length(predValues)]);
P=contingency/length(trueLabel);
pi=sum(P,2);
pj=sum(P,1);
mutualInformation=0;

for i=1:size(P,1)
    for j=1:size(P,2)
        if P(i,j)>0
            mutualInformation=mutualInformation+ ...
                P(i,j)*log(P(i,j)/max(pi(i)*pj(j),epsilon));
        end
    end
end

ht=-sum(pi(pi>0).*log(pi(pi>0)));
hp=-sum(pj(pj>0).*log(pj(pj>0)));
denominator=0.5*(ht+hp);

if denominator<=epsilon
    value=1;
else
    value=mutualInformation/denominator;
end
value=max(0,min(1,value));
end


function value=neu_lbp_ari(trueLabel,predictedLabel,epsilon)
trueLabel=trueLabel(:);
predictedLabel=predictedLabel(:);

[trueValues,dummy1,trueIndex]=unique(trueLabel); %#ok<ASGLU>
[predValues,dummy2,predIndex]=unique(predictedLabel); %#ok<ASGLU>

C=accumarray([trueIndex predIndex],1, ...
    [length(trueValues) length(predValues)]);
a=sum(C,2);
b=sum(C,1);
N=length(trueLabel);

sumComb=sum(C(:).*(C(:)-1)/2);
rowComb=sum(a.*(a-1)/2);
colComb=sum(b.*(b-1)/2);
totalComb=N*(N-1)/2;

if totalComb<=epsilon
    value=1;
    return;
end

expected=rowComb*colComb/totalComb;
maximum=0.5*(rowComb+colComb);
denominator=maximum-expected;

if abs(denominator)<=epsilon
    value=1;
else
    value=(sumComb-expected)/denominator;
end
value=max(-1,min(1,value));
end


function C=neu_lbp_confusion(V,y,L,c)
[dummy,prediction]=max(V,[],2); %#ok<ASGLU>
mask=true(length(y),1);
mask(L)=false;
C=accumarray([y(mask) prediction(mask)],1,[c c]);
end


function neu_lbp_write_raw(path,rows,completed,runs,K)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end

fprintf(fid,['label_fraction,seed,method,labeled_count,theta,' ...
    'delta_cv,se_cv,fallback,ACC,NMI,ARI,labeled_ACC,' ...
    'objective_initial,objective_final,objective_nonincrease,' ...
    'run_seconds\n']);

for q=1:size(completed,1)
    for r=1:runs
        if completed(q,r)
            base=((q-1)*runs+(r-1))*K;
            for k=1:K
                S=rows(base+k);
                fprintf(fid,['%.6f,%d,%s,%d,%.15g,%.15g,%.15g,%d,' ...
                    '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%.6f\n'], ...
                    S.label_fraction,S.seed,S.method,S.labeled_count, ...
                    S.theta,S.delta_cv,S.se_cv,S.fallback, ...
                    S.ACC,S.NMI,S.ARI,S.labeled_ACC, ...
                    S.objective_initial,S.objective_final, ...
                    S.objective_nonincrease,S.run_seconds);
            end
        end
    end
end
fclose(fid);
end


function [summary,tests]=neu_lbp_summarize( ...
    rows,fractions,runs,methodNames,tol)

summaryTemplate=struct( ...
    'label_fraction',0,'method','', ...
    'ACC_mean',0,'ACC_sd',0,'NMI_mean',0,'NMI_sd',0, ...
    'ARI_mean',0,'ARI_sd',0,'theta_mean',NaN,'theta_sd',NaN, ...
    'fallbacks',0,'runs',runs);

summary=repmat(summaryTemplate,1,length(fractions)*length(methodNames));
counter=0;

for q=1:length(fractions)
    for k=1:length(methodNames)
        selected=false(1,length(rows));
        for i=1:length(rows)
            selected(i)=abs(rows(i).label_fraction-fractions(q))<1e-14 && ...
                strcmp(rows(i).method,methodNames{k});
        end
        R=rows(selected);
        if length(R)~=runs
            error('Summary selection failed.');
        end

        counter=counter+1;
        S=summaryTemplate;
        S.label_fraction=fractions(q);
        S.method=methodNames{k};
        S.ACC_mean=mean([R.ACC]);
        S.ACC_sd=neu_lbp_std([R.ACC]);
        S.NMI_mean=mean([R.NMI]);
        S.NMI_sd=neu_lbp_std([R.NMI]);
        S.ARI_mean=mean([R.ARI]);
        S.ARI_sd=neu_lbp_std([R.ARI]);

        if strcmp(methodNames{k},'CGC-GOCNMF')
            S.theta_mean=mean([R.theta]);
            S.theta_sd=neu_lbp_std([R.theta]);
            S.fallbacks=sum([R.theta]<=tol);
        end
        summary(counter)=S;
    end
end

metrics={'ACC','NMI','ARI'};
testTemplate=struct( ...
    'label_fraction',0,'metric','', ...
    'delta_mean',0,'delta_sd',0, ...
    'wins',0,'ties',0,'losses',0, ...
    'p_raw',1,'p_holm',1);

tests=repmat(testTemplate,1,length(fractions)*length(metrics));
counter=0;
rawP=zeros(length(tests),1);

for q=1:length(fractions)
    goc=neu_lbp_select_rows(rows,fractions(q),'GOCNMF');
    cgc=neu_lbp_select_rows(rows,fractions(q),'CGC-GOCNMF');

    for j=1:length(metrics)
        counter=counter+1;
        metric=metrics{j};
        baseValues=zeros(runs,1);
        newValues=zeros(runs,1);

        for r=1:runs
            baseValues(r)=goc(r).(metric);
            newValues(r)=cgc(r).(metric);
        end

        difference=newValues-baseValues;
        T=testTemplate;
        T.label_fraction=fractions(q);
        T.metric=metric;
        T.delta_mean=mean(difference);
        T.delta_sd=neu_lbp_std(difference);
        T.wins=sum(difference>tol);
        T.ties=sum(abs(difference)<=tol);
        T.losses=sum(difference<-tol);
        T.p_raw=neu_lbp_exact_wilcoxon(difference,tol);
        rawP(counter)=T.p_raw;
        tests(counter)=T;
    end
end

adjusted=neu_lbp_holm(rawP);
for i=1:length(tests)
    tests(i).p_holm=adjusted(i);
end
end


function R=neu_lbp_select_rows(rows,fraction,method)
selected=false(1,length(rows));
for i=1:length(rows)
    selected(i)=abs(rows(i).label_fraction-fraction)<1e-14 && ...
        strcmp(rows(i).method,method);
end
R=rows(selected);
[dummy,order]=sort([R.seed]); %#ok<ASGLU>
R=R(order);
end


function p=neu_lbp_exact_wilcoxon(difference,tol)
difference=difference(:);
difference=difference(abs(difference)>tol);
number=length(difference);

if number==0
    p=1;
    return;
end

ranks=neu_lbp_tied_ranks(abs(difference));
rank2=round(2*ranks);
totalRank2=sum(rank2);
positiveRank2=sum(rank2(difference>0));
observed=min(positiveRank2,totalRank2-positiveRank2);

counts=zeros(1,totalRank2+1);
counts(1)=1;

for i=1:number
    shift=rank2(i);
    previous=counts;
    counts(shift+1:end)=counts(shift+1:end)+ ...
        previous(1:end-shift);
end

possible=0:totalRank2;
extreme=min(possible,totalRank2-possible)<=observed;
p=sum(counts(extreme))/(2^number);
p=min(1,max(0,p));
end


function ranks=neu_lbp_tied_ranks(values)
[sorted,order]=sort(values);
ranks=zeros(size(values));
i=1;
tolerance=1e-12;

while i<=length(sorted)
    j=i;
    scale=max(1,abs(sorted(i)));
    while j<length(sorted) && ...
            abs(sorted(j+1)-sorted(i))<=tolerance*scale
        j=j+1;
    end
    averageRank=0.5*(i+j);
    ranks(order(i:j))=averageRank;
    i=j+1;
end
end


function adjusted=neu_lbp_holm(raw)
raw=raw(:);
number=length(raw);
[sorted,order]=sort(raw);
adjustedSorted=zeros(number,1);
previous=0;

for i=1:number
    current=min(1,(number-i+1)*sorted(i));
    current=max(previous,current);
    adjustedSorted(i)=current;
    previous=current;
end

adjusted=zeros(number,1);
adjusted(order)=adjustedSorted;
end


function neu_lbp_write_summary(path,summary)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end

fprintf(fid,['feature,label_fraction,method,runs,' ...
    'ACC_mean,ACC_sd,NMI_mean,NMI_sd,ARI_mean,ARI_sd,' ...
    'theta_mean,theta_sd,fallbacks\n']);

for i=1:length(summary)
    S=summary(i);
    fprintf(fid,['LBP59_4x4,%.6f,%s,%d,%.15g,%.15g,' ...
        '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d\n'], ...
        S.label_fraction,S.method,S.runs, ...
        S.ACC_mean,S.ACC_sd,S.NMI_mean,S.NMI_sd, ...
        S.ARI_mean,S.ARI_sd,S.theta_mean,S.theta_sd,S.fallbacks);
end
fclose(fid);
end


function neu_lbp_write_paired_tests(path,tests)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end

fprintf(fid,['feature,label_fraction,comparison,metric,delta_mean,' ...
    'delta_sd,wins,ties,losses,p_raw,p_holm\n']);

for i=1:length(tests)
    T=tests(i);
    fprintf(fid,['LBP59_4x4,%.6f,CGC-GOCNMF_minus_GOCNMF,%s,' ...
        '%.15g,%.15g,%d,%d,%d,%.15g,%.15g\n'], ...
        T.label_fraction,T.metric,T.delta_mean,T.delta_sd, ...
        T.wins,T.ties,T.losses,T.p_raw,T.p_holm);
end
fclose(fid);
end


function neu_lbp_write_confusions( ...
    outputDir,prefix,counts,methodNames,fractions)

for q=1:length(fractions)
    for k=1:length(methodNames)
        C=squeeze(counts(q,k,:,:));
        rowSum=sum(C,2);
        rowSum(rowSum==0)=1;
        normalized=bsxfun(@rdivide,C,rowSum);

        rate=round(100*fractions(q));
        name=strrep(methodNames{k},'-','');
        countPath=fullfile(outputDir,sprintf( ...
            '%s_confusion_counts_label%d_%s.csv',prefix,rate,name));
        normalizedPath=fullfile(outputDir,sprintf( ...
            '%s_confusion_normalized_label%d_%s.csv', ...
            prefix,rate,name));

        csvwrite(countPath,C);
        csvwrite(normalizedPath,normalized);
    end
end
end


function neu_lbp_write_latex(path,summary,tests)
fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end

fprintf(fid,'\\begin{table}[t]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,['\\caption{GOCNMF and CGC-GOCNMF on NEU-CLS using ' ...
    'nonnegative spatial uniform-LBP features.}\\label{tab:neu_lbp}\n']);
fprintf(fid,'\\small\n');
fprintf(fid,'\\begin{tabular}{clccc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Label rate & Method & ACC (\\%%) & NMI (\\%%) & ARI (\\%%) \\\\\n');
fprintf(fid,'\\midrule\n');

fractions=unique([summary.label_fraction]);
for q=1:length(fractions)
    G=neu_lbp_find_summary(summary,fractions(q),'GOCNMF');
    C=neu_lbp_find_summary(summary,fractions(q),'CGC-GOCNMF');

    fprintf(fid,'\\multirow{2}{*}{%d\\%%} ',round(100*fractions(q)));
    fprintf(fid,'& GOCNMF & %.2f$\\pm$%.2f & %.2f$\\pm$%.2f & %.2f$\\pm$%.2f \\\\\n', ...
        100*G.ACC_mean,100*G.ACC_sd,100*G.NMI_mean,100*G.NMI_sd, ...
        100*G.ARI_mean,100*G.ARI_sd);
    fprintf(fid,['& CGC-GOCNMF & %.2f$\\pm$%.2f & %.2f$\\pm$%.2f ' ...
        '& %.2f$\\pm$%.2f \\\\\n'], ...
        100*C.ACC_mean,100*C.ACC_sd,100*C.NMI_mean,100*C.NMI_sd, ...
        100*C.ARI_mean,100*C.ARI_sd);

    if q<length(fractions)
        fprintf(fid,'\\midrule\n');
    end
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n\n');
fprintf(fid,'\\vspace{1mm}\n');
fprintf(fid,'\\parbox{0.98\\linewidth}{\\footnotesize\n');
fprintf(fid,['Values are mean percentages $\\pm$ sample standard ' ...
    'deviations over 20 paired trials on unlabeled samples. ' ...
    'Paired Wilcoxon tests are reported separately.}\n']);
fprintf(fid,'\\end{table}\n');

fclose(fid);
end


function S=neu_lbp_find_summary(summary,fraction,method)
found=0;
S=summary(1);

for i=1:length(summary)
    if abs(summary(i).label_fraction-fraction)<1e-14 && ...
            strcmp(summary(i).method,method)
        S=summary(i);
        found=1;
        break;
    end
end

if ~found
    error('Requested summary row was not found.');
end
end


function neu_lbp_write_decision( ...
    path,runs,rows,summary,tests,tol,graphSeconds)

finitePass=1;
thetaPass=1;
objectivePass=1;
labeledPass=1;
fallbackPass=1;

for i=1:length(rows)
    values=[rows(i).theta rows(i).delta_cv rows(i).se_cv ...
        rows(i).ACC rows(i).NMI rows(i).ARI ...
        rows(i).objective_initial rows(i).objective_final];
    if any(~isfinite(values))
        finitePass=0;
    end
    if rows(i).theta<-tol || rows(i).theta>1+tol
        thetaPass=0;
    end
    if rows(i).objective_nonincrease~=1
        objectivePass=0;
    end
    if abs(rows(i).labeled_ACC-1)>tol
        labeledPass=0;
    end
end

fractions=unique([rows.label_fraction]);
for q=1:length(fractions)
    G=neu_lbp_select_rows(rows,fractions(q),'GOCNMF');
    C=neu_lbp_select_rows(rows,fractions(q),'CGC-GOCNMF');
    for r=1:length(G)
        if C(r).theta<=tol
            if abs(C(r).ACC-G(r).ACC)>tol || ...
                    abs(C(r).NMI-G(r).NMI)>tol || ...
                    abs(C(r).ARI-G(r).ARI)>tol
                fallbackPass=0;
            end
        end
    end
end

fid=fopen(path,'w');
if fid<0, error('Cannot create %s.',path); end

fprintf(fid,'NEU-CLS LBP GOCNMF/CGC-GOCNMF audit\n');
fprintf(fid,'Feature: LBP_{8,1}^{u2}, 4-by-4 cells, 944 dimensions\n');
fprintf(fid,'Runs per label rate: %d\n',runs);
fprintf(fid,'Graph construction seconds: %.6f\n',graphSeconds);
fprintf(fid,'Direct row-wise argmax; K-means used: NO\n\n');

for q=1:length(fractions)
    G=neu_lbp_find_summary(summary,fractions(q),'GOCNMF');
    C=neu_lbp_find_summary(summary,fractions(q),'CGC-GOCNMF');

    fprintf(fid,['Label %d%% GOCNMF: ACC=%.6f+-%.6f, ' ...
        'NMI=%.6f+-%.6f, ARI=%.6f+-%.6f\n'], ...
        round(100*fractions(q)),G.ACC_mean,G.ACC_sd, ...
        G.NMI_mean,G.NMI_sd,G.ARI_mean,G.ARI_sd);

    fprintf(fid,['Label %d%% CGC-GOCNMF: ACC=%.6f+-%.6f, ' ...
        'NMI=%.6f+-%.6f, ARI=%.6f+-%.6f, ' ...
        'theta=%.6f+-%.6f, fallback=%d/%d\n'], ...
        round(100*fractions(q)),C.ACC_mean,C.ACC_sd, ...
        C.NMI_mean,C.NMI_sd,C.ARI_mean,C.ARI_sd, ...
        C.theta_mean,C.theta_sd,C.fallbacks,C.runs);

    metrics={'ACC','NMI','ARI'};
    for j=1:length(metrics)
        T=neu_lbp_find_test(tests,fractions(q),metrics{j});
        fprintf(fid,['  %s: delta=%+.6f, wins/ties/losses=%d/%d/%d, ' ...
            'p=%.8g, Holm=%.8g\n'], ...
            metrics{j},T.delta_mean,T.wins,T.ties,T.losses, ...
            T.p_raw,T.p_holm);
    end
    fprintf(fid,'\n');
end

fprintf(fid,'FINITE_OUTPUT_PASS=%d\n',finitePass);
fprintf(fid,'THETA_RANGE_PASS=%d\n',thetaPass);
fprintf(fid,'FIXED_LABELED_ROWS_PASS=%d\n',labeledPass);
fprintf(fid,'OBJECTIVE_NONINCREASE_PASS=%d\n',objectivePass);
fprintf(fid,'EXACT_FALLBACK_IDENTITY_PASS=%d\n',fallbackPass);

if runs==3
    fprintf(fid,'NEU_LBP_THREE_SEED_INTEGRITY_PASS=%d\n', ...
        finitePass && thetaPass && labeledPass && ...
        objectivePass && fallbackPass);
else
    fprintf(fid,'NEU_LBP_TWENTY_SEED_FORMAL_COMPLETE=%d\n', ...
        finitePass && thetaPass && labeledPass && ...
        objectivePass && fallbackPass);
end

fclose(fid);
end


function T=neu_lbp_find_test(tests,fraction,metric)
found=0;
T=tests(1);

for i=1:length(tests)
    if abs(tests(i).label_fraction-fraction)<1e-14 && ...
            strcmp(tests(i).metric,metric)
        T=tests(i);
        found=1;
        break;
    end
end

if ~found
    error('Requested paired test was not found.');
end
end


function value=neu_lbp_std(x)
if length(x)>1
    value=std(x,0);
else
    value=0;
end
end
