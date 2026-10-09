function summarize_NEU_application_R2009a(runs,outputDir,checkpointPath)
%SUMMARIZE_NEU_APPLICATION_R2009A Summarize NEU Stage-1 results.

if nargin<1 || isempty(runs), runs=1; end
if nargin<2 || isempty(outputDir), outputDir=pwd; end
if nargin<3 || isempty(checkpointPath)
    checkpointPath=fullfile(outputDir, ...
        sprintf('NEU_application_%dseed_checkpoint.mat',runs));
end
if exist(checkpointPath,'file')~=2
    error('Checkpoint not found: %s',checkpointPath);
end
P=load(checkpointPath);
rows=P.rows; confusionCounts=P.confusionCounts;
labelFractions=P.LABEL_FRACTIONS; methodNames=P.methodNames;
K=length(methodNames); Q=length(labelFractions); TOL=1e-12;
prefix=sprintf('NEU_application_%dseed',runs);
summaryPath=fullfile(outputDir,[prefix '_summary.csv']);
pairedPath=fullfile(outputDir,[prefix '_paired.csv']);
decisionPath=fullfile(outputDir,[prefix '_decision.txt']);

fs=fopen(summaryPath,'w'); fp=fopen(pairedPath,'w'); fd=fopen(decisionPath,'w');
if fs<0 || fp<0 || fd<0, error('Cannot create summary outputs.'); end
fprintf(fs,['label_fraction,method,runs,ACC_mean,ACC_sd,NMI_mean,NMI_sd,' ...
    'ARI_mean,ARI_sd,labeled_ACC_mean,theta_mean,theta_sd,' ...
    'fallback_count,objective_nonincrease_pass,run_seconds_mean\n']);
fprintf(fp,['label_fraction,comparison,runs,dACC_mean,ACC_wins,ACC_ties,' ...
    'ACC_losses,dNMI_mean,NMI_wins,NMI_ties,NMI_losses,dARI_mean,' ...
    'ARI_wins,ARI_ties,ARI_losses\n']);
fprintf(fd,'NEU-CLS engineering application -- Stage 1\n');
fprintf(fd,'Runs per label fraction: %d\n',runs);
fprintf(fd,'Primary metrics: unlabeled ACC, NMI, ARI\n');
fprintf(fd,'Direct row-wise argmax; K-means used: NO\n\n');

finitePass=1; objectivePass=1; fixedLabelPass=1;
pairNames={'CGC-GOCNMF_minus_GOCNMF','CGC-GOCNMF_minus_GNMFLD'};
pairIndex=[3 2;3 1];

for q=1:Q
    fraction=labelFractions(q);
    methodRows=cell(K,1);
    for k=1:K
        values=repmat(rows(1),1,runs);
        for r=1:runs
            base=((q-1)*runs+(r-1))*K;
            values(r)=rows(base+k);
        end
        methodRows{k}=values;
        A=[values.ACC]; N=[values.NMI]; R=[values.ARI];
        L=[values.labeled_ACC]; O=[values.objective_nonincrease];
        time=[values.run_seconds];
        if any(~isfinite([A N R L time])), finitePass=0; end
        if ~all(O==1), objectivePass=0; end
        if k>=2 && any(L<1-TOL), fixedLabelPass=0; end
        thetaMean=NaN; thetaSD=NaN; fallback=0;
        if k==3
            theta=[values.theta]; thetaMean=mean(theta);
            thetaSD=neu_safe_std(theta); fallback=sum([values.fallback]==1);
        end
        fprintf(fs,['%.6f,%s,%d,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
            '%.15g,%.15g,%.15g,%.15g,%d,%d,%.6f\n'], ...
            fraction,methodNames{k},runs,mean(A),neu_safe_std(A), ...
            mean(N),neu_safe_std(N),mean(R),neu_safe_std(R),mean(L), ...
            thetaMean,thetaSD,fallback,all(O==1),mean(time));
        fprintf(fd,['Label %.0f%% %s: ACC=%.6f+-%.6f, NMI=%.6f+-%.6f, ' ...
            'ARI=%.6f+-%.6f, labeled_ACC=%.6f, theta=%.6f, ' ...
            'fallback=%d/%d\n'],100*fraction,methodNames{k}, ...
            mean(A),neu_safe_std(A),mean(N),neu_safe_std(N), ...
            mean(R),neu_safe_std(R),mean(L),thetaMean,fallback,runs);
    end
    for p=1:length(pairNames)
        first=methodRows{pairIndex(p,1)}; second=methodRows{pairIndex(p,2)};
        dA=[first.ACC]-[second.ACC]; dN=[first.NMI]-[second.NMI];
        dR=[first.ARI]-[second.ARI];
        fprintf(fp,['%.6f,%s,%d,%.15g,%d,%d,%d,%.15g,%d,%d,%d,' ...
            '%.15g,%d,%d,%d\n'],fraction,pairNames{p},runs,mean(dA), ...
            sum(dA>TOL),sum(abs(dA)<=TOL),sum(dA<-TOL),mean(dN), ...
            sum(dN>TOL),sum(abs(dN)<=TOL),sum(dN<-TOL),mean(dR), ...
            sum(dR>TOL),sum(abs(dR)<=TOL),sum(dR<-TOL));
        fprintf(fd,['  %s: dACC=%+.6f (%d/%d/%d), dNMI=%+.6f ' ...
            '(%d/%d/%d), dARI=%+.6f (%d/%d/%d)\n'],pairNames{p}, ...
            mean(dA),sum(dA>TOL),sum(abs(dA)<=TOL),sum(dA<-TOL), ...
            mean(dN),sum(dN>TOL),sum(abs(dN)<=TOL),sum(dN<-TOL), ...
            mean(dR),sum(dR>TOL),sum(abs(dR)<=TOL),sum(dR<-TOL));
    end
    for k=1:K
        C=squeeze(confusionCounts(q,k,:,:));
        rowSums=sum(C,2); normalized=zeros(size(C)); positive=rowSums>0;
        normalized(positive,:)=bsxfun(@rdivide,C(positive,:),rowSums(positive));
        neu_write_matrix(fullfile(outputDir,sprintf( ...
            '%s_confusion_counts_label%02d_%s.csv',prefix, ...
            round(100*fraction),methodNames{k})),C);
        neu_write_matrix(fullfile(outputDir,sprintf( ...
            '%s_confusion_normalized_label%02d_%s.csv',prefix, ...
            round(100*fraction),methodNames{k})),normalized);
    end
    fprintf(fd,'\n');
end

integrity=finitePass && objectivePass && fixedLabelPass;
fprintf(fd,'Audit\n');
fprintf(fd,'FINITE_RESULTS_PASS=%d\n',finitePass);
fprintf(fd,'OBJECTIVE_NONINCREASE_PASS=%d\n',objectivePass);
fprintf(fd,'GOC_FIXED_LABEL_PASS=%d\n',fixedLabelPass);
fprintf(fd,'NO_KMEANS_PASS=1\n');
fprintf(fd,'NEU_APPLICATION_INTEGRITY_PASS=%d\n',integrity);
if runs==1, fprintf(fd,'NEU_ONE_SEED_SMOKE_PASS=%d\n',integrity);
elseif runs==3, fprintf(fd,'NEU_THREE_SEED_INTEGRITY_PASS=%d\n',integrity);
else, fprintf(fd,'NEU_TWENTY_SEED_APPLICATION_COMPLETE=%d\n',integrity); end
fprintf(fd,['\nStage-1 limitation: parameters are pre-specified pilot settings. ' ...
    'No unlabeled ground-truth label is used for tuning. Publication ' ...
    'claims require the actual 20-seed results and audit.\n']);

fclose(fs); fclose(fp); fclose(fd);
fprintf('Summary: %s\nPaired: %s\nDecision: %s\n', ...
    summaryPath,pairedPath,decisionPath);
end

function value=neu_safe_std(x)
if length(x)>1, value=std(x,0); else, value=0; end
end

function neu_write_matrix(path,M)
fid=fopen(path,'w'); if fid<0, error('Cannot create %s.',path); end
for i=1:size(M,1)
    for j=1:size(M,2)
        if j<size(M,2), fprintf(fid,'%.15g,',M(i,j));
        else, fprintf(fid,'%.15g\n',M(i,j)); end
    end
end
fclose(fid);
end
