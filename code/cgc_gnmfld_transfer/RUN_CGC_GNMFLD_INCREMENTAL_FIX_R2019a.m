% RUN_CGC_GNMFLD_INCREMENTAL_FIX_R2019a
% Recomputes ONLY:
%   harmonic calibration -> theta -> CGC-GNMFLD
% and reuses the already completed GNMFLD baseline from the previous raw CSV.
%
% This is the correct incremental rerun after fixing Chinese-MATLAB detection
% of nearly-singular harmonic systems.

clear;
clc;

data_dir = fullfile(pwd, 'data');
old_raw  = fullfile(pwd, 'previous_results', 'CGC_GNMFLD_transfer_raw.csv');
out_dir  = fullfile(pwd, 'results_incremental_fixed');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
if ~exist(old_raw,'file')
    error('Previous raw CSV not found: %s', old_raw);
end

dataset_names = {'PIE','YaleB','COIL20','COIL100','Optdigits','MNIST'};
p_list     = [5, 5, 5, 5, 5, 5];
alpha_list = [1e5, 1e5, 1e4, 1e4, 1e5, 1e5];
beta_list  = [10, 100, 10, 10, 10, 10];

formal_seeds = 20260617:20260636;
label_rate = 0.10;
F = 5;
kappa = 1;
max_iter = 200;
eps_update = 1e-12;
harmonic_stab = 1e-10;
knn_block_size = 512;

oldT = readtable(old_raw);

new_raw = fullfile(out_dir, 'CGC_GNMFLD_transfer_raw_FIXED.csv');
fid = fopen(new_raw,'w');
if fid < 0
    error('Cannot create %s',new_raw);
end
fprintf(fid, ['Dataset,Seed,Theta,Delta,S,Fallback,StabW0,StabW1,' ...
              'ACC_GNMFLD,NMI_GNMFLD,ACC_CGC_GNMFLD,NMI_CGC_GNMFLD,' ...
              'dACC,dNMI\n']);

all_results = cell(numel(dataset_names),1);

for dset_id = 1:numel(dataset_names)
    dataset_name = dataset_names{dset_id};

    fprintf('\n============================================================\n');
    fprintf('Incremental harmonic fix: %s\n',dataset_name);
    fprintf('============================================================\n');

    data_file = fullfile(data_dir,[dataset_name '.mat']);
    [X,y] = load_dataset_cgc_transfer(data_file);
    X = normalize_columns_l2(X);

    n = size(X,2);
    m = size(X,1);
    c = max(y);

    p = p_list(dset_id);
    alpha = alpha_list(dset_id);
    beta = beta_list(dset_id);

    fprintf('Constructing W0/W1 once for %s ...\n',dataset_name);
    [W0,W1] = build_cgc_graph_pair(X,p,knn_block_size,eps_update);

    rows_old = strcmp(oldT.Dataset,dataset_name);
    Told = oldT(rows_old,:);
    if height(Told) ~= numel(formal_seeds)
        fclose(fid);
        error('Expected 20 old rows for %s, found %d.',dataset_name,height(Told));
    end

    R = zeros(numel(formal_seeds),14);

    for s_id = 1:numel(formal_seeds)
        seed = formal_seeds(s_id);

        rr = find(Told.Seed == seed,1);
        if isempty(rr)
            fclose(fid);
            error('Old baseline row missing: %s seed %d.',dataset_name,seed);
        end

        % Recreate the exact formal split.
        [labeled_idx,unlabeled_idx] = make_labeled_split(y,label_rate,seed);

        % Recreate the exact paired initialization used previously.
        rng(seed + 100000,'twister');
        U0 = rand(m,c) + 0.1;
        V0 = rand(n,c) + 0.1;

        % Recompute calibration with robust warning-ID detection and 1e-10I fallback.
        [theta,Delta,svar,diaginfo] = cgc_oof_theta( ...
            W0,W1,labeled_idx,y,c,F,kappa,harmonic_stab,seed+200000);

        Wtheta = (1-theta)*W0 + theta*W1;

        % Reuse completed GNMFLD baseline from old raw CSV.
        acc_base = Told.ACC_GNMFLD(rr) / 100;
        nmi_base = Told.NMI_GNMFLD(rr) / 100;

        % Recompute only CGC-GNMFLD.
        V_cgc = gnmfld_fit(X,y,labeled_idx,c,Wtheta, ...
            alpha,beta,U0,V0,max_iter,eps_update);

        [~,pred_cgc] = max(V_cgc,[],2);
        [acc_cgc,nmi_cgc] = direct_acc_nmi( ...
            y(unlabeled_idx),pred_cgc(unlabeled_idx),c);

        dacc = acc_cgc - acc_base;
        dnmi = nmi_cgc - nmi_base;
        fallback = double(theta == 0);

        R(s_id,:) = [seed,theta,Delta,svar,fallback, ...
            diaginfo.stab0,diaginfo.stab1, ...
            acc_base,nmi_base,acc_cgc,nmi_cgc,dacc,dnmi,numel(unlabeled_idx)];

        fprintf(fid,'%s,%d,%.16g,%.16g,%.16g,%d,%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
            dataset_name,seed,theta,Delta,svar,fallback, ...
            diaginfo.stab0,diaginfo.stab1, ...
            100*acc_base,100*nmi_base,100*acc_cgc,100*nmi_cgc, ...
            100*dacc,100*dnmi);

        fprintf(['  seed=%d theta=%.4f stab=%d/%d' ...
                 ' ACC %+.3f pp NMI %+.3f pp\n'], ...
            seed,theta,diaginfo.stab0,diaginfo.stab1,100*dacc,100*dnmi);

        clear labeled_idx unlabeled_idx U0 V0 Wtheta V_cgc pred_cgc
    end

    all_results{dset_id} = R;
end

fclose(fid);

write_transfer_summary_fixed(dataset_names,all_results,out_dir);

fprintf('\nIncremental correction complete.\n');
fprintf('Upload the three files in:\n  %s\n',out_dir);
