function summarize_GOCNMF_terminal_validation_R2009a
% Aggregate all six completed dataset result files.
%
% Main output:
%   GOCNMF_terminal_six_dataset_summary.txt

close all;
clc;

config = GOCNMF_terminal_validation_config_R2009a;

num_datasets = length(config.dataset_names);
num_seeds = length(config.seeds);

% Dataset summary columns:
%  1 base ACC, 2 base NMI, 3 clean ACC, 4 clean NMI,
%  5 dALL ACC, 6 dALL NMI, 7 relative error reduction,
%  8 all-control specificity ACC, 9 all-control specificity NMI,
% 10 all dual wins, 11 rescue, 12 harm,
% 13 dLU ACC, 14 dLU NMI, 15 dUU ACC, 16 dUU NMI,
% 17 interaction ACC, 18 interaction NMI,
% 19 phiLU ACC, 20 phiLU NMI, 21 phiUU ACC, 22 phiUU NMI,
% 23 LU share ACC, 24 LU share NMI,
% 25 LU causal support, 26 UU causal support,
% 27 LU AUC, 28 LU lift, 29 LU recall, 30 LU detect support,
% 31 UU AUC, 32 UU lift, 33 UU recall, 34 UU detect support,
% 35 purity, 36 monotone.
summary_values = NaN*ones(num_datasets,36);

for dataset_id = 1:num_datasets
    dataset_name = config.dataset_names{dataset_id};
    result_file = fullfile( ...
        ['terminal_results_' dataset_name], ...
        ['terminal_results_' dataset_name '.mat']);

    if exist(result_file,'file')~=2
        error('Missing result file: %s.',result_file);
    end

    data = load(result_file);

    if ~all(data.completed)
        error('Dataset %s is incomplete: %d/%d splits.', ...
            dataset_name,sum(data.completed),num_seeds);
    end

    scores = data.condition_scores;
    edges = data.edge_stats;
    reliability = data.reliability_metrics;
    rescue_harm = data.rescue_harm;

    idx_BASE = find_condition_R2009a( ...
        config.condition_names,'BASE');
    idx_LU = find_condition_R2009a( ...
        config.condition_names,'REMOVE_WRONG_LU');
    idx_UU = find_condition_R2009a( ...
        config.condition_names,'REMOVE_WRONG_UU');
    idx_ALL = find_condition_R2009a( ...
        config.condition_names,'REMOVE_WRONG_ALL');

    idx_CTRL_LU = zeros(3,1);
    idx_CTRL_UU = zeros(3,1);
    idx_CTRL_ALL = zeros(3,1);

    for control_id = 1:3
        idx_CTRL_LU(control_id) = find_condition_R2009a( ...
            config.condition_names, ...
            ['CTRL_CORRECT_LU_' num2str(control_id)]);
        idx_CTRL_UU(control_id) = find_condition_R2009a( ...
            config.condition_names, ...
            ['CTRL_CORRECT_UU_' num2str(control_id)]);
        idx_CTRL_ALL(control_id) = find_condition_R2009a( ...
            config.condition_names, ...
            ['CTRL_CORRECT_ALL_' num2str(control_id)]);
    end

    BASE = squeeze(scores(:,idx_BASE,1:2));
    LU = squeeze(scores(:,idx_LU,1:2));
    UU = squeeze(scores(:,idx_UU,1:2));
    ALL = squeeze(scores(:,idx_ALL,1:2));

    CTRL_LU = zeros(num_seeds,2);
    CTRL_UU = zeros(num_seeds,2);
    CTRL_ALL = zeros(num_seeds,2);

    for control_id = 1:3
        CTRL_LU = CTRL_LU+squeeze( ...
            scores(:,idx_CTRL_LU(control_id),1:2))/3;
        CTRL_UU = CTRL_UU+squeeze( ...
            scores(:,idx_CTRL_UU(control_id),1:2))/3;
        CTRL_ALL = CTRL_ALL+squeeze( ...
            scores(:,idx_CTRL_ALL(control_id),1:2))/3;
    end

    dLU = LU-BASE;
    dUU = UU-BASE;
    dALL = ALL-BASE;
    interaction = dALL-dLU-dUU;

    phiLU = 0.5*(dLU+dALL-dUU);
    phiUU = 0.5*(dUU+dALL-dLU);

    mean_dALL = mean(dALL,1);
    relative_error_reduction = mean( ...
        dALL(:,1)./max(1-BASE(:,1),1e-12*ones(num_seeds,1)));

    specificity_ALL = mean(ALL-CTRL_ALL,1);
    specificity_LU = mean(LU-CTRL_LU,1);
    specificity_UU = mean(UU-CTRL_UU,1);

    all_dual_wins = sum(dALL(:,1)>0 & dALL(:,2)>0);
    LU_dual_wins = sum(dLU(:,1)>0 & dLU(:,2)>0);
    UU_dual_wins = sum(dUU(:,1)>0 & dUU(:,2)>0);

    LU_causal_support = ...
        all(mean(dLU,1)>0) && ...
        all(specificity_LU>=0.002) && ...
        LU_dual_wins>=14;

    UU_causal_support = ...
        all(mean(dUU,1)>0) && ...
        all(specificity_UU>=0.002) && ...
        UU_dual_wins>=14;

    mean_phiLU = mean(phiLU,1);
    mean_phiUU = mean(phiUU,1);
    LU_share = mean_phiLU./max(abs(mean_dALL),1e-12*ones(1,2));

    LU_values = squeeze(reliability(:,1,:));
    UU_values = squeeze(reliability(:,2,:));

    LU_AUC = mean(LU_values(:,1));
    LU_lift = mean(LU_values(:,4));
    LU_recall = mean(LU_values(:,3));
    LU_stable = sum(LU_values(:,1)>=0.65);

    UU_AUC = mean(UU_values(:,1));
    UU_lift = mean(UU_values(:,4));
    UU_recall = mean(UU_values(:,3));
    UU_stable = sum(UU_values(:,1)>=0.65);

    LU_detect_support = ...
        LU_AUC>=0.70 && ...
        LU_stable>=14 && ...
        LU_lift>=2.0 && ...
        LU_recall>=0.40 && ...
        min(LU_values(:,9))>=5 && ...
        min(LU_values(:,10))>=5;

    UU_detect_support = ...
        UU_AUC>=0.70 && ...
        UU_stable>=14 && ...
        UU_lift>=2.0 && ...
        UU_recall>=0.40 && ...
        min(UU_values(:,9))>=5 && ...
        min(UU_values(:,10))>=5;

    all_increases = squeeze(scores(:,:,4));
    monotone = max(all_increases(:))==0;

    summary_values(dataset_id,:) = [ ...
        mean(BASE,1),mean(ALL,1),mean_dALL, ...
        relative_error_reduction,specificity_ALL, ...
        all_dual_wins,mean(rescue_harm(:,1)), ...
        mean(rescue_harm(:,2)), ...
        mean(dLU,1),mean(dUU,1),mean(interaction,1), ...
        mean_phiLU,mean_phiUU,LU_share, ...
        LU_causal_support,UU_causal_support, ...
        LU_AUC,LU_lift,LU_recall,LU_detect_support, ...
        UU_AUC,UU_lift,UU_recall,UU_detect_support, ...
        mean(edges(:,8)),monotone];
end

% -------------------------------------------------------------------------
% Locked terminal decision
% -------------------------------------------------------------------------
oracle_support = false(num_datasets,1);
joint_causal = false(num_datasets,1);
UU_dominant = false(num_datasets,1);
LU_detect = false(num_datasets,1);
UU_detect = false(num_datasets,1);

for d = 1:num_datasets
    oracle_support(d) = ...
        summary_values(d,5)>0 && ...
        summary_values(d,6)>0 && ...
        summary_values(d,7)>=0.20 && ...
        summary_values(d,8)>=0.005 && ...
        summary_values(d,9)>=0.005 && ...
        summary_values(d,10)>=16 && ...
        summary_values(d,11)>=0.20 && ...
        summary_values(d,12)<=0.05;

    joint_causal(d) = ...
        summary_values(d,25)>=0.5 && ...
        summary_values(d,26)>=0.5;

    UU_dominant(d) = ...
        summary_values(d,21)>summary_values(d,19) && ...
        summary_values(d,22)>summary_values(d,20);

    LU_detect(d) = summary_values(d,30)>=0.5;
    UU_detect(d) = summary_values(d,34)>=0.5;
end

no_oracle_harm = ...
    min(summary_values(:,5))>=-0.002 && ...
    min(summary_values(:,6))>=-0.002;

all_monotone = all(summary_values(:,36)>=0.5);

overall_pass = ...
    sum(oracle_support)>=5 && ...
    sum(joint_causal)>=4 && ...
    sum(UU_dominant)>=4 && ...
    sum(LU_detect)>=4 && ...
    sum(UU_detect)>=1 && ...
    sum(UU_detect)<=3 && ...
    no_oracle_harm && ...
    all_monotone;

summary_file = 'GOCNMF_terminal_six_dataset_summary.txt';
fid = fopen(summary_file,'w');

if fid<0
    error('Cannot create terminal summary.');
end

fprintf(fid,'GOCNMF SIX-DATASET TERMINAL VALIDATION\n');
fprintf(fid,'======================================\n\n');

fprintf(fid,'LOCKED PROTOCOL\n');
fprintf(fid,'---------------\n');
fprintf(fid,'Datasets: COIL20, COIL100, MNIST, Optdigits, PIE, YaleB.\n');
fprintf(fid,'Splits: 20 seeds from 20260617 to 20260636.\n');
fprintf(fid,'Labels: 10 percent independently in each class.\n');
fprintf(fid,'GOCNMF iterations: 50.\n');
fprintf(fid,'Published p and alpha values are fixed.\n');
fprintf(fid,'Three matched correct-edge controls per edge type.\n');
fprintf(fid,'The edge-excluded reliability score is unchanged.\n\n');

for d = 1:num_datasets
    v = summary_values(d,:);

    fprintf(fid,'Dataset: %s\n',config.dataset_names{d});
    fprintf(fid,'  BASE ACC/NMI: %.8f / %.8f\n',v(1),v(2));
    fprintf(fid,'  ORACLE-CLEAN ACC/NMI: %.8f / %.8f\n',v(3),v(4));
    fprintf(fid,'  Oracle dACC/dNMI: %+.8f / %+.8f\n',v(5),v(6));
    fprintf(fid,'  Relative baseline-error reduction: %.8f\n',v(7));
    fprintf(fid,'  Oracle specificity ACC/NMI: %+.8f / %+.8f\n', ...
        v(8),v(9));
    fprintf(fid,'  Oracle dual wins: %d / %d\n',round(v(10)),num_seeds);
    fprintf(fid,'  Rescue/harm rate: %.8f / %.8f\n',v(11),v(12));
    fprintf(fid,'  dLU ACC/NMI: %+.8f / %+.8f\n',v(13),v(14));
    fprintf(fid,'  dUU ACC/NMI: %+.8f / %+.8f\n',v(15),v(16));
    fprintf(fid,'  Interaction ACC/NMI: %+.8f / %+.8f\n',v(17),v(18));
    fprintf(fid,'  Shapley LU ACC/NMI: %+.8f / %+.8f\n',v(19),v(20));
    fprintf(fid,'  Shapley UU ACC/NMI: %+.8f / %+.8f\n',v(21),v(22));
    fprintf(fid,'  LU Shapley share ACC/NMI: %.8f / %.8f\n',v(23),v(24));
    fprintf(fid,'  LU causal support: %s\n', ...
        pass_text_R2009a(v(25)>=0.5));
    fprintf(fid,'  UU causal support: %s\n', ...
        pass_text_R2009a(v(26)>=0.5));
    fprintf(fid,'  LU detect AUC/lift/recall: %.8f / %.8f / %.8f\n', ...
        v(27),v(28),v(29));
    fprintf(fid,'  LU detect support: %s\n', ...
        pass_text_R2009a(v(30)>=0.5));
    fprintf(fid,'  UU detect AUC/lift/recall: %.8f / %.8f / %.8f\n', ...
        v(31),v(32),v(33));
    fprintf(fid,'  UU detect support: %s\n', ...
        pass_text_R2009a(v(34)>=0.5));
    fprintf(fid,'  Original graph purity: %.8f\n',v(35));
    fprintf(fid,'  Objective histories nonincreasing: %s\n', ...
        pass_text_R2009a(v(36)>=0.5));
    fprintf(fid,'  Oracle dataset support: %s\n', ...
        pass_text_R2009a(oracle_support(d)));
    fprintf(fid,'  Joint LU-UU causal support: %s\n', ...
        pass_text_R2009a(joint_causal(d)));
    fprintf(fid,'  UU Shapley dominant: %s\n\n', ...
        pass_text_R2009a(UU_dominant(d)));
end

fprintf(fid,'LOCKED TERMINAL DECISION\n');
fprintf(fid,'------------------------\n');
fprintf(fid,'Oracle support datasets: %d / %d\n', ...
    sum(oracle_support),num_datasets);
fprintf(fid,'Joint LU-UU causal datasets: %d / %d\n', ...
    sum(joint_causal),num_datasets);
fprintf(fid,'UU Shapley-dominant datasets: %d / %d\n', ...
    sum(UU_dominant),num_datasets);
fprintf(fid,'LU detectable datasets: %d / %d\n', ...
    sum(LU_detect),num_datasets);
fprintf(fid,'UU detectable datasets: %d / %d\n', ...
    sum(UU_detect),num_datasets);
fprintf(fid,'No oracle harm beyond 0.002: %s\n', ...
    pass_text_R2009a(no_oracle_harm));
fprintf(fid,'All objectives nonincreasing: %s\n', ...
    pass_text_R2009a(all_monotone));

fprintf(fid,'\nOVERALL TERMINAL PATTERN: %s\n', ...
    pass_text_R2009a(overall_pass));

if overall_pass
    fprintf(fid,['Interpretation: the six-dataset evidence supports the ' ...
        'research conclusion that wrong graph edges are a major causal ' ...
        'bottleneck, UU errors carry most aggregate damage, LU errors are ' ...
        'more observable, and UU detectability is heterogeneous and often ' ...
        'insufficient. Proceed to the causal-mechanism and detectability-' ...
        'boundary paper. Do not construct a weighted GOCNMF model.\n']);
else
    fprintf(fid,['Interpretation: the predeclared six-dataset pattern does ' ...
        'not replicate. Stop the GOCNMF improvement project rather than ' ...
        'changing thresholds, scores, or models after observing results.\n']);
end

fclose(fid);

save('GOCNMF_terminal_six_dataset_results.mat', ...
    'summary_values','oracle_support','joint_causal', ...
    'UU_dominant','LU_detect','UU_detect','overall_pass', ...
    'no_oracle_harm','all_monotone','config');

plot_terminal_results_R2009a( ...
    config.dataset_names,summary_values);

fprintf('\nTerminal aggregation completed.\n');
fprintf('Summary: %s\n',summary_file);
fprintf('Decision: %s\n',pass_text_R2009a(overall_pass));

end

function index = find_condition_R2009a(names,target)

index = 0;

for i = 1:length(names)
    if strcmp(names{i},target)
        index = i;
        return;
    end
end

if index==0
    error('Condition %s is missing.',target);
end

end

function text = pass_text_R2009a(flag)

if flag
    text = 'PASS';
else
    text = 'FAIL';
end

end

function plot_terminal_results_R2009a(names,values)

x = 1:length(names);

figure;
bar(x,[values(:,1),values(:,3)]);
set(gca,'XTick',x);
set(gca,'XTickLabel',names);
ylabel('ACC');
legend('Original GOCNMF','Oracle-clean graph','Location','Best');
title('Six-dataset oracle graph intervention: ACC');
grid on;
box on;
saveas(gcf,'terminal_oracle_ACC.png');

figure;
bar(x,[values(:,2),values(:,4)]);
set(gca,'XTick',x);
set(gca,'XTickLabel',names);
ylabel('NMI');
legend('Original GOCNMF','Oracle-clean graph','Location','Best');
title('Six-dataset oracle graph intervention: NMI');
grid on;
box on;
saveas(gcf,'terminal_oracle_NMI.png');

figure;
bar(x,[values(:,19),values(:,21)]);
set(gca,'XTick',x);
set(gca,'XTickLabel',names);
ylabel('Shapley contribution to ACC gain');
legend('LU','UU','Location','Best');
title('LU-UU causal contribution decomposition');
grid on;
box on;
saveas(gcf,'terminal_Shapley_ACC.png');

figure;
bar(x,[values(:,27),values(:,31)]);
set(gca,'XTick',x);
set(gca,'XTickLabel',names);
ylabel('Correct-edge detection AUC');
legend('LU','UU','Location','Best');
title('Frozen edge-reliability detectability');
grid on;
box on;
saveas(gcf,'terminal_reliability_AUC.png');

end
