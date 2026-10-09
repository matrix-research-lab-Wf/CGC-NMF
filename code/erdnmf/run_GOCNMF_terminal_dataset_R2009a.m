function run_GOCNMF_terminal_dataset_R2009a(dataset_identifier)
% RUN_GOCNMF_TERMINAL_DATASET_R2009A
%
% Six-dataset terminal validation worker.
%
% Examples:
%   run_GOCNMF_terminal_dataset_R2009a(1)          % COIL20
%   run_GOCNMF_terminal_dataset_R2009a('YaleB')
%
% The worker is resumable and saves one checkpoint after every split.
%
% It performs:
%   - original GOCNMF;
%   - oracle removal of wrong LU, wrong UU, and all active wrong edges;
%   - three matched correct-edge controls for LU, UU, and all active edges;
%   - LU/UU Shapley causal decomposition;
%   - the frozen edge-excluded reliability audit.
%
% It does not train a new weighted model.

close all;
clc;

config = GOCNMF_terminal_validation_config_R2009a;
dataset_id = resolve_dataset_id_R2009a( ...
    dataset_identifier,config.dataset_names);

dataset_name = config.dataset_names{dataset_id};
aliases = config.dataset_aliases{dataset_id};
p = config.p_values(dataset_id);
alpha = config.alpha_values(dataset_id);

source_file = locate_dataset_file_R2009a(aliases,dataset_name);

fprintf('Dataset source: %s\n',source_file);

[X,y] = load_dataset_R2009a(source_file);
[X,y] = preprocess_dataset_R2009a(X,y,config.eps0);

[m,n] = size(X);
c = length(unique(y));
num_seeds = length(config.seeds);
num_conditions = length(config.condition_names);

fprintf('Dataset %s: %d features x %d samples x %d classes\n', ...
    dataset_name,m,n,c);
fprintf('p=%d, alpha=%.16g, splits=%d\n', ...
    p,alpha,num_seeds);

result_folder = ['terminal_results_' dataset_name];

if exist(result_folder,'dir')~=7
    mkdir(result_folder);
end

checkpoint_file = fullfile(result_folder, ...
    ['terminal_checkpoint_' dataset_name '.mat']);

% condition_scores(seed,condition,metric)
% metric: 1 ACC, 2 NMI, 3 RE, 4 objective increases, 5 runtime.
condition_scores = NaN*ones(num_seeds,num_conditions,5);

% edge_stats:
% 1 total, 2 wrong_LL, 3 wrong_LU, 4 wrong_UU,
% 5 correct_LU, 6 correct_UU, 7 correct_active, 8 purity.
edge_stats = NaN*ones(num_seeds,8);

% reliability_metrics(seed,type,metric), type 1 LU, 2 UU.
% metric:
% 1 AUC, 2 low-score wrong precision, 3 wrong recall,
% 4 enrichment, 5 prevalence, 6 selected fraction,
% 7 zero fraction, 8 edge count, 9 wrong count, 10 correct count.
reliability_metrics = NaN*ones(num_seeds,2,10);

% rescue_harm(seed,1:2): rescue and harm after all-wrong removal.
rescue_harm = NaN*ones(num_seeds,2);

completed = false(num_seeds,1);

if exist(checkpoint_file,'file')==2
    old = load(checkpoint_file);

    required = {'condition_scores','edge_stats', ...
        'reliability_metrics','rescue_harm','completed'};

    valid_checkpoint = true;

    for r = 1:length(required)
        if ~isfield(old,required{r})
            valid_checkpoint = false;
        end
    end

    if valid_checkpoint
        condition_scores = old.condition_scores;
        edge_stats = old.edge_stats;
        reliability_metrics = old.reliability_metrics;
        rescue_harm = old.rescue_harm;
        completed = old.completed;
        fprintf('Resuming checkpoint: %d/%d splits completed.\n', ...
            sum(completed),num_seeds);
    else
        error('Existing checkpoint is incomplete or incompatible.');
    end
end

raw_file = fullfile(result_folder, ...
    ['terminal_raw_' dataset_name '.csv']);

for seed_id = 1:num_seeds
    if completed(seed_id)
        fprintf('Skipping completed seed %d.\n', ...
            config.seeds(seed_id));
        continue;
    end

    seed = config.seeds(seed_id);

    fprintf('\n----------------------------------------------------\n');
    fprintf('Dataset %s, split %d/%d, seed %d\n', ...
        dataset_name,seed_id,num_seeds,seed);
    fprintf('----------------------------------------------------\n');

    labeled = make_labeled_split_R2009a( ...
        y,c,config.label_fraction,seed);

    labeled_mask = false(n,1);
    labeled_mask(labeled) = true;
    unlabeled = find(~labeled_mask);

    fprintf('Building exact blockwise graph...\n');
    W = build_binary_graph_blockwise_R2009a( ...
        X,p,config.graph_block_size);

    [edge_row,edge_col,edge_type,edge_wrong,stats] = ...
        classify_edges_R2009a(W,y,labeled_mask);

    edge_stats(seed_id,:) = stats;

    wrong_LU = find(edge_wrong & edge_type==2);
    wrong_UU = find(edge_wrong & edge_type==3);
    wrong_active = [wrong_LU;wrong_UU];

    correct_LU = find(~edge_wrong & edge_type==2);
    correct_UU = find(~edge_wrong & edge_type==3);
    correct_active = [correct_LU;correct_UU];

    fprintf(['Edges total=%d, wrong LL/LU/UU=%d/%d/%d, ' ...
        'purity=%.6f\n'], ...
        round(stats(1)),round(stats(2)),round(stats(3)), ...
        round(stats(4)),stats(8));

    [U0,V0] = initialize_GOCNMF_R2009a( ...
        X,y,labeled,c,seed+100000*dataset_id,config.eps0);

    graphs = cell(num_conditions,1);
    graphs{1} = W;
    graphs{2} = remove_edge_ids_R2009a( ...
        W,edge_row,edge_col,wrong_LU);
    graphs{3} = remove_edge_ids_R2009a( ...
        W,edge_row,edge_col,wrong_UU);
    graphs{4} = remove_edge_ids_R2009a( ...
        W,edge_row,edge_col,wrong_active);

    for control_id = 1:config.num_controls
        control_seed = seed+1000000*dataset_id ...
            +10000*control_id;

        chosen = choose_matched_edges_R2009a( ...
            correct_LU,length(wrong_LU),control_seed+11);
        graphs{4+control_id} = remove_edge_ids_R2009a( ...
            W,edge_row,edge_col,chosen);

        chosen = choose_matched_edges_R2009a( ...
            correct_UU,length(wrong_UU),control_seed+22);
        graphs{7+control_id} = remove_edge_ids_R2009a( ...
            W,edge_row,edge_col,chosen);

        chosen = choose_matched_edges_R2009a( ...
            correct_active,length(wrong_active),control_seed+33);
        graphs{10+control_id} = remove_edge_ids_R2009a( ...
            W,edge_row,edge_col,chosen);
    end

    predictions = cell(num_conditions,1);
    baseline_U = [];
    baseline_V = [];
    baseline_history = [];

    for condition_id = 1:num_conditions
        fprintf('  %-22s ',config.condition_names{condition_id});

        start_time = tic;

        [U,V,prediction,history,reconstruction_norm] = ...
            train_GOCNMF_fixed_initial_R2009a( ...
            X,y,labeled,c,graphs{condition_id},alpha, ...
            config.max_iterations,U0,V0,config.eps0);

        runtime_value = toc(start_time);

        ACC = mean(prediction(unlabeled)==y(unlabeled));
        NMI = nmi_score_R2009a( ...
            y(unlabeled),prediction(unlabeled));
        RE = reconstruction_norm/max(norm(X,'fro'),config.eps0);
        increases = count_increases_R2009a( ...
            history,config.objective_tolerance);

        condition_scores(seed_id,condition_id,:) = ...
            [ACC,NMI,RE,increases,runtime_value];

        predictions{condition_id} = prediction;

        if condition_id==1
            baseline_U = U;
            baseline_V = V;
            baseline_history = history;
        end

        fprintf('ACC %.6f, NMI %.6f\n',ACC,NMI);
    end

    base_correct = predictions{1}(unlabeled)==y(unlabeled);
    clean_correct = predictions{4}(unlabeled)==y(unlabeled);

    rescue_count = sum(~base_correct & clean_correct);
    harm_count = sum(base_correct & ~clean_correct);

    rescue_harm(seed_id,1) = ...
        rescue_count/max(sum(~base_correct),1);
    rescue_harm(seed_id,2) = ...
        harm_count/max(sum(base_correct),1);

    [LU_score,LU_wrong,UU_score,UU_wrong] = ...
        compute_edge_excluded_scores_R2009a( ...
        X,baseline_U,baseline_V,W,y,labeled_mask, ...
        alpha,c,config.eps0);

    LU_values = evaluate_score_R2009a( ...
        LU_score,LU_wrong,config.low_score_fraction);
    UU_values = evaluate_score_R2009a( ...
        UU_score,UU_wrong,config.low_score_fraction);

    reliability_metrics(seed_id,1,:) = LU_values;
    reliability_metrics(seed_id,2,:) = UU_values;

    fprintf(['Reliability LU AUC %.4f, lift %.4f, recall %.4f; ' ...
        'UU AUC %.4f, lift %.4f, recall %.4f\n'], ...
        LU_values(1),LU_values(4),LU_values(3), ...
        UU_values(1),UU_values(4),UU_values(3));

    completed(seed_id) = true;

    save(checkpoint_file, ...
        'dataset_name','source_file','p','alpha','config', ...
        'condition_scores','edge_stats','reliability_metrics', ...
        'rescue_harm','completed','baseline_history');

    write_dataset_raw_R2009a( ...
        raw_file,dataset_name,config.seeds, ...
        config.condition_names,condition_scores,edge_stats, ...
        reliability_metrics,rescue_harm,completed);

    drawnow;
end

result_file = fullfile(result_folder, ...
    ['terminal_results_' dataset_name '.mat']);

save(result_file, ...
    'dataset_name','source_file','p','alpha','config', ...
    'condition_scores','edge_stats','reliability_metrics', ...
    'rescue_harm','completed');

write_dataset_raw_R2009a( ...
    raw_file,dataset_name,config.seeds, ...
    config.condition_names,condition_scores,edge_stats, ...
    reliability_metrics,rescue_harm,completed);

fprintf('\nDataset %s completed: %d/%d splits.\n', ...
    dataset_name,sum(completed),num_seeds);
fprintf('Result file: %s\n',result_file);

end

% =========================================================================
% Configuration and file loading
% =========================================================================

function dataset_id = resolve_dataset_id_R2009a( ...
    identifier,dataset_names)

if isnumeric(identifier)
    dataset_id = round(identifier(1));

    if dataset_id<1 || dataset_id>length(dataset_names)
        error('Dataset index must be from 1 to %d.', ...
            length(dataset_names));
    end

elseif ischar(identifier)
    dataset_id = 0;

    for i = 1:length(dataset_names)
        if strcmpi(identifier,dataset_names{i})
            dataset_id = i;
            break;
        end
    end

    if dataset_id==0
        error('Unknown dataset name: %s.',identifier);
    end
else
    error('Dataset identifier must be an index or character name.');
end

end

function file_name = locate_dataset_file_R2009a( ...
    aliases,dataset_name)

current_folder = pwd;
parent_folder = fileparts(current_folder);

file_name = recursive_find_R2009a( ...
    current_folder,aliases,4);

if isempty(file_name) && ...
        ~isempty(parent_folder) && ...
        ~strcmp(parent_folder,current_folder)

    file_name = recursive_find_R2009a( ...
        parent_folder,aliases,3);
end

if isempty(file_name)
    filter_text = '*.mat';

    [selected,path_name] = uigetfile( ...
        filter_text,['Select the full ' dataset_name ' MAT file']);

    if isequal(selected,0)
        error('No MAT file selected for %s.',dataset_name);
    end

    file_name = fullfile(path_name,selected);
end

end

function found_file = recursive_find_R2009a( ...
    folder,target_names,remaining_depth)

found_file = '';

for i = 1:length(target_names)
    candidate = fullfile(folder,target_names{i});

    if exist(candidate,'file')==2
        found_file = candidate;
        return;
    end
end

if remaining_depth<=0
    return;
end

entries = dir(folder);

for i = 1:length(entries)
    name = entries(i).name;

    if entries(i).isdir && ...
            ~strcmp(name,'.') && ~strcmp(name,'..')

        child = fullfile(folder,name);

        found_file = recursive_find_R2009a( ...
            child,target_names,remaining_depth-1);

        if ~isempty(found_file)
            return;
        end
    end
end

end

function [X,y] = load_dataset_R2009a(file_name)

data = load(file_name);
field_names = fieldnames(data);

feature_priority = {'fea','X','data','features','trainData'};
label_priority = {'gnd','labels','label','y','Y','truth'};

X = [];
y = [];

for i = 1:length(feature_priority)
    if isfield(data,feature_priority{i})
        value = data.(feature_priority{i});

        if isnumeric(value) && ndims(value)==2 && ...
                min(size(value))>1
            X = value;
            break;
        end
    end
end

for i = 1:length(label_priority)
    if isfield(data,label_priority{i})
        value = data.(label_priority{i});

        if isnumeric(value) && isvector(value)
            y = value(:);
            break;
        end
    end
end

if isempty(y)
    for i = 1:length(field_names)
        value = data.(field_names{i});

        if isnumeric(value) && isvector(value) && ...
                length(value)>=10
            y = value(:);
            break;
        end
    end
end

if isempty(X) && ~isempty(y)
    best_elements = 0;

    for i = 1:length(field_names)
        value = data.(field_names{i});

        if isnumeric(value) && ndims(value)==2 && ...
                min(size(value))>1 && ...
                (size(value,1)==length(y) || ...
                 size(value,2)==length(y))

            elements = numel(value);

            if elements>best_elements
                X = value;
                best_elements = elements;
            end
        end
    end
end

if isempty(X) || isempty(y)
    fprintf('Fields found in %s:\n',file_name);

    for i = 1:length(field_names)
        value = data.(field_names{i});
        dimensions = size(value);
        fprintf('  %s: ',field_names{i});
        fprintf('%dx',dimensions(1:end-1));
        fprintf('%d\n',dimensions(end));
    end

    error('Could not identify feature and label arrays.');
end

X = double(X);
y = double(y(:));

if size(X,2)==length(y)
    % Features by samples.
elseif size(X,1)==length(y)
    X = X';
else
    error('Feature dimensions do not match label count.');
end

end

function [X,y] = preprocess_dataset_R2009a(X,y,eps0)

valid = isfinite(y) & all(isfinite(X),1)';
X = X(:,valid);
y = y(valid);

minimum_value = min(X(:));

if minimum_value<0
    X = X-minimum_value;
end

X(X<0) = 0;

norms = sqrt(sum(X.^2,1));
keep = norms>eps0;

X = X(:,keep);
y = y(keep);
norms = norms(keep);

X = X./repmat(norms,size(X,1),1);

class_values = unique(y);
new_y = zeros(length(y),1);

for k = 1:length(class_values)
    new_y(y==class_values(k)) = k;
end

y = new_y;

end

function labeled = make_labeled_split_R2009a( ...
    y,c,fraction,seed)

rand('twister',seed);
labeled = [];

for k = 1:c
    indices = find(y==k);
    count = max(2,ceil(fraction*length(indices)));
    order = randperm(length(indices));
    labeled = [labeled;indices(order(1:count))]; %#ok<AGROW>
end

labeled = sort(labeled);

end

% =========================================================================
% Graph construction and interventions
% =========================================================================

function W = build_binary_graph_blockwise_R2009a( ...
    X,p,block_size)

n = size(X,2);
maximum_entries = n*p;

rows = zeros(maximum_entries,1);
columns = zeros(maximum_entries,1);
position = 0;

for first = 1:block_size:n
    last = min(n,first+block_size-1);
    block = first:last;

    similarities = X(:,block)'*X;

    for local = 1:length(block)
        global_index = block(local);
        similarities(local,global_index) = -Inf;

        [dummy_value,order] = sort( ...
            similarities(local,:),'descend'); %#ok<ASGLU>

        neighbors = order(1:min(p,n-1));
        count = length(neighbors);

        rows(position+1:position+count) = global_index;
        columns(position+1:position+count) = neighbors(:);
        position = position+count;
    end

    fprintf('  graph block %d:%d of %d\n',first,last,n);
    drawnow;
end

rows = rows(1:position);
columns = columns(1:position);

W = sparse(rows,columns,1,n,n);
W = spones(W+W');

end

function [edge_row,edge_col,edge_type,edge_wrong,stats] = ...
    classify_edges_R2009a(W,y,labeled_mask)

[edge_row,edge_col] = find(triu(W,1));
num_edges = length(edge_row);

edge_type = zeros(num_edges,1);
edge_wrong = false(num_edges,1);

for e = 1:num_edges
    i = edge_row(e);
    j = edge_col(e);

    if labeled_mask(i) && labeled_mask(j)
        edge_type(e) = 1;
    elseif xor(labeled_mask(i),labeled_mask(j))
        edge_type(e) = 2;
    else
        edge_type(e) = 3;
    end

    edge_wrong(e) = y(i)~=y(j);
end

wrong_LL = sum(edge_wrong & edge_type==1);
wrong_LU = sum(edge_wrong & edge_type==2);
wrong_UU = sum(edge_wrong & edge_type==3);

correct_LU = sum(~edge_wrong & edge_type==2);
correct_UU = sum(~edge_wrong & edge_type==3);

stats = [ ...
    num_edges,wrong_LL,wrong_LU,wrong_UU, ...
    correct_LU,correct_UU,correct_LU+correct_UU, ...
    mean(~edge_wrong)];

end

function W_new = remove_edge_ids_R2009a( ...
    W,edge_row,edge_col,edge_ids)

W_new = W;

for t = 1:length(edge_ids)
    e = edge_ids(t);
    i = edge_row(e);
    j = edge_col(e);

    W_new(i,j) = 0;
    W_new(j,i) = 0;
end

W_new = spones(W_new);

end

function chosen = choose_matched_edges_R2009a( ...
    candidates,target_count,seed)

if target_count<=0 || isempty(candidates)
    chosen = [];
    return;
end

count = min(target_count,length(candidates));

rand('twister',seed);
order = randperm(length(candidates));
chosen = candidates(order(1:count));

end

% =========================================================================
% GOCNMF
% =========================================================================

function [U0,V0] = initialize_GOCNMF_R2009a( ...
    X,y,labeled,c,seed,eps0)

rand('twister',seed);

[m,n] = size(X);
U0 = rand(m,c)+0.1;
V0 = rand(n,c)+0.1;

for t = 1:length(labeled)
    i = labeled(t);
    V0(i,:) = 0;
    V0(i,y(i)) = 1;
end

U0(U0<eps0) = eps0;
V0(V0<eps0) = eps0;

for t = 1:length(labeled)
    i = labeled(t);
    V0(i,:) = 0;
    V0(i,y(i)) = 1;
end

end

function [U,V,prediction,history,reconstruction_norm] = ...
    train_GOCNMF_fixed_initial_R2009a( ...
    X,y,labeled,c,W,alpha,max_iter,U0,V0,eps0)

U = U0;
V = V0;

n = size(X,2);

labeled_mask = false(n,1);
labeled_mask(labeled) = true;
unlabeled = find(~labeled_mask);

C = zeros(length(labeled),c);

for t = 1:length(labeled)
    C(t,y(labeled(t))) = 1;
end

degree = full(sum(W,2));
D = spdiags(degree,0,n,n);
L = D-W;

history = zeros(max_iter+1,1);
history(1) = goc_objective_R2009a(X,U,V,L,alpha);

for iter = 1:max_iter
    U = U.*((X*V)./(U*(V'*V)+eps0));

    numerator = X'*U+alpha*W*V;
    denominator = V*(U'*U)+alpha*D*V+eps0;

    V(unlabeled,:) = V(unlabeled,:).* ...
        (numerator(unlabeled,:)./denominator(unlabeled,:));

    V(labeled,:) = C;

    history(iter+1) = goc_objective_R2009a( ...
        X,U,V,L,alpha);

    if ~isfinite(history(iter+1))
        error('Nonfinite GOCNMF objective.');
    end
end

[dummy_value,prediction] = max(V,[],2); %#ok<ASGLU>

residual = X-U*V';
reconstruction_norm = norm(residual,'fro');

end

function value = goc_objective_R2009a(X,U,V,L,alpha)

residual = X-U*V';
graph_value = sum(sum(V.*(L*V)));

value = 0.5*sum(residual(:).^2) ...
    +0.5*alpha*graph_value;

end

function count = count_increases_R2009a(history,tolerance)

count = 0;

for i = 1:length(history)-1
    if history(i+1)>history(i) ...
            +tolerance*max(1,abs(history(i)))
        count = count+1;
    end
end

end

% =========================================================================
% Frozen edge-excluded reliability score
% =========================================================================

function [LU_score,LU_wrong,UU_score,UU_wrong] = ...
    compute_edge_excluded_scores_R2009a( ...
    X,U,V,W,y,labeled_mask,alpha,c,eps0)

data_evidence = X'*U;
graph_evidence = alpha*(W*V);

[edge_i,edge_j] = find(triu(W,1));
num_edges = length(edge_i);

LU_score = zeros(num_edges,1);
LU_wrong = false(num_edges,1);
LU_count = 0;

UU_score = zeros(num_edges,1);
UU_wrong = false(num_edges,1);
UU_count = 0;

chance = 1/c;
normalizer = 1-chance;

for e = 1:num_edges
    i = edge_i(e);
    j = edge_j(e);

    if xor(labeled_mask(i),labeled_mask(j))
        LU_count = LU_count+1;

        if labeled_mask(i)
            labeled_node = i;
            unlabeled_node = j;
        else
            labeled_node = j;
            unlabeled_node = i;
        end

        evidence = data_evidence(unlabeled_node,:) ...
            +graph_evidence(unlabeled_node,:) ...
            -alpha*full(W(unlabeled_node,labeled_node)) ...
                *V(labeled_node,:);

        probability = normalize_nonnegative_row_R2009a( ...
            evidence,eps0);

        agreement = probability(y(labeled_node));

        LU_score(LU_count) = ...
            max((agreement-chance)/normalizer,0);
        LU_wrong(LU_count) = ...
            y(unlabeled_node)~=y(labeled_node);

    elseif ~labeled_mask(i) && ~labeled_mask(j)
        UU_count = UU_count+1;

        evidence_i = data_evidence(i,:) ...
            +graph_evidence(i,:) ...
            -alpha*full(W(i,j))*V(j,:);

        evidence_j = data_evidence(j,:) ...
            +graph_evidence(j,:) ...
            -alpha*full(W(j,i))*V(i,:);

        probability_i = normalize_nonnegative_row_R2009a( ...
            evidence_i,eps0);
        probability_j = normalize_nonnegative_row_R2009a( ...
            evidence_j,eps0);

        agreement = probability_i*probability_j';

        UU_score(UU_count) = ...
            max((agreement-chance)/normalizer,0);
        UU_wrong(UU_count) = y(i)~=y(j);
    end
end

LU_score = LU_score(1:LU_count);
LU_wrong = LU_wrong(1:LU_count);

UU_score = UU_score(1:UU_count);
UU_wrong = UU_wrong(1:UU_count);

end

function probability = normalize_nonnegative_row_R2009a( ...
    evidence,eps0)

evidence = max(evidence,0);
total = sum(evidence);

if total<=eps0
    probability = ones(size(evidence))/length(evidence);
else
    probability = evidence/total;
end

end

function values = evaluate_score_R2009a( ...
    score,is_wrong,low_fraction)

score = score(:);
is_wrong = logical(is_wrong(:));

num_edges = length(score);
num_wrong = sum(is_wrong);
num_correct = num_edges-num_wrong;

if num_edges==0 || num_wrong==0 || num_correct==0
    values = [NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
        num_edges,num_wrong,num_correct];
    return;
end

AUC = auc_with_ties_R2009a(score,~is_wrong);

sorted_score = sort(score,'ascend');
quantile_index = max(1,ceil(low_fraction*num_edges));
threshold = sorted_score(quantile_index);

selected = score<=threshold;
selected_count = sum(selected);
selected_wrong = sum(is_wrong & selected);

precision = selected_wrong/max(selected_count,1);
recall = selected_wrong/max(num_wrong,1);
prevalence = num_wrong/num_edges;
enrichment = precision/max(prevalence,eps);
selected_fraction = selected_count/num_edges;
zero_fraction = mean(score<=1e-14);

values = [ ...
    AUC,precision,recall,enrichment,prevalence, ...
    selected_fraction,zero_fraction,num_edges, ...
    num_wrong,num_correct];

end

function AUC = auc_with_ties_R2009a(score,is_positive)

score = score(:);
is_positive = logical(is_positive(:));

n = length(score);
[sorted_score,order] = sort(score,'ascend');

ranks = zeros(n,1);
position = 1;

while position<=n
    last = position;

    while last<n && sorted_score(last+1)==sorted_score(position)
        last = last+1;
    end

    average_rank = 0.5*(position+last);
    ranks(order(position:last)) = average_rank;
    position = last+1;
end

num_positive = sum(is_positive);
num_negative = n-num_positive;

rank_sum_positive = sum(ranks(is_positive));

AUC = (rank_sum_positive ...
    -num_positive*(num_positive+1)/2) ...
    /(num_positive*num_negative);

end

% =========================================================================
% Metrics and output
% =========================================================================

function value = nmi_score_R2009a(true_labels,predicted_labels)

true_labels = true_labels(:);
predicted_labels = predicted_labels(:);

true_classes = unique(true_labels);
predicted_classes = unique(predicted_labels);
n = length(true_labels);

mutual_information = 0;
entropy_true = 0;
entropy_predicted = 0;

for i = 1:length(true_classes)
    count = sum(true_labels==true_classes(i));
    probability = count/n;

    if probability>0
        entropy_true = entropy_true-probability*log(probability);
    end
end

for j = 1:length(predicted_classes)
    count = sum(predicted_labels==predicted_classes(j));
    probability = count/n;

    if probability>0
        entropy_predicted = ...
            entropy_predicted-probability*log(probability);
    end
end

for i = 1:length(true_classes)
    true_mask = true_labels==true_classes(i);
    true_count = sum(true_mask);

    for j = 1:length(predicted_classes)
        predicted_mask = ...
            predicted_labels==predicted_classes(j);
        predicted_count = sum(predicted_mask);
        joint_count = sum(true_mask & predicted_mask);

        if joint_count>0
            mutual_information = mutual_information ...
                +(joint_count/n)*log( ...
                (joint_count*n)/(true_count*predicted_count));
        end
    end
end

if entropy_true+entropy_predicted<=eps
    value = 1;
else
    value = 2*mutual_information/ ...
        (entropy_true+entropy_predicted);
end

end

function write_dataset_raw_R2009a( ...
    file_name,dataset_name,seeds,condition_names, ...
    condition_scores,edge_stats,reliability_metrics, ...
    rescue_harm,completed)

fid = fopen(file_name,'w');

if fid<0
    error('Cannot create dataset raw CSV.');
end

fprintf(fid,['dataset,seed,condition,ACC,NMI,RE,' ...
    'objective_increases,runtime,total_edges,wrong_LL,' ...
    'wrong_LU,wrong_UU,correct_LU,correct_UU,' ...
    'correct_active,edge_purity,rescue_rate,harm_rate,' ...
    'LU_AUC,LU_precision,LU_recall,LU_enrichment,' ...
    'UU_AUC,UU_precision,UU_recall,UU_enrichment\n']);

for seed_id = 1:length(seeds)
    if ~completed(seed_id)
        continue;
    end

    LU = squeeze(reliability_metrics(seed_id,1,:));
    UU = squeeze(reliability_metrics(seed_id,2,:));

    for condition_id = 1:length(condition_names)
        values = squeeze( ...
            condition_scores(seed_id,condition_id,:));

        fprintf(fid,['%s,%d,%s,%.16g,%.16g,%.16g,' ...
            '%d,%.16g,%d,%d,%d,%d,%d,%d,%d,%.16g,' ...
            '%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,' ...
            '%.16g,%.16g,%.16g,%.16g\n'], ...
            dataset_name,seeds(seed_id), ...
            condition_names{condition_id}, ...
            values(1),values(2),values(3),round(values(4)), ...
            values(5),round(edge_stats(seed_id,1)), ...
            round(edge_stats(seed_id,2)), ...
            round(edge_stats(seed_id,3)), ...
            round(edge_stats(seed_id,4)), ...
            round(edge_stats(seed_id,5)), ...
            round(edge_stats(seed_id,6)), ...
            round(edge_stats(seed_id,7)), ...
            edge_stats(seed_id,8), ...
            rescue_harm(seed_id,1),rescue_harm(seed_id,2), ...
            LU(1),LU(2),LU(3),LU(4), ...
            UU(1),UU(2),UU(3),UU(4));
    end
end

fclose(fid);

end
