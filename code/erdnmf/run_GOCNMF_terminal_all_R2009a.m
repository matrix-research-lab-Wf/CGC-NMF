function run_GOCNMF_terminal_all_R2009a
% Run all six datasets sequentially, then aggregate.
%
% The worker is resumable. Re-running this master continues unfinished
% splits rather than restarting completed ones.

close all;
clc;

config = GOCNMF_terminal_validation_config_R2009a;

for dataset_id = 1:length(config.dataset_names)
    fprintf('\n====================================================\n');
    fprintf('Starting dataset %d/%d: %s\n', ...
        dataset_id,length(config.dataset_names), ...
        config.dataset_names{dataset_id});
    fprintf('====================================================\n');

    run_GOCNMF_terminal_dataset_R2009a(dataset_id);
end

summarize_GOCNMF_terminal_validation_R2009a;

end
