try
    root = fileparts(mfilename('fullpath'));
    dataset = getenv('SNMFWLP_DATASET');
    outputDir = getenv('SNMFWLP_OUTPUT');
    alpha = str2double(getenv('SNMFWLP_ALPHA'));
    beta = str2double(getenv('SNMFWLP_BETA'));
    seedStart = str2double(getenv('SNMFWLP_SEED_START'));
    if isempty(dataset) || isempty(outputDir) || ~isfinite(alpha) || ...
            ~isfinite(beta) || ~isfinite(seedStart)
        error('Missing or invalid SNMFWLP pilot environment variable.');
    end
    run_SNMFWLP_common_stop_R2019a(1,fullfile(root,'data_ascii'), ...
        outputDir,dataset,alpha,beta,true,seedStart);
    exit(0);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
