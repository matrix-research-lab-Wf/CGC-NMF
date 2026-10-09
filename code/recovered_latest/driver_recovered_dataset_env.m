try
    codeDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(fileparts(codeDir));
    dataDir = fullfile(rootDir,'data','main_six');
    outputRoot = fullfile(rootDir,'results','recovered_six_dataset_3seed');
    dataset = getenv('RECOVERED_DATASET');
    alpha = str2double(getenv('RECOVERED_SNM_ALPHA'));
    beta = str2double(getenv('RECOVERED_SNM_BETA'));
    if isempty(dataset) || ~isfinite(alpha) || ~isfinite(beta)
        error('Missing or invalid recovered-run environment variables.');
    end
    run_common_stop_core_R2019a(3,dataDir, ...
        fullfile(outputRoot,['common_stop_' dataset]),dataset);
    run_SNMFWLP_common_stop_R2019a(3,dataDir, ...
        fullfile(outputRoot,['snmfwlp_' dataset]),dataset, ...
        alpha,beta,false,20260617);
    exit(0);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
