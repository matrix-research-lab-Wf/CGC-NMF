try
    root = fileparts(mfilename('fullpath'));
    dataset = getenv('SNMFWLP_DATASET');
    outputDir = getenv('SNMFWLP_OUTPUT');
    alpha = str2double(getenv('SNMFWLP_ALPHA'));
    beta = str2double(getenv('SNMFWLP_BETA'));
    if isempty(dataset) || isempty(outputDir) || ~isfinite(alpha) || ~isfinite(beta)
        error('Missing or invalid SNMFWLP formal-run environment variable.');
    end
    run_SNMFWLP_common_stop_R2019a(3,fullfile(root,'data_ascii'), ...
        outputDir,dataset,alpha,beta,false,20260617);
    exit(0);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
