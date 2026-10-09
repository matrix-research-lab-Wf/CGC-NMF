try
    codeDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(fileparts(codeDir));
    run_SNMFWLP_common_stop_R2019a(3,fullfile(rootDir,'data','main_six'), ...
        fullfile(rootDir,'results','recovered_six_dataset_3seed','snmfwlp_COIL100'), ...
        'COIL100',10,100,false,20260617);
    exit(0);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
