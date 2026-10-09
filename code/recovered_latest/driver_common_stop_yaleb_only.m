try
    codeDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(fileparts(codeDir));
    run_common_stop_core_R2019a(3,fullfile(rootDir,'data','main_six'), ...
        fullfile(rootDir,'results','recovered_six_dataset_3seed','common_stop_YaleB'), ...
        'YaleB');
    exit(0);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
