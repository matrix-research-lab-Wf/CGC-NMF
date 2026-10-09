try
    codeDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(fileparts(codeDir));
    dataDir = fullfile(rootDir,'data','main_six');
    outputRoot = fullfile(rootDir,'results','recovered_latest_smoke_v2');
    if exist(outputRoot,'dir') ~= 7
        mkdir(outputRoot);
    end

    run_common_stop_core_R2019a(1,dataDir, ...
        fullfile(outputRoot,'common_stop_COIL20'),'COIL20');
    run_SNMFWLP_common_stop_R2019a(1,dataDir, ...
        fullfile(outputRoot,'snmfwlp_COIL20'),'COIL20',1,1,false,20260617);

    fid = fopen(fullfile(outputRoot,'RECOVERED_EXECUTION_COMPLETE.txt'),'w');
    fprintf(fid,'RECOVERED_SCRIPTS_EXECUTED=1\n');
    fprintf(fid,'NOTE=This marker confirms execution only; inspect each decision file for scientific acceptance gates.\n');
    fclose(fid);
    exit(0);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
