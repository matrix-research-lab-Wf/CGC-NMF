try
    codeDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(fileparts(codeDir));
    dataDir = fullfile(rootDir,'data','main_six');
    outputRoot = fullfile(rootDir,'results','recovered_six_dataset_3seed');
    if exist(outputRoot,'dir') ~= 7
        mkdir(outputRoot);
    end

    datasets = {'PIE','YaleB','COIL20','Optdigits','MNIST','COIL100'};
    snmAlpha = [100,10,10,10,10,10];
    snmBeta = [100,100,100,100,100,100];

    for i = 1:numel(datasets)
        name = datasets{i};
        run_common_stop_core_R2019a(3,dataDir, ...
            fullfile(outputRoot,['common_stop_' name]),name);
        run_SNMFWLP_common_stop_R2019a(3,dataDir, ...
            fullfile(outputRoot,['snmfwlp_' name]),name, ...
            snmAlpha(i),snmBeta(i),false,20260617);
    end

    fid = fopen(fullfile(outputRoot,'EXECUTION_COMPLETE.txt'),'w');
    fprintf(fid,'SIX_DATASET_THREE_SEED_EXECUTION_COMPLETE=1\n');
    fprintf(fid,'NOTE=Inspect all decision files; execution completion is not a convergence claim.\n');
    fclose(fid);
    exit(0);
catch ME
    disp(getReport(ME,'extended'));
    exit(1);
end
