function prepare_all_official_data_R2019a(rawRoot,outputRoot)
%PREPARE_ALL_OFFICIAL_DATA_R2019A Rebuild all legally obtainable datasets.
% COIL20, COIL100 and MNIST download automatically. NEU-CLS attempts the
% official Google Drive download. PIE and Extended Yale B require provider
% access/permission and are converted when extracted folders are present at
% rawRoot/PIE and rawRoot/YaleB.

if nargin<1||isempty(rawRoot),rawRoot=fullfile(pwd,'official_raw');end
if nargin<2||isempty(outputRoot),outputRoot=fullfile(pwd,'data');end
mainDir=fullfile(outputRoot,'main_six');neuDir=fullfile(outputRoot,'neu_cls');
if exist(mainDir,'dir')~=7,mkdir(mainDir);end;if exist(neuDir,'dir')~=7,mkdir(neuDir);end

prepare_COIL_from_official_R2019a('COIL20',fullfile(rawRoot,'coil20'),fullfile(mainDir,'COIL20_Obj.mat'));
prepare_COIL_from_official_R2019a('COIL100',fullfile(rawRoot,'coil100'),fullfile(mainDir,'COIL100_Obj.mat'));
prepare_MNIST_Han_R2019a(fullfile(rawRoot,'mnist'),fullfile(mainDir,'MNIST_Han.mat'));

pieRoot=fullfile(rawRoot,'PIE');
if exist(pieRoot,'dir')==7,prepare_PIE_YaleB_R2019a('PIE',pieRoot,fullfile(mainDir,'PIE.mat'));
else,fprintf('PIE_PENDING: obtain CMU PIE through the provider and extract to %s\n',pieRoot);end
yaleRoot=fullfile(rawRoot,'YaleB');
if exist(yaleRoot,'dir')==7,prepare_PIE_YaleB_R2019a('YaleB',yaleRoot,fullfile(mainDir,'YaleB.mat'));
else,fprintf('YALEB_PENDING: obtain Cropped Extended Yale B and extract to %s\n',yaleRoot);end

try,prepare_NEU_CLS_official_R2019a(fullfile(rawRoot,'neu_cls'),neuDir);
catch ME,warning('NEU-CLS preparation pending: %s',ME.message);end
fprintf('OFFICIAL_DATA_PIPELINE_FINISHED=1\n');
end
