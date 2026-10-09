function prepare_NEU_CLS_official_R2019a(workRoot,outputRoot)
%PREPARE_NEU_CLS_OFFICIAL_R2019A Download/extract NEU-CLS and build MAT files.
% The official provider uses Google Drive and a RAR archive. If automated
% download or 7-Zip extraction is blocked, this function opens the official
% provider page and reports the exact manual step instead of using a mirror.

if nargin<1 || isempty(workRoot),workRoot=fullfile(pwd,'official_raw','neu_cls');end
if nargin<2 || isempty(outputRoot),outputRoot=fullfile(pwd,'data','neu_cls');end
if exist(workRoot,'dir')~=7,mkdir(workRoot);end
if exist(outputRoot,'dir')~=7,mkdir(outputRoot);end

imageRoot=find_image_root(workRoot);
if isempty(imageRoot)
    archive=fullfile(workRoot,'NEU-CLS.rar');
    if exist(archive,'file')~=2
        id='1NGlXT9sIaQpyxUoT6MLKm1Pr6x8oxOvc';
        url=['https://drive.usercontent.google.com/download?id=' id '&export=download&confirm=t'];
        try, websave(archive,url); catch, open_official_and_fail(workRoot); end
        info=dir(archive); if isempty(info)||info.bytes<1000000,delete(archive);open_official_and_fail(workRoot);end
    end
    seven=find_7zip();
    if isempty(seven)
        error(['NEU-CLS.rar was downloaded to %s. Install 7-Zip or extract it ' ...
            'under %s, then rerun this function.'],archive,workRoot);
    end
    command=sprintf('"%s" x -y -o"%s" "%s"',seven,workRoot,archive);
    status=system(command); if status~=0,error('7-Zip failed to extract %s.',archive);end
    imageRoot=find_image_root(workRoot);
end
if isempty(imageRoot),error('Could not locate the 1800 extracted NEU-CLS images.');end

base=fileparts(mfilename('fullpath')); neuCode=fullfile(base,'..','neu_cls'); addpath(neuCode);
prepare_NEU_CLS_R2009a(imageRoot,fullfile(outputRoot,'NEU_CLS_32x32.mat'));
prepare_NEU_CLS_LBP_R2009a(imageRoot,fullfile(outputRoot,'NEU_CLS_LBP59_4x4.mat'));
fprintf('NEU_OFFICIAL_PIPELINE_PASS=1\n');
end

function root=find_image_root(base)
root=''; files=collect(base); if numel(files)==1800,root=base;return;end
d=dir(base); for i=1:numel(d),if d(i).isdir&&~any(strcmp(d(i).name,{'.','..'})),p=fullfile(base,d(i).name);if numel(collect(p))==1800,root=p;return;end,end,end
end
function files=collect(root)
files=cell(0,1);d=dir(root);for i=1:numel(d),if any(strcmp(d(i).name,{'.','..'})),continue;end;p=fullfile(root,d(i).name);if d(i).isdir,files=[files;collect(p)];else,[~,~,e]=fileparts(p);if any(strcmpi(e,{'.bmp','.png','.jpg','.jpeg','.tif','.tiff'})),files{end+1,1}=p;end,end,end %#ok<AGROW>
end
function p=find_7zip()
p=''; candidates={'C:\Program Files\7-Zip\7z.exe','C:\Program Files (x86)\7-Zip\7z.exe'};
for i=1:numel(candidates),if exist(candidates{i},'file')==2,p=candidates{i};return;end,end
end
function open_official_and_fail(workRoot)
web('https://faculty.neu.edu.cn/songkc/en/zdylm/263265/list/','-browser');
error(['Google Drive blocked the automatic download. Download NEU-CLS.rar from ' ...
    'the official NEU page, place it in %s, and rerun.'],workRoot);
end
