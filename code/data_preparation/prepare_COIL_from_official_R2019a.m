function prepare_COIL_from_official_R2019a(datasetName, workRoot, outputMat)
%PREPARE_COIL_FROM_OFFICIAL_R2019A Download and convert COIL-20/COIL-100.
% Requires MATLAB R2019a. Image Processing Toolbox is optional: this file
% implements nearest-neighbour resizing and RGB-to-gray conversion locally.
%
% Examples:
%   prepare_COIL_from_official_R2019a('COIL20',pwd,'COIL20_Obj.mat')
%   prepare_COIL_from_official_R2019a('COIL100',pwd,'COIL100_Obj.mat')

if nargin<1 || isempty(datasetName), error('datasetName is required.'); end
if nargin<2 || isempty(workRoot), workRoot=fullfile(pwd,'official_raw'); end
if nargin<3 || isempty(outputMat), outputMat=[upper(datasetName) '_Obj.mat']; end
if exist(workRoot,'dir')~=7, mkdir(workRoot); end

name=upper(strrep(datasetName,'-',''));
if strcmp(name,'COIL20')
    url='https://www.cs.columbia.edu/CAVE/databases/SLAM_coil-20_coil-100/coil-20/coil-20-proc.zip';
    archive=fullfile(workRoot,'coil-20-proc.zip'); expectedClasses=20;
elseif strcmp(name,'COIL100')
    url='https://www.cs.columbia.edu/CAVE/databases/SLAM_coil-20_coil-100/coil-100/coil-100.zip';
    archive=fullfile(workRoot,'coil-100.zip'); expectedClasses=100;
else
    error('datasetName must be COIL20 or COIL100.');
end

if exist(archive,'file')~=2
    fprintf('Downloading %s from Columbia CAVE...\n',name);
    websave(archive,url);
end
extractRoot=fullfile(workRoot,lower(name));
if exist(extractRoot,'dir')~=7, mkdir(extractRoot); unzip(archive,extractRoot); end
files=collect_images(extractRoot);
if numel(files)~=expectedClasses*72
    error('Expected %d images, found %d under %s.',expectedClasses*72,numel(files),extractRoot);
end

labels=zeros(numel(files),1); angles=zeros(numel(files),1);
for i=1:numel(files)
    [~,base]=fileparts(files{i}); token=regexp(base,'^obj([0-9]+)__([0-9]+)$','tokens','once');
    if isempty(token), error('Unexpected COIL filename: %s.',files{i}); end
    labels(i)=str2double(token{1}); angles(i)=str2double(token{2});
end
[~,order]=sortrows([labels angles],[1 2]); files=files(order); labels=labels(order);

fea=zeros(numel(files),1024); gnd=labels;
for i=1:numel(files)
    A=double(imread(files{i}));
    if ndims(A)==3
        % ITU-R BT.601 luma, implemented without Image Processing Toolbox.
        A=0.2989360213*A(:,:,1)+0.5870430745*A(:,:,2)+0.1140209043*A(:,:,3);
    end
    B=resize_nearest(A,32,32);
    fea(i,:)=round(min(255,max(0,B(:)')));
end
if any(~isfinite(fea(:))) || any(fea(:)<0), error('Invalid COIL features.'); end
source_url=url; preprocessing='grayscale, nearest-neighbour resize to 32x32, MATLAB column-major vectorization'; %#ok<NASGU>
save(outputMat,'fea','gnd','source_url','preprocessing','-v7');
fprintf('Output: %s\n',outputMat);
fprintf('%s_PREPARATION_PASS=1\n',name);
if strcmp(name,'COIL100')
    warning(['The archived COIL100_Obj.mat was produced by an older undocumented ' ...
        'grayscale pipeline. Dimensions and labels match, but byte identity is not claimed.']);
end
end

function files=collect_images(root)
files=cell(0,1); d=dir(root);
for i=1:numel(d)
    if strcmp(d(i).name,'.') || strcmp(d(i).name,'..'), continue; end
    p=fullfile(root,d(i).name);
    if d(i).isdir, files=[files;collect_images(p)]; %#ok<AGROW>
    else
        [~,~,e]=fileparts(d(i).name);
        if any(strcmpi(e,{'.png','.ppm','.pgm','.jpg','.jpeg','.bmp'})), files{end+1,1}=p; end %#ok<AGROW>
    end
end
end

function B=resize_nearest(A,nr,nc)
r=min(size(A,1),max(1,round(((1:nr)-0.5)*size(A,1)/nr+0.5)));
c=min(size(A,2),max(1,round(((1:nc)-0.5)*size(A,2)/nc+0.5)));
B=A(r,c);
end
