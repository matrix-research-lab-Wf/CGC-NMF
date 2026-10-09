function prepare_PIE_YaleB_R2019a(datasetName,imageRoot,outputMat)
%PREPARE_PIE_YALEB_R2019A Convert lawfully obtained official face images.
% PIE: use the frontal c27 camera; exactly 42 images per each of 68 people.
% YaleB: use Cropped Extended Yale B frontal images, exclude Ambient images.
% Raw downloads require the provider's access/permission procedure and are
% deliberately not bypassed by this program.

if nargin<2 || isempty(imageRoot), error('imageRoot is required.'); end
if exist(imageRoot,'dir')~=7, error('Image directory not found: %s.',imageRoot); end
name=upper(datasetName);
if nargin<3 || isempty(outputMat), outputMat=fullfile(pwd,[name '.mat']); end
files=collect_images(imageRoot);

if strcmp(name,'PIE')
    keep=false(numel(files),1); labels=zeros(numel(files),1);
    for i=1:numel(files)
        low=lower(strrep(files{i},'\','/'));
        tok=regexp(low,'(?:^|/)([0-9]{5})(?:/|_)','tokens','once');
        keep(i)=~isempty(strfind(low,'c27')) && ~isempty(tok); %#ok<STREMP>
        if keep(i), labels(i)=str2double(tok{1}); end
    end
    files=files(keep); labels=labels(keep);
    people=unique(labels); if numel(people)~=68,error('Expected 68 PIE subjects; found %d.',numel(people));end
    selected=cell(0,1); gnd=zeros(0,1);
    for k=1:numel(people)
        f=sort(files(labels==people(k)));
        if numel(f)~=42,error('PIE subject %d has %d c27 images; expected 42.',people(k),numel(f));end
        selected=[selected;f(:)]; gnd=[gnd;repmat(k,42,1)]; %#ok<AGROW>
    end
    files=selected; source_url='https://www.cs.cmu.edu/afs/cs/project/vision/vasc/idb/www/html/face/';
elseif strcmp(name,'YALEB')
    keep=false(numel(files),1); labels=zeros(numel(files),1);
    for i=1:numel(files)
        low=lower(strrep(files{i},'\','/'));
        tok=regexp(low,'yaleb([0-9]{2})','tokens','once');
        keep(i)=isempty(strfind(low,'ambient')) && ~isempty(tok); %#ok<STREMP>
        if keep(i), labels(i)=str2double(tok{1}); end
    end
    files=files(keep); labels=labels(keep); [labels,ord]=sort(labels); files=files(ord);
    people=unique(labels); if numel(people)~=38,error('Expected 38 Extended Yale B subjects; found %d.',numel(people));end
    gnd=zeros(size(labels)); for k=1:numel(people),gnd(labels==people(k))=k;end
    if numel(files)~=2414,error('Expected 2414 non-Ambient cropped YaleB images; found %d.',numel(files));end
    source_url='https://vision.ucsd.edu/datasets/extended-yale-face-database-b-b';
else, error('datasetName must be PIE or YaleB.'); end

fea=zeros(numel(files),1024);
for i=1:numel(files)
    A=double(imread(files{i})); if ndims(A)==3,A=0.2989360213*A(:,:,1)+0.5870430745*A(:,:,2)+0.1140209043*A(:,:,3);end
    B=resize_nearest(A,32,32); fea(i,:)=round(min(255,max(0,B(:)')));
end
save(outputMat,'fea','gnd','source_url','-v7');
fprintf('Output: %s\n%s_PREPARATION_PASS=1\n',outputMat,name);
warning('Dimensions are verified; byte identity with legacy MAT preprocessing must be checked separately.');
end

function files=collect_images(root)
files=cell(0,1); d=dir(root);
for i=1:numel(d)
    if strcmp(d(i).name,'.')||strcmp(d(i).name,'..'),continue;end
    p=fullfile(root,d(i).name);
    if d(i).isdir,files=[files;collect_images(p)]; %#ok<AGROW>
    else,[~,~,e]=fileparts(d(i).name);if any(strcmpi(e,{'.pgm','.png','.jpg','.jpeg','.bmp','.tif','.tiff'})),files{end+1,1}=p;end,end %#ok<AGROW>
end
end
function B=resize_nearest(A,nr,nc)
r=min(size(A,1),max(1,round(((1:nr)-0.5)*size(A,1)/nr+0.5)));
c=min(size(A,2),max(1,round(((1:nc)-0.5)*size(A,2)/nc+0.5))); B=A(r,c);
end
