function prepare_NEU_CLS_LBP_R2009a(imageRoot,outputMat)
%PREPARE_NEU_CLS_LBP_R2009A
% Extract a nonnegative spatial uniform-LBP representation from NEU-CLS.
%
% Representation:
%   - each image is converted to grayscale and resized to 34-by-34;
%   - LBP_{8,1}^{u2} is computed on the valid 32-by-32 center region;
%   - the region is divided into 4-by-4 nonoverlapping cells;
%   - each cell uses a 59-bin uniform-LBP histogram;
%   - feature dimension = 4*4*59 = 944;
%   - each cell histogram is L1 normalized;
%   - the final vector is nonnegative and is L2 normalized by the
%     experiment program.
%
% The program does not require the Image Processing Toolbox imresize
% function and is compatible with MATLAB R2009a.
%
% Usage:
%   prepare_NEU_CLS_LBP_R2009a( ...
%       'D:\datasets\NEU-CLS', ...
%       fullfile(pwd,'NEU_CLS_LBP59_4x4.mat'));

if nargin<1 || isempty(imageRoot)
    imageRoot=uigetdir(pwd,'Select the NEU-CLS image folder');
    if isequal(imageRoot,0)
        error('Image-folder selection was cancelled.');
    end
end
if nargin<2 || isempty(outputMat)
    outputMat=fullfile(pwd,'NEU_CLS_LBP59_4x4.mat');
end
if exist(imageRoot,'dir')~=7
    error('Image directory not found: %s',imageRoot);
end

class_codes={'Cr','In','Pa','PS','RS','Sc'};
class_names={'Crazing','Inclusion','Patches','Pitted surface', ...
    'Rolled-in scale','Scratches'};

files=neu_lbp_collect_images(imageRoot);
if isempty(files)
    error('No NEU-CLS images were found under %s.',imageRoot);
end

n=length(files);
feature_dimension=4*4*59;
fea=zeros(n,feature_dimension);
gnd=zeros(n,1);
image_files=cell(n,1);
counts=zeros(6,1);
uniform_mapping=neu_lbp_uniform_mapping();

fprintf('\nNEU-CLS LBP FEATURE PREPARATION\n');
fprintf('Images found: %d\n',n);
fprintf('Feature: LBP_{8,1}^{u2}, 4-by-4 cells, 59 bins/cell.\n');
fprintf('Feature dimension: %d\n\n',feature_dimension);

for i=1:n
    label=neu_lbp_infer_label(files{i});
    A=double(imread(files{i}));

    if ndims(A)==3
        A=mean(A,3);
    elseif ndims(A)~=2
        error('Unsupported image shape: %s',files{i});
    end
    if isempty(A) || any(~isfinite(A(:)))
        error('Invalid image: %s',files{i});
    end

    A=A-min(A(:));
    maximumValue=max(A(:));
    if maximumValue>0
        A=A/maximumValue;
    end

    B=neu_lbp_resize_bilinear(A,34,34);
    B=min(1,max(0,B));
    feature=neu_lbp_spatial_histogram(B,uniform_mapping);

    if length(feature)~=feature_dimension
        error('Unexpected feature dimension for %s.',files{i});
    end
    if any(~isfinite(feature)) || any(feature<0)
        error('Invalid LBP feature for %s.',files{i});
    end

    fea(i,:)=feature;
    gnd(i)=label;
    image_files{i}=files{i};
    counts(label)=counts(label)+1;

    if mod(i,100)==0 || i==n
        fprintf('  processed %d/%d\n',i,n);
    end
end

if n~=1800
    error('Expected 1800 NEU-CLS images, found %d.',n);
end
if any(counts~=300)
    disp([class_codes(:) num2cell(counts)]);
    error('Each class must contain exactly 300 images.');
end
if any(~isfinite(fea(:))) || min(fea(:))<-1e-12
    error('Prepared features are not finite and nonnegative.');
end
if any(sum(fea,2)<=0)
    error('At least one image has a zero LBP feature vector.');
end

feature_name='spatial uniform LBP_{8,1}^{u2}';
image_size=[34 34];
valid_lbp_size=[32 32];
cell_grid=[4 4];
bins_per_cell=59;

save(outputMat,'fea','gnd','class_names','class_codes', ...
    'image_files','counts','feature_name','feature_dimension', ...
    'image_size','valid_lbp_size','cell_grid','bins_per_cell','-v7');

fprintf('\nOutput: %s\n',outputMat);
fprintf('NEU_LBP_DIMENSION=%d\n',feature_dimension);
fprintf('NEU_LBP_NONNEGATIVE_PASS=1\n');
fprintf('NEU_LBP_CLASS_BALANCE_PASS=1\n');
fprintf('NEU_LBP_PREPARATION_COMPLETE=1\n');
end


function feature=neu_lbp_spatial_histogram(B,mapping)
if any(size(B)~=[34 34])
    error('The resized image must be 34 by 34.');
end

center=B(2:33,2:33);
code=zeros(32,32);

neighbors=cell(8,1);
neighbors{1}=B(1:32,1:32);
neighbors{2}=B(1:32,2:33);
neighbors{3}=B(1:32,3:34);
neighbors{4}=B(2:33,3:34);
neighbors{5}=B(3:34,3:34);
neighbors{6}=B(3:34,2:33);
neighbors{7}=B(3:34,1:32);
neighbors{8}=B(2:33,1:32);

for b=1:8
    code=code+(neighbors{b}>=center)*2^(b-1);
end

mapped=reshape(mapping(code(:)+1),32,32);
feature=zeros(1,4*4*59);
position=0;

for rowCell=1:4
    rowIndex=(rowCell-1)*8+(1:8);
    for colCell=1:4
        colIndex=(colCell-1)*8+(1:8);
        values=mapped(rowIndex,colIndex);
        histogram=accumarray(values(:),1,[59 1]);
        total=sum(histogram);
        if total>0
            histogram=histogram/total;
        end
        feature(position+(1:59))=histogram';
        position=position+59;
    end
end
end


function mapping=neu_lbp_uniform_mapping()
mapping=59*ones(256,1);
nextBin=1;

for code=0:255
    bits=zeros(1,8);
    for b=1:8
        bits(b)=bitget(code,b);
    end
    shifted=bits([2:8 1]);
    transitions=sum(bits~=shifted);
    if transitions<=2
        mapping(code+1)=nextBin;
        nextBin=nextBin+1;
    end
end

if nextBin~=59
    error('Uniform-LBP mapping must contain exactly 58 uniform bins.');
end
end


function files=neu_lbp_collect_images(folder)
files=cell(0,1);
entries=dir(folder);

for i=1:length(entries)
    name=entries(i).name;
    if strcmp(name,'.') || strcmp(name,'..')
        continue;
    end

    path=fullfile(folder,name);
    if entries(i).isdir
        files=[files;neu_lbp_collect_images(path)]; %#ok<AGROW>
    else
        [dummy1,dummy2,extension]=fileparts(name); %#ok<ASGLU>
        extension=lower(extension);
        if strcmp(extension,'.bmp') || strcmp(extension,'.png') || ...
                strcmp(extension,'.jpg') || strcmp(extension,'.jpeg') || ...
                strcmp(extension,'.tif') || strcmp(extension,'.tiff')
            files{end+1,1}=path; %#ok<AGROW>
        end
    end
end

if ~isempty(files)
    files=sort(files);
end
end


function label=neu_lbp_infer_label(path)
[parent,base,dummy]=fileparts(path); %#ok<ASGLU>
[dummy,parentName]=fileparts(parent); %#ok<ASGLU>
base=upper(base);
parentName=lower(parentName);
label=0;

if neu_lbp_has_prefix(base,'CR')
    label=1;
elseif neu_lbp_has_prefix(base,'IN')
    label=2;
elseif neu_lbp_has_prefix(base,'PA')
    label=3;
elseif neu_lbp_has_prefix(base,'PS')
    label=4;
elseif neu_lbp_has_prefix(base,'RS')
    label=5;
elseif neu_lbp_has_prefix(base,'SC')
    label=6;
end

if label==0
    if ~isempty(strfind(parentName,'craz'))
        label=1;
    elseif ~isempty(strfind(parentName,'inclu'))
        label=2;
    elseif ~isempty(strfind(parentName,'patch'))
        label=3;
    elseif ~isempty(strfind(parentName,'pitted')) || strcmp(parentName,'ps')
        label=4;
    elseif ~isempty(strfind(parentName,'rolled')) || ...
            ~isempty(strfind(parentName,'scale')) || strcmp(parentName,'rs')
        label=5;
    elseif ~isempty(strfind(parentName,'scratch')) || strcmp(parentName,'sc')
        label=6;
    end
end

if label==0
    error('Cannot infer the NEU-CLS class from %s.',path);
end
end


function tf=neu_lbp_has_prefix(name,code)
tf=0;
if length(name)<length(code) || ~strcmp(name(1:length(code)),code)
    return;
end
if length(name)==length(code)
    tf=1;
    return;
end
character=name(length(code)+1);
tf=character=='_' || character=='-' || ...
    (character>='0' && character<='9');
end


function B=neu_lbp_resize_bilinear(A,newRows,newCols)
[oldRows,oldCols]=size(A);
if oldRows==newRows && oldCols==newCols
    B=A;
    return;
end

xNew=linspace(1,oldCols,newCols);
yNew=linspace(1,oldRows,newRows)';
[Xnew,Ynew]=meshgrid(xNew,yNew);
B=interp2(1:oldCols,(1:oldRows)',A,Xnew,Ynew,'linear');

if any(~isfinite(B(:)))
    error('Bilinear resizing produced nonfinite values.');
end
end
