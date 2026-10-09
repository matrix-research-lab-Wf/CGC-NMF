function prepare_NEU_CLS_R2009a(imageRoot, outputMat)
%PREPARE_NEU_CLS_R2009A Convert NEU-CLS images to a 32x32 MAT file.
% Supports a flat folder with Cr_1.bmp, In_1.bmp, Pa_1.bmp, PS_1.bmp,
% RS_1.bmp, Sc_1.bmp, or class subfolders.

if nargin < 1 || isempty(imageRoot), imageRoot = pwd; end
if nargin < 2 || isempty(outputMat)
    outputMat = fullfile(pwd,'NEU_CLS_32x32.mat');
end
if exist(imageRoot,'dir') ~= 7
    error('Image directory not found: %s',imageRoot);
end

class_codes = {'Cr','In','Pa','PS','RS','Sc'};
class_names = {'Crazing','Inclusion','Patches','Pitted surface', ...
    'Rolled-in scale','Scratches'};
files = neu_collect_images(imageRoot);
if isempty(files), error('No NEU-CLS images found under %s.',imageRoot); end

n = length(files);
fea = zeros(n,1024);
gnd = zeros(n,1);
image_files = cell(n,1);
counts = zeros(6,1);

fprintf('Images found: %d\n',n);
for i = 1:n
    label = neu_infer_label(files{i});
    A = double(imread(files{i}));
    if ndims(A)==3
        A = mean(A,3);
    elseif ndims(A)~=2
        error('Unsupported image shape: %s',files{i});
    end
    if isempty(A) || any(~isfinite(A(:)))
        error('Invalid image: %s',files{i});
    end
    A = A-min(A(:));
    maximumValue = max(A(:));
    if maximumValue>0, A = A/maximumValue; end
    B = neu_resize_bilinear(A,32,32);
    B = min(1,max(0,B));
    fea(i,:) = B(:)';
    gnd(i) = label;
    image_files{i} = files{i};
    counts(label) = counts(label)+1;
    if mod(i,100)==0 || i==n
        fprintf('  processed %d/%d\n',i,n);
    end
end

if n~=1800, error('Expected 1800 images, found %d.',n); end
if any(counts~=300)
    disp([class_codes(:) num2cell(counts)]);
    error('Each class must contain exactly 300 images.');
end
if any(~isfinite(fea(:))) || min(fea(:)) < -1e-12
    error('Prepared features are invalid.');
end

save(outputMat,'fea','gnd','class_names','class_codes', ...
    'image_files','counts','-v7');
fprintf('Output: %s\n',outputMat);
fprintf('NEU_DATA_PREPARATION_PASS=1\n');
end

function files = neu_collect_images(folder)
files = cell(0,1);
entries = dir(folder);
for i = 1:length(entries)
    name = entries(i).name;
    if strcmp(name,'.') || strcmp(name,'..'), continue; end
    path = fullfile(folder,name);
    if entries(i).isdir
        files = [files; neu_collect_images(path)]; %#ok<AGROW>
    else
        [dummy1,dummy2,ext] = fileparts(name); %#ok<ASGLU>
        ext = lower(ext);
        if strcmp(ext,'.bmp') || strcmp(ext,'.png') || ...
                strcmp(ext,'.jpg') || strcmp(ext,'.jpeg') || ...
                strcmp(ext,'.tif') || strcmp(ext,'.tiff')
            files{end+1,1} = path; %#ok<AGROW>
        end
    end
end
if ~isempty(files), files = sort(files); end
end

function label = neu_infer_label(path)
[parent,base,dummy] = fileparts(path); %#ok<ASGLU>
[dummy,parentName] = fileparts(parent); %#ok<ASGLU>
base = upper(base);
parentName = lower(parentName);
label = 0;
if neu_has_prefix(base,'CR'), label = 1;
elseif neu_has_prefix(base,'IN'), label = 2;
elseif neu_has_prefix(base,'PA'), label = 3;
elseif neu_has_prefix(base,'PS'), label = 4;
elseif neu_has_prefix(base,'RS'), label = 5;
elseif neu_has_prefix(base,'SC'), label = 6;
end
if label==0
    if ~isempty(strfind(parentName,'craz')), label = 1;
    elseif ~isempty(strfind(parentName,'inclu')), label = 2;
    elseif ~isempty(strfind(parentName,'patch')), label = 3;
    elseif ~isempty(strfind(parentName,'pitted')) || strcmp(parentName,'ps'), label = 4;
    elseif ~isempty(strfind(parentName,'rolled')) || ...
            ~isempty(strfind(parentName,'scale')) || strcmp(parentName,'rs'), label = 5;
    elseif ~isempty(strfind(parentName,'scratch')) || strcmp(parentName,'sc'), label = 6;
    end
end
if label==0, error('Cannot infer class: %s',path); end
end

function tf = neu_has_prefix(name,code)
tf = 0;
if length(name)<length(code) || ~strcmp(name(1:length(code)),code), return; end
if length(name)==length(code), tf = 1; return; end
ch = name(length(code)+1);
tf = ch=='_' || ch=='-' || (ch>='0' && ch<='9');
end

function B = neu_resize_bilinear(A,newRows,newCols)
[oldRows,oldCols] = size(A);
if oldRows==newRows && oldCols==newCols, B=A; return; end
xNew = linspace(1,oldCols,newCols);
yNew = linspace(1,oldRows,newRows)';
[Xnew,Ynew] = meshgrid(xNew,yNew);
B = interp2(1:oldCols,(1:oldRows)',A,Xnew,Ynew,'linear');
if any(~isfinite(B(:))), error('Resize produced non-finite values.'); end
end
