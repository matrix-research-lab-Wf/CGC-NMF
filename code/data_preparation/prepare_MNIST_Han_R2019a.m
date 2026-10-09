function prepare_MNIST_Han_R2019a(workRoot,outputMat,indexCsv)
%PREPARE_MNIST_HAN_R2019A Download official MNIST and reproduce MNIST_Han.
% The committed index list was recovered by exact pixel-and-label matching
% against the author-side matrix; it contains no MNIST image data.

if nargin<1 || isempty(workRoot), workRoot=fullfile(pwd,'official_raw','mnist'); end
if nargin<2 || isempty(outputMat), outputMat=fullfile(pwd,'MNIST_Han.mat'); end
if nargin<3 || isempty(indexCsv), indexCsv=fullfile(fileparts(mfilename('fullpath')),'mnist_han_indices.csv'); end
if exist(workRoot,'dir')~=7, mkdir(workRoot); end

names={'train-images-idx3-ubyte','train-labels-idx1-ubyte', ...
       't10k-images-idx3-ubyte','t10k-labels-idx1-ubyte'};
primary='http://yann.lecun.com/exdb/mnist/';
authorizedMirror='https://storage.googleapis.com/cvdf-datasets/mnist/';
for i=1:numel(names)
    raw=fullfile(workRoot,names{i}); gz=[raw '.gz'];
    if exist(raw,'file')~=2
        if exist(gz,'file')~=2
            try, websave(gz,[primary names{i} '.gz']);
            catch, warning('Primary MNIST host unavailable; using the CVDF authorized mirror.'); websave(gz,[authorizedMirror names{i} '.gz']); end
        end
        gunzip(gz,workRoot);
    end
end
[X,y]=read_all(workRoot,names);
indices=dlmread(indexCsv);
if numel(indices)~=6996 || any(indices<1 | indices>70000 | indices~=round(indices))
    error('Invalid MNIST_Han index list: %s.',indexCsv);
end
fea=double(X(indices,:)); gnd=double(y(indices))+1;
if ~isequal(size(fea),[6996 784]), error('Unexpected MNIST_Han dimensions.'); end
source_url='http://yann.lecun.com/exdb/mnist/'; %#ok<NASGU>
selection_index=indices(:); %#ok<NASGU>
save(outputMat,'fea','gnd','source_url','selection_index','-v7');
fprintf('Output: %s\nMNIST_HAN_PREPARATION_PASS=1\n',outputMat);
end

function [X,y]=read_all(root,n)
X=[read_images(fullfile(root,n{1}));read_images(fullfile(root,n{3}))];
y=[read_labels(fullfile(root,n{2}));read_labels(fullfile(root,n{4}))];
end
function X=read_images(path)
f=fopen(path,'rb','ieee-be'); if f<0,error('Cannot open %s.',path);end; c=onCleanup(@()fclose(f));
magic=fread(f,1,'uint32'); n=fread(f,1,'uint32'); r=fread(f,1,'uint32'); q=fread(f,1,'uint32');
if magic~=2051 || r~=28 || q~=28,error('Invalid MNIST image file.');end
v=fread(f,double(n*r*q),'*uint8'); if numel(v)~=double(n*r*q),error('Truncated MNIST image file.');end
X=reshape(v,double(r*q),double(n))';
end
function y=read_labels(path)
f=fopen(path,'rb','ieee-be'); if f<0,error('Cannot open %s.',path);end; c=onCleanup(@()fclose(f));
magic=fread(f,1,'uint32'); n=fread(f,1,'uint32'); if magic~=2049,error('Invalid MNIST label file.');end
y=fread(f,double(n),'*uint8'); if numel(y)~=double(n),error('Truncated MNIST label file.');end
end
