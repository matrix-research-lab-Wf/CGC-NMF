function indices = derive_MNIST_Han_indices_R2019a(authorMat, rawRoot, outputCsv)
%DERIVE_MNIST_HAN_INDICES_R2019A Recover the public MNIST row indices.
% This author-side utility matches every row of MNIST_Han.mat against the
% official 60,000 training and 10,000 test images. The resulting index list
% contains no image data and may be distributed with the reproducibility
% package. Indices 1:60000 refer to training images and 60001:70000 to test.

if nargin < 1 || isempty(authorMat), error('authorMat is required.'); end
if nargin < 2 || isempty(rawRoot), error('rawRoot is required.'); end
if nargin < 3 || isempty(outputCsv)
    outputCsv = fullfile(fileparts(mfilename('fullpath')),'mnist_han_indices.csv');
end

S = load(authorMat,'fea','gnd');
if ~isfield(S,'fea') || ~isfield(S,'gnd') || ~isequal(size(S.fea),[6996 784])
    error('Expected author matrix with fea 6996-by-784 and gnd 6996-by-1.');
end
[allFea,allLabels] = read_mnist_pair(rawRoot);
target = uint8(S.fea);
targetLabels = double(S.gnd(:))-1;
if any(double(target(:)) ~= S.fea(:)) || any(targetLabels < 0 | targetLabels > 9)
    error('Author matrix must contain uint8-valued pixels and labels 1..10.');
end

fprintf('Indexing 70,000 official MNIST images...\n');
lookup = containers.Map('KeyType','char','ValueType','double');
for i = 1:size(allFea,1)
    key = row_key(allFea(i,:),allLabels(i));
    if ~isKey(lookup,key), lookup(key)=i; end
end

indices = zeros(size(target,1),1);
for i = 1:size(target,1)
    key = row_key(target(i,:),targetLabels(i));
    if ~isKey(lookup,key)
        % Some MAT files store each 28-by-28 image after MATLAB column-wise
        % vectorization. Try that documented alternative before failing.
        transposed = reshape(reshape(target(i,:),28,28)',1,784);
        key = row_key(transposed,targetLabels(i));
    end
    if ~isKey(lookup,key), error('No official MNIST match for author row %d.',i); end
    indices(i)=lookup(key);
    if mod(i,500)==0 || i==size(target,1), fprintf('  matched %d/%d\n',i,size(target,1)); end
end

dlmwrite(outputCsv,indices,'precision','%d');
fprintf('Wrote %s\n',outputCsv);
fprintf('MNIST_INDEX_DERIVATION_PASS=1\n');
end

function key = row_key(row,label)
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(row(:)),'int8'));
d = typecast(md.digest(),'uint8');
key = sprintf('%d_%02x',label,d);
end

function [fea,labels] = read_mnist_pair(rootDir)
trainImages = ensure_unzipped(rootDir,'train-images-idx3-ubyte');
trainLabels = ensure_unzipped(rootDir,'train-labels-idx1-ubyte');
testImages = ensure_unzipped(rootDir,'t10k-images-idx3-ubyte');
testLabels = ensure_unzipped(rootDir,'t10k-labels-idx1-ubyte');
fea = [read_idx_images(trainImages); read_idx_images(testImages)];
labels = [read_idx_labels(trainLabels); read_idx_labels(testLabels)];
end

function path = ensure_unzipped(rootDir,baseName)
path = fullfile(rootDir,baseName);
if exist(path,'file')==2, return; end
gzPath = [path '.gz'];
if exist(gzPath,'file')~=2, error('Missing %s or %s.',path,gzPath); end
gunzip(gzPath,rootDir);
end

function X = read_idx_images(path)
fid=fopen(path,'rb','ieee-be'); if fid<0, error('Cannot open %s.',path); end
c=onCleanup(@() fclose(fid));
magic=fread(fid,1,'uint32'); n=fread(fid,1,'uint32'); r=fread(fid,1,'uint32'); q=fread(fid,1,'uint32');
if magic~=2051 || r~=28 || q~=28, error('Invalid MNIST image file: %s.',path); end
raw=fread(fid,double(n)*double(r)*double(q),'*uint8');
if numel(raw)~=double(n)*double(r)*double(q), error('Truncated MNIST image file.'); end
X=reshape(raw,double(r)*double(q),double(n))';
end

function y = read_idx_labels(path)
fid=fopen(path,'rb','ieee-be'); if fid<0, error('Cannot open %s.',path); end
c=onCleanup(@() fclose(fid));
magic=fread(fid,1,'uint32'); n=fread(fid,1,'uint32');
if magic~=2049, error('Invalid MNIST label file: %s.',path); end
y=fread(fid,double(n),'*uint8');
if numel(y)~=double(n), error('Truncated MNIST label file.'); end
end
