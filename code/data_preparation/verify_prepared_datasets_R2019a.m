function verify_prepared_datasets_R2019a(dataRoot,referenceRoot)
%VERIFY_PREPARED_DATASETS_R2019A Validate shapes, labels and optional equality.
if nargin<1||isempty(dataRoot),dataRoot=fullfile(pwd,'data');end
if nargin<2,referenceRoot='';end
spec={ 'main_six','PIE.mat',2856,1024,68; 'main_six','YaleB.mat',2414,1024,38; ...
 'main_six','COIL20_Obj.mat',1440,1024,20; 'main_six','COIL100_Obj.mat',7200,1024,100; ...
 'main_six','MNIST_Han.mat',6996,784,10; 'neu_cls','NEU_CLS_32x32.mat',1800,1024,6; ...
 'neu_cls','NEU_CLS_LBP59_4x4.mat',1800,944,6};
for i=1:size(spec,1)
 p=fullfile(dataRoot,spec{i,1},spec{i,2});if exist(p,'file')~=2,fprintf('MISSING %s\n',p);continue;end
 S=load(p,'fea','gnd');ok=isequal(size(S.fea),[spec{i,3} spec{i,4}])&&numel(S.gnd)==spec{i,3}&&numel(unique(S.gnd))==spec{i,5}&&all(isfinite(S.fea(:)))&&all(S.fea(:)>=0);
 fprintf('%s STRUCTURE_%s\n',spec{i,2},passfail(ok));if ~ok,error('Validation failed for %s.',p);end
 if ~isempty(referenceRoot)
  q=find_reference(referenceRoot,spec{i,2});if ~isempty(q),R=load(q,'fea','gnd');exact=isequal(S.fea,R.fea)&&isequal(S.gnd(:),R.gnd(:));fprintf('%s AUTHOR_EQUALITY_%s\n',spec{i,2},passfail(exact));end
 end
end
end
function s=passfail(x),if x,s='PASS';else,s='FAIL';end,end
function q=find_reference(root,name),q='';d=dir(fullfile(root,'**',name));if ~isempty(d),q=fullfile(d(1).folder,d(1).name);end,end
