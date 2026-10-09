clear all 
clc
 load('COIL20.mat');	
 nClass = length(unique(gnd));
 dataset='COIL20';      
 fea=double(fea);
 
 %Sort samples by labels
 [labels,index]=sort(gnd,'ascend');
 gnd=labels;
 fea=fea(index,:);
 fea=double(fea);
 
 %Unit Euclidean length normalization
 fea = NormalizeFea(fea); 

 LabelsRatio=0.2;
 K=nClass; %full-size dataset as input data  
 meanAC=[];
 meanMI=[];
 
  TempAC=[];
  TempMI=[];
  
  % hyperparameter  alpha setting
  %datasets 'Yale','PIE',   'MSRA25', 'AR',  'COIL20',  'COIL100', 'Optdigits','USPS10'
     %EDDNMF   1	 50	        50	   0.1	   0.01	        10	       0.001	0.1
     %ERDNMF   10	 400	    10	   200	    800	        1000	     600	600
 
  
  for h=1:20
     %Select the data points of K classes as the input of the algorithm
     [X,Smpgnd,count]=CreatSampleDatasets(fea,K,gnd,nClass,LabelsRatio);
     % X:selected data points belong to K Classes
     % Smpgnd: the lables of selected data points
     % count: the number of selected data points
     
     Options.maxIter=200;
     Options.gndSmpNum=count;
     Options.alpha=0.01;
     Options.labels=Smpgnd;
     Options.KClass=K;
     [~,V]=EDDNMF(X',Options);
     [~, label] = max(V');
    % label = litekmeans(V,KClass,'Replicates',20);
     newL=bestMap(Smpgnd,label);
     AC=Accuracy(newL,Smpgnd);
     MIhat = MutualInfo(Smpgnd,label);  
     TempAC(1,h)=AC;
     TempMI(1,h)=MIhat;


     %seminmf_div
     Options.maxIter=200;
     Options.alpha=800;
     Options.gndSmpNum=count;
     Options.labels=Smpgnd;
     Options.KClass=K;
    [~, V] =ERDNMF(X',Options);
    [~, label] = max(V');
    % label = litekmeans(V',KClass,'Replicates',20);
     newL=bestMap(Smpgnd,label);
     AC=Accuracy(newL,Smpgnd);
     MIhat = MutualInfo(Smpgnd,label);  
     TempAC(2,h)=AC;
     TempMI(2,h)=MIhat;

     
 end
    meanAC=mean(TempAC,2);
    meanMI=mean(TempMI,2);
  