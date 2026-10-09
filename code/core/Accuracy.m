function [AC]=Accuracy(label,gnd)

% Compute clustering accuracy after aligning predicted and true labels.
% label: predicted cluster labels.
% gnd: ground-truth class labels.
% 
[m,n]=size(label);
[M,N]=size(gnd);
delta=0;
if M==m   
    for i=1:m
        if label(i,1)==gnd(i,1)
            delta=delta+1;
        end
    end
else
    error('The predicted and ground-truth label vectors have different lengths.');
end
AC=delta/m;
