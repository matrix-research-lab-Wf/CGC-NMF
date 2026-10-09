 function   [U,V]=RationNMF(X,Options)
      
 %..............................................
       [M, N]= size(X);

%........................................................
        Data_num = length(Options.labels);
        n=Options.gndSmpNum;
        labels=Options.labels;
        labels(n+1:end)=0;
        label_matrix = zeros(Data_num, Options.KClass);
 %...........................................       
        for i =1:Data_num
            if labels(i)>0
              label_matrix(i,labels(i)) = 1;
            end
        end
        
        Q =label_matrix;
        Q1=abs(Q-1);
        Q1(n+1:end,:) = 0;
        Q=Q';
        Q1= Q1';
        U = rand(M, Options.KClass);
        V = rand(N, Options.KClass);
        for iters = 1:Options.maxIter
          U =U.*(X*V)./(U*V'*V+eps);
          F1=X'*U+Options.alpha* (sum(sum(Q1.*V'))/sum(sum(Q.*V'))^2)*(Q'); 
          F2=V*U'*U +Options.alpha*(1/sum(sum(Q.*V')))*(Q1');
          V =V.*(F1./ (F2+eps) );%             
        end
end
       
    