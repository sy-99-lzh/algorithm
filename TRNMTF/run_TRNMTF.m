function [U_T, S_T, V_T, obj_TRNMTF]=run_TRNMTF(data,labels, U, V, Iter)
%     addpath ../evaluation_index
%     addpath ../classic_clustering_algoritms/AP_clustering
%     addpath ../classic_clustering_algoritms/DP_clustering

    c = length(unique(labels));
    c1=c;
    c2=c;

    options = [];
    options.Metric = 'Euclidean';
    options.NeighborMode = 'KNN';
    options.k =4;
    options.WeightMode = 'Binary';
    options.t = 1; 

    WV = constructW(data, options);
    for i=1:size(WV,1)
        WV(i,i)=0;
    end

    options = [];
    options.Metric = 'Euclidean';
    options.NeighborMode = 'KNN';
    options.k =4;
    options.WeightMode = 'Binary';
    options.t = 1;

    WU = constructW(data',options);
    for i=1:size(WU,1)
        WU(i,i)=0;
    end

    %% % DNMF
%     Iter = 5;
    lambda = 0.1;
    mu = lambda;

    % [U, V, obj] = DNMF(data', WU, WV, c, lambda,mu, Iter);
%     [U, V, obj_DNMF] = DNMF(data', WU, WV, c, lambda,mu, Iter, U, V);
    [U_T, S_T, V_T, obj_TRNMTF] = TRNMTF(data', WU, WV, c1, c2, lambda, mu, Iter,U,V,labels);
    
%     figure (1)
%     subplot(2,2,3)
%     plot(log(obj_DNMTF)) 
%     title('DNMTF')
%     xlabel ('Number of iterations')
%     ylabel ('Log objective')





end
