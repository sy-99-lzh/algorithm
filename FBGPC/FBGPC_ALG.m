function [row_cluster_label ,column_cluster_label, time] =  FBGPC_ALG(data, row_class_number, column_class_number)
%% parameter
    tic
    stop_iter = 10;
    sparse_ratio = 0.5;
    pruning = 0.001;
    
 %% Running   
    tic
    [row_cluster_label,column_cluster_label] = FBGPC(data, row_class_number, column_class_number, stop_iter, sparse_ratio, pruning);
    time = toc;
end