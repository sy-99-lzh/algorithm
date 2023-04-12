function [row_cluster_label, column_cluster_label, time] = EDVWs_ALG(data, row_class_number)
    %% Running
    tic
    alg = 'alg1'; 
    alpha = 0.5;
    [row_cluster_label, column_cluster_label] = edvw(data, row_class_number, alg, alpha);
    time = toc;
end