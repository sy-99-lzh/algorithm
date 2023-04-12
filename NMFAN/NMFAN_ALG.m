function [row_cluster_label, column_cluster_label, time] = NMFAN_ALG(data, row_class_number, column_class_number)
    %% Running
    tic
    k = 5;% k-NN, average number of neighbors 
    lambda = 100;
    ITE = 50 ;
    [V, U, objV] = NMFAN(data', row_class_number, k, lambda, ITE);
    row_cluster_label = kmeans(V,row_class_number);
    column_cluster_label = kmeans(U, column_class_number);
    time = toc;
end
