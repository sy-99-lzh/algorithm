function [row_cluster_label, column_cluster_label, time] = SOBG_ALG(data, row_class_number)
    %% Running
    tic
    [row_cluster_label, column_cluster_label, S] = coclustering_bipartite_fast(data, row_class_number);
    time = toc;
end  