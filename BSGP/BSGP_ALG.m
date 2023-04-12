function [row_cluster_label, column_cluster_label, time] = BSGP_ALG(data, row_class_number)
    %% running
    tic
    [row_cluster_label, column_cluster_label] = SpectralCoClustering(data, row_class_number);
    time = toc;
end