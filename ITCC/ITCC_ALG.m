function [row_cluster_label, column_cluster_label, time] = ITCC_ALG(data, row_class_number, column_class_number)
    %% Running
    tic
    [row_cluster_label, column_cluster_label] = ITL(data, row_class_number, column_class_number);
    time = toc;
end