function [row_cluster_label, column_cluster_label, time] = BCC_ALG(data, row_class_number, column_class_number)
    %% Running
    tic
    [row_cluster_label, column_cluster_label] = Bregcc(data, row_class_number, column_class_number);
    time = toc;
end    