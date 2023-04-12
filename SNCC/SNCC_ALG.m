function [row_cluster_label, column_cluster_label, time] = SNCC_ALG(data, row_class_number, column_class_number)
    %% Running
    tic
    SNCC_data = data;
    for SNCC_num = 1:size(SNCC_data)
        SNCC_data(SNCC_num,:) = SNCC_data(SNCC_num,:) ./ max(1e-12, norm(SNCC_data(SNCC_num,:))); %Normalize
    end
    [row_cluster_label, column_cluster_label] = SNCC(SNCC_data, row_class_number, column_class_number);
    time = toc;
end