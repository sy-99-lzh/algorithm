function [row_cluster_label, column_cluster_label, time] = TRNMTF_ALG(data, row_class_number, column_class_number)
    %% Running
    tic
    [row_size, column_size] = size(data);
    iter = 300;
    V = abs(randn(row_size,row_class_number));
    U = abs(randn(column_size,column_class_number));
    
    [U, S, V, obj] = run_TRNMTF(data, row_class_number, U, V, iter);
    
    V_data = V{1, iter} + 0.1;
    U_data = U{1, iter} + 0.1;
    [row_cluster_label, ~, ~, ~] = kmeans(V{1, iter}, row_class_number);
    [column_cluster_label, ~, ~, ~] = kmeans(U{1, iter}, column_class_number);
    time = toc;
end  