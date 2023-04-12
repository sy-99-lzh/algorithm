function [row_cluster_label ,column_cluster_label, time] = DRCC_ALG(data, row_class_number, column_class_number, row_size, column_size)
    
    %% parameter
    tic
    [row_size, column_size] = size(data);
    ITE = 100;
    eps = 1e-9;
    norms = sqrt(sum(data.^2, 2));
    DRCC_data = data./(repmat(norms, 1, column_size)+eps);
    K = 5;
    lambda = 0.01;
    mu = 0.01;

    option.Metric = 'Euclidean';
    option.NeighborMode = 'KNN';
    option.k = K;
    option.WeightMode = 'Binary';
    option.t = 1;
    
    %% running
    tic
    [F, S, G, residue] = DRCC(DRCC_data', row_class_number , column_class_number, lambda, mu, ITE, option);
    time = toc;

    %% predict label
    row_cluster_label = zeros(row_size,1);
    for num = 1:row_size
        [tmp, row_cluster_label(num)] = max(F(num,:));
    end
    column_cluster_label = zeros(column_size,1);
    for num = 1:column_size
        [tmp, column_cluster_label(num)] = max(G(num,:));
    end
end
