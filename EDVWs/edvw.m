function [row_labels, col_labels] = edvw(X, k, alg, alpha)
    normalized_features = star_laplacian(X, k, alg, alpha);
    idx = kmeans(normalized_features, k);
    [row, col] = size(X);
    row_labels = idx(1: row, :);
    col_labels = idx( row + 1: row + col, :);
end

function VV = star_laplacian(X, k, alg, alpha)
    R = X';
    hyperedge_weights = comp_hyperedge_weights(R);
    W = comp_W(R, hyperedge_weights);
    P_VE = comp_P_VE(W);
    P_EV = comp_P_EV(R);
    P = comp_P(P_VE, P_EV);
    P_alpha = comp_P_alpha(P, alpha);
    pi = comp_pi(P_alpha);
    A_bar = comp_A_bar(P_VE, P_EV, pi);
    V = comp_V(A_bar, k);
   if strcmp(alg, 'alg2')
       VV = alg2(V);
   else
       VV = alg1(V, pi);
   end
end

function hyperedge_weights = comp_hyperedge_weights(R)
    hyperedge_weights = std(R, 0, 2)';
end

function W = comp_W(R, hyperedge_weights)
    [~, num_nodes] = size(R);
    hw_mat = repmat(hyperedge_weights, num_nodes, 1);
    RR = R';
    RR(RR <= 0) = 0;
    RR(RR > 0) = 1;
    W = RR .* hw_mat;
end

function P_VE = comp_P_VE(W)
    P_VE = W ./ sum(W, 2);
end

function P_EV = comp_P_EV(R)
    P_EV = R ./ sum(R, 2);
end

function P = comp_P(P_VE, P_EV)
    [num_edges, num_nodes] = size(P_EV);
    P = [zeros(num_nodes), P_VE; P_EV, zeros(num_edges)];
end

function P_alpha = comp_P_alpha(P, alpha)
    P_alpha = (1-alpha) * eye(size(P, 1)) + alpha * P;
end

function pi = comp_pi(P)
    [~, ~, W] = eig(P);
    pi = real(W(:,1));
    pi = pi / sum(pi);
end

function A_bar = comp_A_bar(P_VE, P_EV, pi)
    [num_edges, num_nodes] = size(P_EV);
    pi_1 = pi(1: num_nodes);
    pi_2 = pi(num_nodes + 1: num_nodes + num_edges);
    A_bar = 0.5 * (diag(pi_1 .^ 0.5) * P_VE * diag(pi_2 .^ 0.5) + diag(pi_1 .^ 0.5) * P_EV' * diag(pi_2 .^ 0.5));
end

function V = comp_V(A_bar, k)
    [U_bar, ~, V_bar] = svd(A_bar);
    V_V = U_bar(:, 1: k);
    V_E = V_bar(:, 1: k);
    V = real([V_V; V_E]);
end

function VV = alg1(V, pi) 
    VV = diag(pi .^ -0.5) * V;
end

function VV = alg2(V)
    VV = V ./ (sum(V .* V, 2));
end