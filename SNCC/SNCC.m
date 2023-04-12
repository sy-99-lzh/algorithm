function [FClu, GClu] = SNCC(X, c1, c2)

% Set parameters
eps = 1e-10; 
maxIte = 50;
alpha = 0.1;
beta = 0.1;

% Construct W1
options = [];
options.NeighborMode = 'KNN';
options.k = 10;
% options.WeightMode = 'Binary';
% options.WeightMode = 'Cosine';options.bNormalized=1;
options.WeightMode = 'HeatKernel';
W1 = constructW(X,options);
% dd = 1 ./ sum(W1);
% dd = sqrt(dd);
% W1 = bsxfun(@times, W1, dd);
% W1 = W1';
% W1 = bsxfun(@times, W1, dd);
% W1 = (W1 + W1') / 2;

% Construct W2
options = [];
options.NeighborMode = 'KNN';
options.k = 10;
% options.WeightMode = 'Binary';
% options.WeightMode = 'Cosine';options.bNormalized=1;
options.WeightMode = 'HeatKernel';
W2 = constructW(X',options);
% dd = 1 ./ sum(W2);
% dd = sqrt(dd);
% W2 = bsxfun(@times, W2, dd);
% W2 = W2';
% W2 = bsxfun(@times, W2, dd);
% W2 = (W2 + W2') / 2;

% Initialization
[d, n] = size(X); 
idx1 = litekmeans(X, c1, 'MaxIter', 100);
idx2 = litekmeans(X', c2, 'MaxIter', 100);
% idx1 = kmeansPP(X', c1); 
% idx2 = kmeansPP(X, c2); 
F = sparse(1:d, idx1, ones(d,1), d, c1) + 0.2;
G = sparse(1:n, idx2, ones(n,1), n, c2) + 0.2;


S = (F'*F)\F'*X*G/(G'*G);

% Update F, S, G, Z1, Z2
for i = 1:maxIte
    %i
    % Z1, Z2
    Z1 = W1'*F /(F'*F+eps); 
    Z2 = W2'*G /(G'*G+eps);
    % F, M, N 
    M = W1*Z1;
    N = Z1'*Z1;
    F = F.*sqrt((X*G*S'+alpha*(0.5*(sign(abs(M))+sign(M)).*M+F*(0.5*(sign(N)-sign(abs(N))).*N)))./...
        (F*S*(G'*G)*S'+alpha*(0.5*(sign(M)-sign(abs(M))).*M+F*(0.5*(sign(abs(N))+sign(N)).*N))+eps)); 

    % G, P, Q
    P = W2*Z2;
    Q = Z2'*Z2;
    G = G.*sqrt((X'*F*S+beta*(0.5*(sign(abs(P))+sign(P)).*P+G*(0.5*(sign(Q)-sign(abs(Q))).*Q)))./...
        (G*S'*(F'*F)*S+beta*(0.5*(sign(P)-sign(abs(P))).*P+G*(0.5*(sign(abs(Q))+sign(Q)).*Q))+eps)); 

    % S    
    S = S.*sqrt((F'*X*G)./(F'*F*S*(G'*G)+eps));

end

% Generate labels
[~,FClu] = max(F');
[~,GClu] = max(G');

end