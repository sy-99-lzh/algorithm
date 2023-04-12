function [y1,y2] = FBGPC(X, k, l, iter ,sparse_ratio,pruning)
%Code for the following paper:
%Wei Chen, Hongjun Wang, Zhiguo Long, Tianrui Li.
%Fast Flexible Bipartite Graph Model for Co-Clustering 
%IEEE Transaction on Knowledge and Data Engineering

%input��
%   X��sample dataset
%   k:  number of row clusters
%   l:   number of column clusters
%   iter: number of inflation iterations. And the default is 10
%   pruning: Pruning parameter
%   sparse_ratio: Termination condition parameter

%% Initialization 
if ~exist('l', 'var') || isempty(l)
    l = k;
end
if ~exist('iter', 'var') || isempty(iter)
    iter = 10;
end
if ~exist('sparse_ratio', 'var') || isempty(sparse_ratio)
    sparse_ratio = 0.8;
end
if ~exist('pruning', 'var') || isempty(pruning)
    pruning = 0.001;
end

[N,M] = size(X);
count_num = N*M;


%% Construct Bipartite Graph
g=X;                                                    %Bipartite Graph
i = 1;                                                    %Iteration
g = bsxfun(@rdivide, g, max(g,[],1));       %normalization 

%% Inflation Operation
while i < iter
    g1 = g;
    zeros_flags = g1==0;
    zeros_ratio = sum(zeros_flags(:))/count_num;
    
    if zeros_ratio>=sparse_ratio
        break;
    end
    
    %Determine Sparsity
    g1(g1 < pruning) = 0;
    
    %Inflation
    inflation_g = g1.*g1;
    grow = bsxfun(@rdivide, inflation_g, sum(inflation_g,1));
    gcol = bsxfun(@rdivide, inflation_g, sum(inflation_g,2));
    g2 = 0.5*grow+0.5*gcol;
    g2 = bsxfun(@rdivide, g2, max(g,[],1));   
    
    g=g2;
    i = i + 1;
end

%% Zero-one normalization
g(g~=0) = 1;

%% Extracting Co-Occurrence Structure and Co-Clustering Results
distance_row = pdist(g, 'euclidean');
cluster_row = linkage(distance_row, 'ward');
id_row = cluster(cluster_row, 'maxclust', k);



distance_col = pdist(g', 'euclidean');
cluster_col = linkage(distance_col, 'ward');
id_col = cluster(cluster_col, 'maxclust', l);

y1 = id_row;
y2 = id_col;
end

