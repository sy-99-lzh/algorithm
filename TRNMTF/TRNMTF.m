function [UU, SS, VV, OBJ] = TRNMTF(X, WU, WV, c1, c2, lambda, mu, T,U,V,labels)
%
% Graph Dual Regularized Tri-Nonnegative Matrix Factorization
% Written by Fanhua Shang
% Key Laboratory of Intelligent Perception and Image Understanding of Ministry of Education of China,
% Xidian University, Xi'an 710071, China.
%
% ||X-U*S*V'||^2+lambda*tr(U'*(DU-WU)*U)+mu*tr(V'*(DV-WV)*V)
%
% X:         NxM input data matrix
% WU:        NxN adjacency matrix of the data graph
% WV:        MxM adjacency matrix of the feature graph
% c1,c2:     number of hidden factors
% lambda,mu: Regularization parameters
% T:         number of iterations


% References: 
% Fanhua Shang, Licheng Jiao, and Fei Wang,
% Graph Dual Regularization Non-negative Matrix Factorization 
% for Co-clustering, Pattern Recognition, 45(6): 2237-2250, 2012.

if ~isempty(WU)
    DU = diag(sum(WU));
    LU = DU - WU;
    LU_P = (abs(LU) + LU)/2;
    LU_N = (abs(LU) - LU)/2;
end

if ~isempty(WV)
    DV = diag(sum(WV));
    LV = DV - WV;
    LV_P = (abs(LV) + LV)/2;
    LV_N = (abs(LV) - LV)/2;
end

[M N] = size(X); 
% U = abs(randn(M, c1));
% V = abs(randn(N, c2));
S = abs(randn(c1, c2));


feaSum = full(sum(X, 2));
D_half = (X'*feaSum).^.5;
for i = 1: N
    X(:,i) = X(:,i)/D_half(i);
end

eps = 1e-4; 

% ALPHA_W=0.9;
% % V
% ALPHA_H=0.9;
% 
% OP=0.001; 
% % V
% OX=0.001;  
%ALPHA_W:L2U;ALPHA_H:L2V;OP:L1U;OX:L1V
% ALPHA_W=10^-xx;
% ALPHA_H=10^-xx;
% OP=10^-yy; 
% OX=10^yy;
ALPHA_W = 0.0001;%U的2范数
ALPHA_H = 0.0001;%V的2范数
OP=0.00005;%U的1范数
OX=0.00005;%V的1范数





c = length(unique(labels));
obj = [];
for iter = 1:T
    
    iter;
    % ===================== update S ========================
    S = S.*(U'*X*V)./(U'*(U*S*V')*V);
    SS{iter} = S;
    
    % ===================== update U ========================
    U = U.*(X*V*S' + lambda*LU_N*U)./((U*S*V')*V*S' + lambda*LU_P*U +ALPHA_W*U+ eps);

         for i=1:size(U,1)
             for j=1:size(U,2)
                  if(U(i,j)>OP)
                     U(i,j)=U(i,j)-OP;      
                  else
                     U(i,j)=0; 
                 end
             end
         end 
    % Normalization of U
    norms = sqrt(sum(U.^2, 1));
    U = U./repmat(norms, M, 1);
    S = S.*repmat(norms', 1, c2);
    
    UU{iter} = U;
    % ===================== update V ========================
    V = V.*((X)'*U*S + mu*LV_N*V)./((U*S*V')'*U*S + mu*LV_P*V +ALPHA_H*V+ eps);
    
%     comptime{iter} = toc;
    
%     tic
    label_kmeans_TRNMTF = litekmeans(V, c, 'Replicates', 20);
    
    
%     LABEL{iter} = label_kmeans_TRNMTF;
%     [kmeans_results_TRNMTF_nmi{iter}, kmeans_results_TRNMTF_purity{iter}, kmeans_results_TRNMTF_acc{iter}]=evaluation_index(labels,label_kmeans_TRNMTF);
%     [AR_TRNMTF{iter},RI_TRNMTF{iter},MI_TRNMTF{iter},HI_TRNMTF{iter}]=RandIndex(labels,label_kmeans_TRNMTF);
%     Index = {kmeans_results_TRNMTF_nmi;kmeans_results_TRNMTF_purity;kmeans_results_TRNMTF_acc;AR_TRNMTF;RI_TRNMTF;MI_TRNMTF;HI_TRNMTF};
    
    for i=1:size(V,1)
         for j=1:size(V,2)
             if(V(i,j)>OX)
                 V(i,j)=V(i,j)-OX; 
             else
                 V(i,j)=0;      
             end         
         end
    end
     
    % Normalization of V
    norms = sqrt(sum(V.^2, 1));
    V = V./repmat(norms, N, 1);
    S = S.*repmat(norms, c1, 1);
   
    VV{iter} = V;
    % The objective function
    temp = sum(sum((X - U*S*V').^2)) + lambda*trace(U'*LU*U) + mu*trace(V'*LV*V) + OP*norm(U,1) + OX*norm(V',1) + ALPHA_W*norm(U,2) + ALPHA_H*norm(V',2) ;
    if temp < eps
        break;
    end
     
    obj = [obj, temp];
    OBJ{iter} = obj;
end