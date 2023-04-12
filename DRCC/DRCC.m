function [F S G residue] = DRCC(X,m1,m2,lambda,mu,ITE,option)

%Dual Regularized Co-Clustering
%written by Quanquan Gu
%2009.2.6
%min ||X-G*S*F'||_F^2+lambda*tr(F'*LF*F)+mu*tr(G'*LG*G)


[d n]=size(X); 

%construct data graph
WF = constructW(X',option);
DF = diag(sum(WF));
LF = DF - WF;

PLF = (abs(LF) + LF)/2;
NLF = (abs(LF) - LF)/2;

%construct feature graph
WG = constructW(X,option);
DG = diag(sum(WG));
LG = DG - WG;

PLG = (abs(LG) + LG)/2;
NLG = (abs(LG) - LG)/2;

%initialized G using Kmeans
res = kmeans(X,m1,'emptyaction','singleton');
G = zeros(d,m1);
for i = 1:d
    G(i,res(i)) = 1;
end
G = G+0.2;

%G = abs(rand(d,m1)); % randomly initialize F

%initialized F using Kmeans
res = kmeans(X',m2,'emptyaction','singleton');
F = zeros(n,m2);
for i = 1:n
    F(i,res(i)) = 1;
end
F = F+0.2;

%F = abs(rand(n,m2)); % randomly initialize G
eps=1e-9; % set your own tolerance

residue = [];
for ite = 1:ITE
    %ite
    % update S
    S = inv(G'*G)*G'*X*F*inv(F'*F);
    
    % update F
    A = X'*G*S;
    PA = (abs(A) + A)/2;
    NA = (abs(A) - A)/2;
    
    B = S'*G'*G*S;
    PB = (abs(B) + B)/2;
    NB = (abs(B) - B)/2;
    
    
    
    F = F.*sqrt((lambda*NLF*F+PA+F*NB)./(lambda*PLF*F+NA+F*PB + eps));
    
    % Renormalize so colloums of F have constant energy
    norms = sqrt(sum(F.^2,1));
    F = F./repmat(norms,n,1);
    S = S.*repmat(norms,m1,1);
    
    %update G
    P = X*F*S';
    PP = (abs(P) + P)/2;
    NP = (abs(P) - P)/2;
    
    Q = S*F'*F*S';
    PQ = (abs(Q) + Q)/2;
    NQ = (abs(Q) - Q)/2;
    
    
    
    G = G.*sqrt((mu*NLG*G+PP+G*NQ)./(mu*PLG*G+NP+G*PQ + eps));
    
    % Renormalize so colloums of G have constant energy
    norms = sqrt(sum(G.^2,1));
    G = G./repmat(norms,d,1);
    S = S.*repmat(norms',1,m2);

        
    %compute residue
    tmp = lambda*trace(F'*LF*F)+ mu*trace(G'*LG*G)+ sum(sum((X-G*S*F').^2)) ;
    residue = [residue tmp];

end

end


