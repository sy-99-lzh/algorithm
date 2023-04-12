function [r,c]=Bregcc(data,k,l);
%load yeast.mat
[m,n]=size(data);
W=ones(m,n);
[R,C]=coclust_euc(data,W,k,l,6,1e-3);
for i=1:k
    R(:,i)=R(:,i)*i;
end
r=sum(R,2);
for j=1:l
    C(:,j)=C(:,j)*j;
end
c=sum(C,2);

    



    