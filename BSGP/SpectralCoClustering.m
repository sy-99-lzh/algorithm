function [row_clust_idx, col_clust_idx]=SpectralCoClustering(A,k,display,names)
%原为：[row_clust_idx,
%col_clust_idx,y_index,x_index,clust_idx]=SpectralCoClustering(A,k,display,names) 
[m,n]=size(A);
%W=ones(m,n);
%b=zeros(m,1);
A(find(A==0))=0.1;
if (isempty(A))
    error('The input matrix is empty');
end;
if (length(find(sum(abs(A),1)==0))>0)
    error('data matrix contains empty features. Please remove them and re-run');
end;
if (length(find(sum(abs(A),2)==0))>0)
    error('data matrix contains empty instances. Please remove them and re-run');
end;
if nargin<2
    k=2;
end;
if nargin<3
    display=0;
end;
D1=diag(sum(A,2));
D2=diag(sum(A,1));
D1(find(D1==0))=eps;
D2(find(D2==0))=eps;
D1_root=abs((D1)^(-0.5));
D2_root=abs((D2)^(-0.5));
An=D1_root*A*D2_root;
[U,S,V]=svd(An,'econ');
num_of_pcs=ceil(log2(k));
z=[D1_root*U(:,2:num_of_pcs+1);D2_root*V(:,2:num_of_pcs+1)];
%z=[D1_root*U(:,2:k);D2_root*V(:,2:k)];
clust_idx=kmeansp(z,k);     %原本为：clust_idx=kmeans(z,k);
row_clust_idx=clust_idx(1:size(A,1));
col_clust_idx=clust_idx(size(A,1)+1:end);
x_index=[];
y_index=[];
for clust=1:k
    clustk_y_idx=find(clust_idx(1:size(A,1))==clust);
    clustk_x_idx=find(clust_idx(size(A,1)+1:end)==clust);
    data_clust=A(clustk_y_idx,clustk_x_idx);
    if (length(clustk_x_idx)>2 & length(clustk_y_idx)>2 & (sum(var(data_clust,0,1))>0 | sum(var(data_clust,0,2))>0))
        if (length(find(sum(data_clust,1)==0))>0 | length(find(sum(data_clust,2)==0))>0) 
            % this section is actually obsolete, since no empty rows or features are allowed
            data_clust_y_zeros=find(sum(data_clust,2)==0);
            data_clust_x_zeros=find(sum(data_clust,1)==0)';
            data_clust_y_no_zeros=find(sum(data_clust,2));
            data_clust_x_no_zeros=find(sum(data_clust,1))';
            data_clust_no_zeros=data_clust(data_clust_y_no_zeros,data_clust_x_no_zeros);
            [clust_x_sorted,clust_y_sorted]=ArrangeData2Clusters(data_clust_no_zeros,'cosine','single',0);
            x_index=[x_index;clustk_x_idx([data_clust_x_no_zeros(clust_x_sorted);data_clust_x_zeros])];
            y_index=[y_index;clustk_y_idx([data_clust_y_no_zeros(clust_y_sorted);data_clust_y_zeros])];
        else
            [clust_x_sorted,clust_y_sorted]=ArrangeData2Clusters(data_clust,'cosine','single',0);
            x_index=[x_index;clustk_x_idx(clust_x_sorted')];
            y_index=[y_index;clustk_y_idx(clust_y_sorted')];
        end;
    else
        x_index=[x_index;clustk_x_idx];
        y_index=[y_index;clustk_y_idx];
    end
end

% displaying the biclusters
if display
    figure;colormap(gray);
    imagesc(-A(y_index,x_index));
    set(gca,'FontSize',14,'FontWeight','bold');
    if nargin>3
        if length(names)==length(y_index)
            names_clusters=names;
            for k=1:length(names)
                name_cluster=[names{k},' (',num2str(clust_idx(k)),')'];
                names_clusters{k}=name_cluster;
            end
            set(gca,'YTick',1:length(y_index),'YTickLabel',names_clusters(y_index'),'FontSize',10,'FontWeight','bold');
        else
            set(gca,'YTick',1:length(y_index),'YTickLabel',clust_idx(y_index'),'FontSize',10,'FontWeight','bold');
        end
    end;
end;
