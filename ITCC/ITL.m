function [row_cluster,column_cluster] = ITL(matrix,RowClusterNo,ColClusterNo)
% Function to find biclusters using Information Theoretic learning
% The algorithm was developed for simultaneously clustering the rows and
% columns of contingency table. It views contingency tables as joint
% probability distribution of the two discrete random variables. Hence the
% problem of simultaneously clustering reduces to maximizing the mutual
% information between the clustered random variables or minimizing the loss
% in mutual information as a result of clustering.
% The algorithm follows the paper by:
% Dhillon, Inderjit S., Subramanyam Mallela, and Dharmendra S. Modha.
% "Information-theoretic co-clustering." In Proceedings of the ninth ACM SIGKDD 
% international conference on Knowledge discovery and data mining, pp. 89-98. ACM, 2003.
% 
% Inputs:
%        matrix          =   Input matrix which represents the joint probability
%                            distribution of the two variables
%        RowClusterNo    =   Desired number of row clusters
%        ColClustNo      =   Desired number of column clusters
%
% Output:
%        row_cluster     =   Vector containing the cluster number of each row.
%        column_cluster  =   Vector containing the cluster number of each column.   
%
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, Indian Institute of Technology, Kanpur, India

%% Input argument check
if nargin < 1
    error('No Input argument specified.');
end
if nargin < 2
    RowClusterNo = 10;
    ColClusterNo = 10;
end
if nargin < 3
    ColClusterNo = RowClusterNo;
end
%% Calling the main program
[row_cluster,column_cluster] = InformationTheoreticLearning(matrix,RowClusterNo,ColClusterNo);
end
