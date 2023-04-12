
clear;
clc;
%% Example
% data
data1 =  0.9*randn(2000,4000);
data2 = 0.3*randn(2000,16000);
data3 = 0.3*randn(3000,4000);
data4 = 0.9*randn(3000,6000);
data5 = 0.3*randn(3000,5000);
data6 = 0.9*randn(3000,5000);
data7 = 0.3*randn(5000,4000);
data8 = 0.3*randn(5000,6000);
data9 = 0.9*randn(5000,10000);
% row labels
rlabel1 = ones(2000,1);
rlabel2 = 2*ones(3000,1);
rlabel3 = 3*ones(5000,1);
% column labels
clabel1 = ones(4000,1);
clabel2 = 2*ones(6000,1);
clabel3 = 3*ones(5000,1);
clabel4 = 4*ones(5000,1);


rlabels = [rlabel1;rlabel2;rlabel3];
clabels = [clabel1;clabel2;clabel3;clabel4];
data = [data1,data2;data3,data4,data5,data6;data7,data8,data9];

%???????Fast Flexible Bipartite Graph Model for Co-Clustering?????
 stop_iter = 10;
 sparse_ratio = 0.7;
 pruning = 0.0001; 
 c1 = length(unique(rlabels));
 c2 = length(unique(clabels));
 
 tic
 [r,c] = FBGPC(data, c1, c2, stop_iter, sparse_ratio, pruning);
 run_time = toc
 rAccuracy = calculateAccuracy(rlabels',r')
 cAccuracy = calculateAccuracy(clabels',c')

 