% hypoclystering:
% 1) computes the reachability plot for the investigated first order
%    cluster

X = load('1_cluster_first_order.txt');  %load the dataset

maxEpsilon = 1.6;     %depends on your dataset; it must be larger than epsilon
minNumPoints = 200;

[order1, reachdist1]=clusterDBSCAN.discoverClusters(X,maxEpsilon,minNumPoints)

figure
bar(reachdist1(order1))
hold on 
yline(1)
axis([0 length(X(:,1)) 0 1.6])
xlabel('order index')
ylabel('reachability distance')
title('cluster C_1')