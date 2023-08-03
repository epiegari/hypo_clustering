
% DBSCAN_soc: 
% 1) runs DBSCAN using epsilon value found by OPTICS analysis and
%    uses second order clusters as input data
% 2) saves third order clusters in different files;
% 3) calls the function pca_planes and checks if a second-order cluster
%    can be associated to a planar surface;
% 4) makes a 3D plot with cluster data and the associated plane, 
%    for each cluster that satisfies the condition for planar geometry;
% 5) saves in the matrix At the planar surface parameters (dip angle, 
%    strike angle, length, height) and the coordinates of the baricenter,
%    which are computed by the function pca_planes.


%% Read data 
% The data file contains the UTM coordinates of hypocenters in the following 
% order: Easting, Nothing and depth. They are measured in km.

data=load('2_cluster_second_order.txt');
x=data(:,1);  % UTM East (km)
y=data(:,2);  % UTM North (km)
depth=data(:,3); % hypocenter depth (km) - negative values
N = length(data); % total number of hypocenters 


%% Run DBSCAN Clustering Algorithm

X=[x y depth]; %non-scaled data
 
% set the input parameters
epsilon=input('enter the initial value for epsilon (km) '); %neighbourhood radius (km)
Z=input('enter the initial value for Z '); %minimum # of points in the eps-neighbouhood

idx=dbscan(X,epsilon,Z);
[GC,GR]=groupcounts(idx);
g = [GC GR];
      
%% save second order clusters in different files
    
%%%%%%%%%%% cluster separation %%%%%%%%%%%%%%
utm_E_c=X(:,1);
utm_N_c=X(:,2);
dep=X(:,3);
g(find(GR==-1),:)=[]; 
s=sortrows(g,'descend');

ntoc_pca=input('enter the number of third order clusters for PCA analysis ')

figure

for i=1:ntoc_pca
savecto = zeros(length(idx),3);
for n=1:length(idx)
    if (idx(n) == s(i,2))
        dataf=[utm_E_c(n), utm_N_c(n), dep(n)];
    savecto(n,:)=dataf;
    end
end
new_savecto=savecto(any(savecto,2),:);
save(sprintf('%d_%s',i,'cluster_third_order.txt'),'new_savecto','-ascii');
filename = [num2str(i),'_cluster_third_order.txt'];
cto=load(filename);

As(ntoc_pca,7)=0; %inizialize and set dimensions of A

[dip,strike,L,H,p,xp,yp,zp] = pca_planes(cto); %call the function pca_planes
At(i,:)=[dip strike L H p]   %saves in At the fault parameters and the
                            %coordinates of the cluster barycenter
% if xp~=0
if dip > 0
    surf(xp,yp,zp)  % plot the plane
    shading interp
    alpha(0.3)
    colormap winter
end
%%plot the first order cluster
hold on
grid on
scatter3(cto(:,1),cto(:,2),cto(:,3),1,'filled') 
xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('depth (km)')
axis equal

end


    