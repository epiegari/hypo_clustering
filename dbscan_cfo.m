
% DBSCAN_cfo: 
% 1) checks if a data scaling with min-max normalization is required;
% 2) searchs for a cluster solution in the Crossover Region(CR)and makes a 
%    3D plot of such a solution that identifies first order clusters; 
% 3) saves first order clusters in different files;
% 4) calls the function pca_planes and checks if a first-order cluster
%    can be associated to a planar surface;
% 5) makes a 3D plot with cluster data and the associated plane, 
%    for each cluster that satisfies the condition for planar geometry;
% 6) saves in the matrix A the planar surface parameters (dip angle, 
%    strike angle, length, height) and the coordinates of the baricenter,
%    which are computed by the function pca_planes.

%Note: in the search for a CR solution, to reduce the computation time a
% maximum number of iterations is introduced; if such a threshold is reached
% unsuccessfully, a change in the input parameter epsilon is suggested
%(doubling it could be helpful). 

%Note: to have an indication about the initial value of Z, run hypodens.m
%that computes the maximum density of hypocenters and makes plots to
%visualize the regions with the highest density values. We initially set Z 
% to about half the maximum density. This choice is different from that
%suggested in Piegari et al., Geophys. J. Int. (2022) 230, 2073–2088.
% The reason is that here we are interested in finding a solution 
%in the lower part of the CR, where the condition to find core points is
%not so strict, (Z/epsilon is not so large). This allows us to find
%regions with the highest values of hypocenter density as nested
%structures included in first order clusters, and then illuminate them  
%as second (or third) order clusters.


%% Read data 
% The data file contains the UTM coordinates of hypocenters in the following 
% order: Easting, Nothing and depth. They are measured in km.

data=load('Taiwandataset_UTM.txt');
x=data(:,1);  % UTM East (km)
y=data(:,2);  % UTM North (km)
depth=-data(:,3); % hypocenter depth (km) - insert positive values
N = length(data); % total number of hypocenters 


%% 1) checks if a data scaling with min-max normalization is required;
horext_x = max(x)-min(x); % horizontal extension (Easting)
horext_y = max(y)-min(y); % horizontal extension  (Northing)
vertext = max(depth)-min(depth); % vertical extension

cond_iso = (vertext < 0.5*horext_x) | (vertext < 0.5*horext_y);

if (cond_iso==1) 
 % scale the data in the depth range using the min-max transformation
 % set the initial DBSCAN input parameters

 disp 'data scaling is required'
 % set the depth range
 dmax=input('enter an initial value for dmax (km) '); %maximum boundary
 dmin=input('enter an initial value for dmin (km) '); %minimum boundary
 
 OldMax_E=max(x);
 OldMin_E=min(x);
 NewMax_E=dmax;
 NewMin_E=dmin;
 OldRange_E = (OldMax_E - OldMin_E);  
 NewRange_E = (NewMax_E - NewMin_E);  
 NewValue_E = (((x - OldMin_E) * NewRange_E) / OldRange_E) + NewMin_E;

 OldMax_N=max(y);
 OldMin_N=min(y);
 NewMax_N=dmax;
 NewMin_N=dmin;
 OldRange_N = (OldMax_N - OldMin_N);  
 NewRange_N = (NewMax_N - NewMin_N);  
 NewValue_N = (((y - OldMin_N) * NewRange_N) / OldRange_N) + NewMin_N;

 X=[NewValue_E(:,1) NewValue_N(:,1) -depth]; %scaled data
 
else
 
 disp 'data scaling is not required'
    
 X=[x y -depth]; %non-scaled data
 
end

% set the input parameters
epsilon=input('enter the initial value for epsilon (km) '); %neighbourhood radius (km)
Z=input('enter the initial value for Z '); %minimum # of points in the eps-neighbouhood

%% Run DBSCAN Clustering Algorithm

idx=dbscan(X,epsilon,Z);
[GC,GR]=groupcounts(idx);
g = [GC GR];


%% 2) identifies a cluster solution in the Crossover Region (CR)

% CR is defined by the two conditions: 
% # of noise points < 60% of data and # of points belonging to the
% biggest cluster < 60% of data

N_noise = GC(find(GR==-1))  %number of noise points
N_Cb = max(GC)     %number of points belonging to the biggest cluster
cond = (N_noise<0.6*N) & (N_Cb<0.6*N);  %condition for belonging to CR

if (cond==1)
    disp 'a solution in CR is found at first iteration'
else
    disp 'a change in input parameters is required'
    nmax=11; %maximum number of iterations for CR solution search
    niter=0; %number of iteration
    while(cond==0 & niter<nmax)
        Zstep=50;   % step increment to change the Z value
        if (N_Cb>0.6*N)
            Z = Z+Zstep;
            disp 'approaching CR from the right'
        end
        if (N_noise>0.6*N)
            Z = Z-Zstep;
            disp 'approaching CR from above'
        end
        idx=dbscan(X,epsilon,Z);
        [GC,GR]=groupcounts(idx);
        g = [GC GR];
        N_noise = GC(find(GR==-1));  
        N_Cb = max(GC);     
        cond = (N_noise<0.6*N) & (N_Cb<0.6*N);  
        niter = niter + 1;
    end
end 
    
    
% makes a plot of the solution that identifies first order clusters
figure
if (cond_iso==1)
    %back-mapping to the real space
    OldMax_xr=dmax;
    OldMin_xr=dmin;
    NewMax_xr=max(x);
    NewMin_xr=min(x);
    OldRange_xr = (OldMax_xr - OldMin_xr);  
    NewRange_xr = (NewMax_xr - NewMin_xr);  
    NewValue_xr = (((X(:,1) - OldMin_xr) * NewRange_xr) / OldRange_xr) + NewMin_xr;

    OldMax_yr=dmax;
    OldMin_yr=dmin;
    NewMax_yr=max(y);
    NewMin_yr=min(y);
    OldRange_yr = (OldMax_yr - OldMin_yr);  
    NewRange_yr = (NewMax_yr - NewMin_yr);  
    NewValue_yr = (((X(:,2)- OldMin_yr) * NewRange_yr) / OldRange_yr) + NewMin_yr;

    scatter3(NewValue_xr,NewValue_yr,X(:,3),1,idx,'filled'); 
    X(:,1)=NewValue_xr; X(:,2)=NewValue_yr;
else
 
    scatter3(X(:,1),X(:,2),X(:,3),2,idx,'filled');
end
    grid on
    xlabel('Easting (km)')
    ylabel('Northing (km)')
    zlabel('depth (km)')
    title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ',',...
    ' Z = ' num2str(Z) ',', ...
    ' Ncluster = ' num2str(max(idx)) ')'] );
    axis equal
    
%% save first order clusters in different files
    
%%%%%%%%%%% cluster separation %%%%%%%%%%%%%%
utm_E_c=X(:,1);
utm_N_c=X(:,2);
dep=X(:,3);
g(find(GR==-1),:)=[]; 
s=sortrows(g,'descend');

nfoc_pca=input('enter the number of first order clusters for PCA analysis ')
figure

for i=1:nfoc_pca
savecfo = zeros(length(idx),3);
for n=1:length(idx)
    if (idx(n) == s(i,2))
        dataf=[utm_E_c(n), utm_N_c(n), dep(n)];
    savecfo(n,:)=dataf;
    end
end
new_savecfo=savecfo(any(savecfo,2),:);
save(sprintf('%d_%s',i,'cluster_first_order.txt'),'new_savecfo','-ascii');
filename = [num2str(i),'_cluster_first_order.txt'];
cfo=load(filename);

A(nfoc_pca,7)=0; %inizialize and set dimensions of A

[dip,strike,L,H,p,xp,yp,zp] = pca_planes(cfo); %call the function pca_planes
A(i,:)=[dip strike L H p]   %saves in A the fault parameters and the
                            %coordinates of the cluster barycenter
if xp~=0
    surf(xp,yp,zp)  % plot the plane
    shading interp
    alpha(0.3)
    colormap winter
end
%%plot the first order cluster
hold on
grid on
scatter3(cfo(:,1),cfo(:,2),cfo(:,3),1,'filled') 
axis equal

end


    