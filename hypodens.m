
%hypodens:
% 1) compute the density of hypocenters, i.e. the number of hypocenters in  
%    a sphere of radius R 
% 2) makes a scatter plot of hypocenter locations with colors representing
%    density of hypocenters
% 3) asks for a threshold on density and makes a scatter plot of hypocenter
%    locations with colors representing density of hypocenters above the
%    threshold


%% Read data
data=load('Aquiladataset_UTM.txt');
xE=data(:,7);  % UTM East (km)
yN=data(:,8);  % UTM North (km)
d=data(:,9); % hypocenter depth (km)


%% hypocenter density computation
xyzp1=[xE yN d];

ptCloud1 = pointCloud(xyzp1);

x1 = ptCloud1.Location(:,1);
y1 = ptCloud1.Location(:,2);
z1 = ptCloud1.Location(:,3);

R = 1; % depends on your data (in km)
for idx1=1:length(xyzp1(:,1))
Distances1 = sqrt((x1-xyzp1(idx1,1)).^2 + (y1-xyzp1(idx1,2)).^2 + ...
                 (z1-xyzp1(idx1,3)).^2);
Ninside1   = length( find(Distances1<=R) ); 
density1(idx1) = Ninside1/(4*pi*R.^3/3); 
end

c1d = [xyzp1 density1']; 
densmax=max(density1); %maximum value of density hypoceneter
c1_E=c1d(:,1); c1_N=c1d(:,2); c1_dep=c1d(:,3); dd1=c1d(:,4);


%% make a scatter plot with colors representing density values

figure
scatter3(c1_E,c1_N,c1_dep,3,dd1','filled')

xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('depth (km)')
axis equal
colorbar

%% make the same plot as previous one but with a threshold on density

thd=input('enter a threshold value for a better visualization ') 
 
%identify only hypoceneters for which the density value is above the
%threshold
c1_Ef=c1_E(find(dd1>thd)); c1_Nf=c1_N(find(dd1>thd)); 
c1_depf=c1_dep(find(dd1>thd)); dd1f=dd1(find(dd1>thd));

figure
scatter3(c1_Ef,c1_Nf,c1_depf,3,dd1f','filled')
xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('depth (km)')
axis equal
colorbar