function [dip,strike,L,H,p,xp,yp,zp] = pca_planes(c11)
%This function checks for planar geometry condition and computes
%the planar surface parameters: dip and strike angle, length and height 

% check the sign of hypocenter depth
if c11(:,3)>0
    c11(:,3) = -c11(:,3);
end

% coordiantes of the cluster barycenter;
p = [mean(c11(:,1)), mean(c11(:,2)), mean(c11(:,3))]; 

% first step of PCA analysis: subtracting the mean
c_rs = [c11(:,1) c11(:,2) c11(:,3)];
c_rs = (c_rs - mean(c_rs)); 

% getting the covariance matrix 
covariance_rs = cov([c_rs(:,1),c_rs(:,2),c_rs(:,3)]);

% getting the eigenvalues and the eigenvectors 
[eigen_vector_rs, eigen_values_rs] = eig(covariance_rs);

eigen_value_1 = eigen_values_rs(1,1); %the smallest one
eigen_vector_1 =eigen_vector_rs(:,1);
eigen_value_2 = eigen_values_rs(2,2);
eigen_vector_2 =eigen_vector_rs(:,2);
eigen_value_3 = eigen_values_rs(3,3);
eigen_vector_3 =eigen_vector_rs(:,3);

%check the condition for planar geometry
cond_plane = ((2.0*eigen_value_1 < eigen_value_2) & ...
    (2.0*eigen_value_1 < eigen_value_3));

if cond_plane==1

 L = sqrt(12*eigen_value_3); %length of the plane
 H = sqrt(12*eigen_value_2); %height of the plane
 
% draw planes 
% get plane's normal vector
v3 = cross(eigen_vector_2,eigen_vector_3); 

% Points on plane
[ xp , yp ] = meshgrid( p(1)+((-L/2:L/2)) , p(2)+((-H/2:H/2)) );

% Equation for the plane
zp = p(3) - (v3(1)*(xp-p(1)) + v3(2)*(yp-p(2)))/v3(3);

% computation of the dip angle 
 norman=sqrt(eigen_vector_1(1,1)^2 + eigen_vector_1(2,1)^2 ...
 +eigen_vector_1(3,1)^2 );
 dip = acos(abs(eigen_vector_1(3,1))/norman)*180/pi;
 
% computation of strike angle 
 strike = atan(eigen_vector_1(1,1)/eigen_vector_1(2,1))*180/pi - 90;

 
else
    disp('condition for planar geometry is not satistied')
    dip=0; strike=0; L=0; H=0; xp=0; yp=0; zp=0;
end
