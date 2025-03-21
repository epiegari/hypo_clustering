# hypo_clustering
The scripts in this repository facilitate imaging of segmented fault surfaces illuminated by the hypocenter clustering method presented in Piegari et al. (2024), published in Earth and Space Science and available at https://doi.org/10.1029/2023EA003267. The script related to using Principal Component Analysis (PCA) for plane representation has been updated with a valuable contribution from Cliff Thurber, to whom I extend my sincere thanks.

## Scientific publication
If you use the codes in this repository, please cite the publication: https://doi.org/10.1029/2023EA003267

## Installation
To make these MATLAB scripts work on your machine you need to install the Statistics and Machine Learning Toolbox of MATLAB. 

## Modules
The MATLAB codes stored in the folder are:
- DBSCAN + PCA for identifying first order fault surfaces: dbscan_cfo.m
- DBSCAN + PCA for identifying second order fault surfaces: dbscan_cso.m
- DBSCAN + PCA for identifying third order fault surfaces: dbscan_cto.m
- OPTICS for determination of reachability plot: optics_hypoclustering.m
- PCA analysis and fault parameter estimation: pca_planes.m
- PCA analysis and fault parameter estimation: pca_planes_updated.m
- hypocenter density computation and visualization: hypodens.m

## Usage
The file dbscan_cfo.m searchs for a cluster solution in the Crossover Region (as presented in Piegari et al., https://doi.org/10.1093/gji/ggac160), while this condition is not required for the identification of second and third order clusters. Since the choice of the input parameters is critical, the file hypodens.m can be run to visualize the hypocenter density and get a tip on the Z value to select. 

Please note:
- The hypocenter input file needs UTM coordinates of hypocenters
