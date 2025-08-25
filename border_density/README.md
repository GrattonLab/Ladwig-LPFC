#Network Density 
1. Networks were derived per individual using https://github.com/cjl2007/PFM-Depression
2. The Conte32k distance matrix was derived (available here) for all vertex-vertex distances.
3. Network density was derived by counting the unique # of network assignments within X mm of each vertex. Association network density was derived by counting the # of unique association networks (generate_network_density.m).
4. The different between individual and group map density was plotted (network_density_plot.m, network_density_plot_yeo.m). 

#Rotation Analysis of rostral CO region network density
1. The rostral CO region was defined by identifying all contiguous CO clusters in the LPFC (make_rostral_pfc_mask.m) and then manually selecting the rostral one.
2. The rostral LPFC regions visualized were defined as all vertices within 15mm of the rostral CO region (for viz only). 
3. Rostral CO regions were inputs to the generate_cifti_rotations.m script which created 10k rotations of those regions.
4. The average association network density of the real location versus comparable spins was tested in network_density_rotation_test.m.

#BOLD signal versus Density 
1. Average mode 1000 BOLD signal maps were derived for each subject and averaged across subjects (make_avg_snr.m)
2. The vertex-wise relationship between mode 1000 BOLD signal and density was plotted (snr_vs_density.m)
