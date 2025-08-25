**Network density**
1. Networks were derived per individual using https://github.com/cjl2007/PFM-Depression.
2. The Conte32k distance matrix was derived (available here) for all vertex-vertex distances.
3. Network density was derived by counting the unique # of network assignments within (6,8,10,12,14) mm of each vertex. Association network density was derived by counting the # of unique association networks in that radius (generate_network_density.m).
4. The difference between individual and group map density was plotted (network_density_plot.m, network_density_plot_yeo.m). 

**Rotation analysis of rostral CO region network density**
1. The rostral CO region was defined by identifying all contiguous CO clusters in the LPFC and then manually selecting the rostral one (make_rostral_pfc_mask.m).
2. Rostral CO regions were rotated 10k times to generate a null model to compare network density (generate_cifti_rotations.m).
3. The average association network density of the real location versus comparable spins was plotted (network_density_rotation_test.m).

**BOLD signal versus network density**
1. Average mode 1000 BOLD signal maps were derived for each subject and averaged across subjects (make_avg_snr.m)
2. The vertex-wise relationship between mode 1000 BOLD signal and density was plotted (snr_vs_density.m)
