library(devtools)
library(ssbrain)
subjects=c("HUBS01", "HUBS02", "HUBS03", "HUBS04", "HUBS05", "HUBS06", "HUBS07", "HUBS08", "HUBS09", "HUBS10")
#set_wbpath('/projects/b1081/Scripts/workbench2/bin_linux64/wb_command')
set_wbpath('singularity exec --env PATH=/opt/workbench/bin/:/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin,WORKBENCH_FFMPEG_DIR=/usr/bin/ -B /projects:/projects /software/singularity/images/connectome-workbench-v1.5.0.sif /opt/workbench/bin/wb_command $@')
hemis=c("left","right")
views=c("lateral","medial")
lh_idx = array(c(rep(1, 32492), rep(0, 32492)))
rh_idx = array(c(rep(0, 32492), rep(1, 32492)))
group_surfR=paste0("/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.R.very_inflated.32k_fs_LR.surf.gii")
group_surfL=paste0("/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.L.very_inflated.32k_fs_LR.surf.gii")
group_brain = ss_surf(surfL = group_surfL, surfR = group_surfR)
lpfc_mask_file = "/projects/b1081/NSF_HUBS/resources/group_lpfc_mask_alt.dscalar.nii"
lpfc_mask=importCifti(lpfc_mask_file,data_only=TRUE)
lpfc_mask = drop(lpfc_mask)
#colors=c("FSL")
colors=c("RedWhiteBlue")

## Rostral High Density Zone Networks
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/lpfc_rostral_hdz/"
hemis = c("left")
if (!dir.exists(output_folder)) dir.create(output_folder)
for (s in 1:length(subjects)) {
  subject=subjects[s]
  networks=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/lynch/alt/Bipartite_PhysicalCommunities+FinalLabeling_rostral_lpfc_alt.dlabel.nii");
  for (h in 1:length(hemis)) {
    hemi=hemis[h]
    border_brain = group_brain + 
      ss_dlabel(filename = networks) +
      ss_view(side='lateral', rotation = c(30,60,90))
    output_file=paste0(output_folder,subject,'_',hemi,'.png')
    captureBrain(border_brain, hemisphere = hemi, filename = output_file)
  }
}

# DN-A seeds for exception subjects
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/articulation_point_seeds/"
if (!dir.exists(output_folder)) dir.create(output_folder)
networks = c("Default_Retrosplenial")
subjects = c("HUBS07")
for (s in 1:length(subjects)) {
  subject=subjects[s]
  border_file_L = paste0('/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-', subject, '/rostral_lpfc_mask_alt_L.border')
  border_file_R = paste0('/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-', subject, '/rostral_lpfc_mask_alt_R.border')
  borders = c(border_file_L,border_file_R)  
  
  networks_borders_l=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/Bipartite_PhysicalCommunities+FinalLabeling_L.border")
  networks_borders_r=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/Bipartite_PhysicalCommunities+FinalLabeling_R.border")
  networks_borders = c(networks_borders_l,networks_borders_r)  
  
  print(subject)
  dconn = importCifti(paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/sub-", subject, "_allsess_smooth1_tmasked_alt.dconn.nii"),data_only=TRUE)
  seed_list = read.csv(paste0("/projects/b1081/NSF_HUBS/resources/", subject, "_seeds_articulation.txt"), header = FALSE, stringsAsFactors = FALSE)
  for (n in 1:length(networks)) {
    network = networks[n]
    print(network)
    seed = seed_list$V1[seed_list$V2 == network]
    for (h in 1:length(hemis)) {
      hemi = hemis[h]
      seed_brain = group_brain + 
        ss_seed(seed_value = seed, colorrange = c(0, 0.5), seed_sphere_size=3,seed_data = dconn) +
        ss_border(filename=networks_borders[h],borders = network, hemisphere = hemi, width=4) +
        ss_border(filename=borders[h],hemisphere = hemi, width=4, colors = c("black"))
      for (v in 1:length(views)) {
        view = views[v]
        print(view)
        view_brain = seed_brain +
          ss_view(side=view, rotation = c(30,60,90))
        output_file=paste0(output_folder,subject,'_',seed,'_',network,'_',hemi,view,'.png')
        captureBrain(view_brain, hemisphere = hemi, filename = output_file)
      }
    }
  }
  rm(dconn)
}