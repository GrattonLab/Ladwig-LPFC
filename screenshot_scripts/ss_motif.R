library(devtools)
library(ssbrain)
subjects=c("HUBS01", "HUBS02", "HUBS03", "HUBS04", "HUBS05", "HUBS06", "HUBS07", "HUBS08", "HUBS09", "HUBS10")
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
colors=c("RedWhiteBlue")


## INDIVIUDAL MOTIFS
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/networks_lpfc_motif/"
if (!dir.exists(output_folder)) dir.create(output_folder)
networks_view = c("DorsalAttention", "Language", "CinguloOpercular_Actionmode")
networks_name = paste(networks_view, collapse = "")
hemis=c("left")

for (s in 1:length(subjects)) {
  subject=subjects[s]
  sub_surfL=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/preproc_fmriprep-23.2.0a2_all/freesurfer/FREESURFER_fs_LR/sub-", subject, "/MNI/fsaverage_LR32k/sub-", subject, ".L.very_inflated.32k_fs_LR.surf.gii")
  sub_surfR=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/preproc_fmriprep-23.2.0a2_all/freesurfer/FREESURFER_fs_LR/sub-", subject, "/MNI/fsaverage_LR32k/sub-", subject, ".R.very_inflated.32k_fs_LR.surf.gii")
  subject_brain = ss_surf(surfL = sub_surfL, surfR = sub_surfR)
  networks=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/lynch/alt/Bipartite_PhysicalCommunities+FinalLabeling_lpfc_alt.dlabel.nii");
  for (h in 1:length(hemis)) {
    hemi=hemis[h]
    border_brain = group_brain + 
      ss_dlabel(filename = networks, labels = networks_view) +
      ss_view(side='lateral', rotation = c(30,60,90))
    output_file=paste0(output_folder,subject,'_', networks_name, '_',hemi,'.png')
    captureBrain(border_brain, hemisphere = hemi, filename = output_file)
  }
}

# GROUP MOTIF
group_networks=paste0("/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc_alt.dlabel.nii")
networks_view = c("DorsalAttention", "Language", "CinguloOpercular_Actionmode")
networks_name = paste(networks_view, collapse = "")
hemis=c("left")

for (h in 1:length(hemis)) {
  hemi=hemis[h]
  for (v in 1:length(views)) {
    view = views[v]
    border_brain = group_brain + 
      ss_dlabel(filename = group_networks, labels = networks_view) +
      ss_view(side=view, rotation = c(30,60,90))
    output_file=paste0(output_folder, 'group_', networks_name, '_',hemi,'_',view,'.png')
    captureBrain(border_brain, hemisphere = hemi, filename = output_file)
  }
}
# INDIVIDUAL SEEDS WITH CO REGION PRESENT
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/seeds_lpfc_motif/"
if (!dir.exists(output_folder)) dir.create(output_folder)
hemi="left"
networks = c("DorsalAttention", "Language")
for (s in 1:length(subjects)) {
  subject=subjects[s]
  label_networks=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/functional_masks/alt/rostral_co_lpfc_alt.dlabel.nii");
  
  print(subject)
  dconn = importCifti(paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/sub-", subject, "_allsess_smooth1_tmasked_alt.dconn.nii"),data_only=TRUE)
  seed_list = read.csv(paste0("/projects/b1081/NSF_HUBS/resources/", subject, "_seeds.txt"), header = FALSE, stringsAsFactors = FALSE)
  
  networks_borders_l=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/lynch/alt/Bipartite_PhysicalCommunities+FinalLabeling_lpfc_L.border")
  networks_borders_r=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/lynch/alt/Bipartite_PhysicalCommunities+FinalLabeling_lpfc_R.border")
  networks_borders = c(networks_borders_l,networks_borders_r)  
  for (n in 1:length(networks)) {
    network = networks[n]
    print(network)
    seed = seed_list$V1[seed_list$V2 == network]
    
    seed_brain = group_brain + 
      ss_dlabel(filename=label_networks,labels="Rostral_LPFC_CO") + 
      ss_seed(seed_value = seed, colorrange = c(0, 0.5), seed_sphere_size=3,seed_data = dconn) + 
      ss_border(networks_borders_l,hemisphere = "left", borders = network, colors=c("black")) + 
      ss_view(side='lateral', rotation = c(30,60,90))
    output_file=paste0(output_folder,subject,'_',seed,'_',network,'_',hemi,'.png')
    captureBrain(seed_brain, hemisphere = hemi, filename = output_file)
  }
  
  rm(dconn)
}

# GROUP SEEDS FOR MOTIFS  WITH CO PRESENT
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/networks_motifs_seeds_co_rostral/"
if (!dir.exists(output_folder)) dir.create(output_folder)
hemi="left"
networks = c("DorsalAttention", "Language")
subject ="group"

label_networks=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/functional_masks/alt/rostral_co_lpfc_alt.dlabel.nii");
dconn = importCifti(paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/sub-", subject, "_allsess_smooth1_tmasked_alt.dconn.nii"),data_only=TRUE)
seed_list = read.csv(paste0("/projects/b1081/NSF_HUBS/resources/", subject, "_seeds.txt"), header = FALSE, stringsAsFactors = FALSE)
group_networks_borders_l=paste0("/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc_L.border")
group_networks_borders_r=paste0("/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc_R.border")
group_networks_borders = c(group_networks_borders_l,group_networks_borders_r)

for (n in 1:length(networks)) {
  network = networks[n]
  print(network)
  seed = seed_list$V1[seed_list$V2 == network]
  
  seed_brain = group_brain + 
    ss_dlabel(filename=label_networks,labels="Rostral_LPFC_CO") + 
    ss_seed(seed_value = seed, thresh = 0.075, colorrange = c(0, 0.25), seed_sphere_size=3,seed_data = dconn) + 
    ss_border(group_networks_borders_l,hemisphere = "left", borders = network, colors=c("black")) + 
    ss_view(side='lateral', rotation = c(30,60,90))
  output_file=paste0(output_folder,subject,'_',seed,'_',network,'_',hemi,'.png')
  captureBrain(seed_brain, hemisphere = hemi, filename = output_file)
}


## Individual Motif TPJ
output_folder="/projects/b1081/NSF_HUBS/images/networks_tpj_motif/"
if (!dir.exists(output_folder)) dir.create(output_folder)
networks_view = c("Language", "DorsalAttention", "CinguloOpercular_Actionmode")
networks_name = paste(networks_view, collapse = "")
hemis =c("left")
#for (s in 1:length(subjects)) {
subject=subjects[s]
sub_surfL=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/preproc_fmriprep-23.2.0a2_all/freesurfer/FREESURFER_fs_LR/sub-", subject, "/MNI/fsaverage_LR32k/sub-", subject, ".L.very_inflated.32k_fs_LR.surf.gii")
sub_surfR=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/preproc_fmriprep-23.2.0a2_all/freesurfer/FREESURFER_fs_LR/sub-", subject, "/MNI/fsaverage_LR32k/sub-", subject, ".R.very_inflated.32k_fs_LR.surf.gii")
subject_brain = ss_surf(surfL = sub_surfL, surfR = sub_surfR)
networks=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/Bipartite_PhysicalCommunities+FinalLabeling_alt.dlabel.nii");
for (h in 1:length(hemis)) {
  hemi=hemis[h]
  border_brain = group_brain + 
    ss_dlabel(filename = networks, labels = networks_view) +
    #ss_dlabel(filename = networks) +
    ss_view(side='lateral', rotation = c(340,70,80))
  output_file=paste0(output_folder,subject,'_', networks_name, '_',hemi,'.png')
  captureBrain(border_brain, hemisphere = hemi, filename = output_file)
}
}
group_networks=paste0("/projects/b1081/NSF_HUBS/resources/lynch_group_prior50_alt.dlabel.nii")

for (h in 1:length(hemis)) {
  hemi=hemis[h]
  for (v in 1:length(views)) {
    view = views[v]
    border_brain = group_brain + 
      ss_dlabel(filename = group_networks, labels = networks_view) +
      ss_view(side='lateral', rotation = c(340,70,80))
    output_file=paste0(output_folder, 'group_', networks_name, '_',hemi,'_',view,'.png')
    captureBrain(border_brain, hemisphere = hemi, filename = output_file)
  }
}


## Individual Motif pre-SMA
output_folder="/projects/b1081/NSF_HUBS/images/networks_sma_motif/"
if (!dir.exists(output_folder)) dir.create(output_folder)
networks_view = c("Language", "DorsalAttention", "CinguloOpercular_Actionmode")
networks_name = paste(networks_view, collapse = "")
hemis =c("left")
#for (s in 1:length(subjects)) {
subject=subjects[s]
sub_surfL=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/preproc_fmriprep-23.2.0a2_all/freesurfer/FREESURFER_fs_LR/sub-", subject, "/MNI/fsaverage_LR32k/sub-", subject, ".L.very_inflated.32k_fs_LR.surf.gii")
sub_surfR=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/preproc_fmriprep-23.2.0a2_all/freesurfer/FREESURFER_fs_LR/sub-", subject, "/MNI/fsaverage_LR32k/sub-", subject, ".R.very_inflated.32k_fs_LR.surf.gii")
subject_brain = ss_surf(surfL = sub_surfL, surfR = sub_surfR)
networks=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/Bipartite_PhysicalCommunities+FinalLabeling_alt.dlabel.nii");
for (h in 1:length(hemis)) {
  hemi=hemis[h]
  border_brain = group_brain + 
    ss_dlabel(filename = networks, labels = networks_view) +
    #ss_dlabel(filename = networks) +
    ss_view(side='lateral', rotation = c(15,-30,-70))
  output_file=paste0(output_folder,subject,'_', networks_name, '_',hemi,'.png')
  captureBrain(border_brain, hemisphere = hemi, filename = output_file)
}
}
group_networks=paste0("/projects/b1081/NSF_HUBS/resources/lynch_group_prior50_alt.dlabel.nii")

for (h in 1:length(hemis)) {
  hemi=hemis[h]
  for (v in 1:length(views)) {
    view = views[v]
    border_brain = group_brain + 
      ss_dlabel(filename = group_networks, labels = networks_view) +
      ss_view(side='lateral', rotation = c(15,-30,-70))
    output_file=paste0(output_folder, 'group_', networks_name, '_',hemi,'_',view,'.png')
    captureBrain(border_brain, hemisphere = hemi, filename = output_file)
  }
}
