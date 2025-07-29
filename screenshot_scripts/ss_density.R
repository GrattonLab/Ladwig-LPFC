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


# INDIVIDUAL ASSOCIATION_NETWORK_DENSITY
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/association_network_density/"
if (!dir.exists(output_folder)) dir.create(output_folder)
hemis=c("left")
views = c("lateral")
#subjects = c("group")
for (s in 1:length(subjects)) {
  subject=subjects[s]
  task_file = paste0('/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-', subject, '/density/lynch/alt/', subject, '_association_network_density_avg_alt.dscalar.nii')
  task_data=importCifti(task_file,data_only=TRUE)
  pos_data = task_data[task_data > 0]
  filtered_data = pos_data[!is.nan(pos_data)]
  hi=quantile(filtered_data, 0.99, na.rm=TRUE)
  thresh_file=file(paste0(output_folder,subject,'_association_colorbar.txt'))
  writeLines(as.character(c(hi, -hi)),thresh_file)
  task_brain = group_brain +
    ss_dscalar(filename=task_file,colorbar='magma',pos_colorrange=c(0,hi), neg_colorrange = c(-hi,0)) 
  for (h in 1:length(hemis)) {
    hemi=hemis[h]
    border_brain = task_brain 
    for (v in 1:length(views)) {
      view = views[v]
      view_brain = border_brain + 
        ss_view(side=view,rotation=c(30,60,90))
      output_file=paste0(output_folder,subject,'_association_network_density_' ,hemi,'_',view,'.png')
      captureBrain(view_brain, hemisphere = hemi, filename = output_file)   
    }
  }
}

# AVG GROUP ASSOCIATION_NETWORK_DENSITY

output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/association_network_density/"
if (!dir.exists(output_folder)) dir.create(output_folder)
hemis=c("left")
views = c("lateral")
types = c("group","avg","avg-group")
subject="group"
for (t in 1:length(types)) {
  task_file = paste0('/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/lynch/alt/', types[t], '_association_network_density_avg_alt.dscalar.nii')
  task_data=importCifti(task_file,data_only=TRUE)
  pos_data = task_data[task_data > 0]
  filtered_data = pos_data[!is.nan(pos_data)]
  if(types[t] == "avg-group") {
    hi=quantile(filtered_data, 0.95, na.rm=TRUE)
  } else {
    hi=quantile(filtered_data, 0.99, na.rm=TRUE)
  }
  
  thresh_file=file(paste0(output_folder,types[t],'_association_colorbar.txt'))
  writeLines(as.character(c(hi, -hi)),thresh_file)
  task_brain = group_brain +
    ss_dscalar(filename=task_file,colorbar='magma',pos_colorrange=c(0,hi), neg_colorrange = c(-hi,0)) 
  for (h in 1:length(hemis)) {
    hemi=hemis[h]
    border_brain = task_brain 
    for (v in 1:length(views)) {
      view = views[v]
      view_brain = border_brain + 
        ss_view(side=view,rotation=c(30,60,90))
      output_file=paste0(output_folder,types[t],'_association_network_density_' ,hemi,'_',view,'.png')
      captureBrain(view_brain, hemisphere = hemi, filename = output_file)   
    }
  }
}


# YEO AVG GROUP ASSOCIATION_NETWORK_DENSITY

output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/association_network_density_yeo/"
if (!dir.exists(output_folder)) dir.create(output_folder)
hemis=c("left")
views = c("lateral")
types = c("group","avg","avg-group")
subject="group"
for (t in 1:length(types)) {
  task_file = paste0('/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/yeo/alt/', types[t], '_association_network_density_avg_alt.dscalar.nii')
  task_data=importCifti(task_file,data_only=TRUE)
  pos_data = task_data[task_data > 0]
  filtered_data = pos_data[!is.nan(pos_data)]
  if(types[t] == "avg-group") {
    hi=quantile(filtered_data, 0.95, na.rm=TRUE)
  } else {
    hi=quantile(filtered_data, 0.99, na.rm=TRUE)
  }
  thresh_file=file(paste0(output_folder,types[t],'_',thresh,'_association_colorbar.txt'))
  writeLines(as.character(c(hi, -hi)),thresh_file)
  task_brain = group_brain +
    ss_dscalar(filename=task_file,colorbar='magma',pos_colorrange=c(0,hi), neg_colorrange = c(-hi,0)) 
  for (h in 1:length(hemis)) {
    hemi=hemis[h]
    border_brain = task_brain 
    for (v in 1:length(views)) {
      view = views[v]
      view_brain = border_brain + 
        ss_view(side=view,rotation=c(30,60,90))
      output_file=paste0(output_folder,types[t],'_association_network_density_' ,hemi,'_',view,'.png')
      captureBrain(view_brain, hemisphere = hemi, filename = output_file)   
    }
  }
}