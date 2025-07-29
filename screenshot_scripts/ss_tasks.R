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


# Task activations X individual-specific networks
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/lpfc_networks_tasks/"
if (!dir.exists(output_folder)) dir.create(output_folder)
#tasks=c("viswm","audwm","msit","vmsit","audattn","visattn","spatialwm","verbalwm")
tasks=c("tom")
#networks_view_borders=c("DorsalAttention","CinguloOpercular_Actionmode","Frontoparietal")
networks_view_borders=c("Default_Parietal")
control_net_colors = c("#0AF014","#5A28A0","#AD9100")
hemis=c("left")
for (s in 1:length(subjects)) {
  subject=subjects[s]
  networks_borders_l=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/lynch/Bipartite_PhysicalCommunities+FinalLabeling_lpfc_L.border")
  networks_borders_r=paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-", subject, "/pfm/lynch/Bipartite_PhysicalCommunities+FinalLabeling_lpfc_R.border")
  networks_borders = c(networks_borders_l,networks_borders_r)
  for (t in 1:length(tasks)) {
    task=tasks[t]
    #task_file = paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-", subject, "/task_summaries/alt/sub-", subject, "_", task, "_zstats_mean_alt.dscalar.nii")
    task_file = paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-", subject, "/domain_summaries/alt/sub-", subject, "_", task, "_zstats_mean_alt.dscalar.nii")
    
    task_data=importCifti(task_file,data_only=TRUE)
    task_data[lpfc_mask==0] = 0
    lpfc_task_data= task_data[lpfc_mask == 1 & rh_idx ==1]
    print(mean(lpfc_task_data <= 1.5) * 100)
    print(quantile(lpfc_task_data,0.75))
    pos_lpfc_task_data=lpfc_task_data[lpfc_task_data > 0]
    
    hi=quantile(lpfc_task_data, 0.98)
    mid_hi=max(quantile(lpfc_task_data, 0.75),0)
    #lo = quantile(lpfc_task_data, 0.02)
    thresh_file=file(paste0(output_folder,subject,'_',task,'.txt'))
    writeLines(as.character(c(hi,mid_hi)),thresh_file)
    task_brain = group_brain +
      ss_dscalar(colorbar=colors, show="pos", pos_thresh = max(0,mid_hi), pos_colorrange=c(max(0,mid_hi),hi), dscalar_data=task_data)

    for (h in 1:length(hemis)) {
      hemi=hemis[h]
      border_brain = task_brain +
        #ss_border(filename = networks_borders[h],border=networks_view_borders, hemisphere = hemi, width=5, colors=control_net_colors) +
        ss_border(filename = networks_borders[2],border=networks_view_borders, hemisphere = hemi, width=5, colors=c('black')) +
        
        ss_view(side="lateral", rotation=c(30,60,90))
      output_file=paste0(output_folder,subject,'_',task,'_',hemi,'_lateral.png')
      captureBrain(border_brain, hemisphere = hemi, filename = output_file)
    }
  }
}


# task activations X GROUP networks
output_folder="/projects/b1081/NSF_HUBS/images/manuscript/verified/lpfc_networks_tasks_groupnet/"
if (!dir.exists(output_folder)) dir.create(output_folder)
#tasks=c("viswm","audwm","msit","vmsit","audattn","visattn","spatialwm","verbalwm")
tasks=c("lang")
#networks_view_borders=c("Frontoparietal","DorsalAttention","CinguloOpercular_Actionmode")
networks_view_borders=c("Language")
control_net_colors = c("#AD9100","#0AF014","#5A28A0")
hemis=c("left")
networks_borders_l=paste0("/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc_L.border")
networks_borders_r=paste0("/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc_R.border")
networks_borders = c(networks_borders_l,networks_borders_r)
for (s in 1:length(subjects)) {
  subject=subjects[s]
  for (t in 1:length(tasks)) {
    task=tasks[t]
    task_file = paste0("/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-", subject, "/domain_summaries/alt/sub-", subject, "_", task, "_zstats_mean_alt.dscalar.nii")
    
    task_data=importCifti(task_file,data_only=TRUE)
    task_data[lpfc_mask==0] = 0
    lpfc_task_data= task_data[lpfc_mask == 1 & lh_idx ==1]
    print(mean(lpfc_task_data <= 1.5) * 100)
    print(quantile(lpfc_task_data,0.75))
    pos_lpfc_task_data=lpfc_task_data[lpfc_task_data > 0]
    
    hi=quantile(lpfc_task_data, 0.98)
    mid_hi=max(quantile(lpfc_task_data, 0.75),0)
    thresh_file=file(paste0(output_folder,subject,'_',task,'.txt'))
    writeLines(as.character(c(hi,mid_hi)),thresh_file)
    task_brain = group_brain +
      ss_dscalar(colorbar=colors, show="pos", pos_thresh = max(0,mid_hi), pos_colorrange=c(max(0,mid_hi),hi), dscalar_data=task_data)

    for (h in 1:length(hemis)) {
      hemi=hemis[h]
      border_brain = task_brain +
      ss_border(filename = networks_borders[h],border=networks_view_borders, hemisphere = hemi, width=5, colors=c("black"))
      #ss_border(filename = networks_borders[h],border=networks_view_borders, hemisphere = hemi, width=4)
      
      ss_view(side="lateral", rotation=c(30,60,90))
      output_file=paste0(output_folder,subject,'_',task,'_',hemi,'_lateral.png')
      captureBrain(border_brain, hemisphere = hemi, filename = output_file)
    }
  }
}