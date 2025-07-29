% LPFC ONLY network files 
clear all;
%template=ft_read_cifti_mod('/projects/b1081/Atlases/brain_template.dscalar.nii');
template=ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');
color_file='/projects/b1081/NSF_HUBS/resources/Lynch_networks_colorfile_edited.txt';
%color_file='/projects/b1081/NSF_HUBS/resources/yeo_networks_colorfile_renamed.txt';
%color_file = '/projects/b1081/NSF_HUBS/resources/kong2019_networks_colorfile_renamed.txt';

group_surface_L = '/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.L.very_inflated.32k_fs_LR.surf.gii';
group_surface_R = '/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.R.very_inflated.32k_fs_LR.surf.gii';
%subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'};
subjects = {'HUBS04', 'HUBS07'};
%subjects = {'group'};
mask_type = 'lpfc';

for s = 1:length(subjects)
    subject=subjects{s};

    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
    %network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/splithalf1/pfm/'];
    %network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/group/pfm/'];

    disp(subject);
    % load networks 
    network_name = 'Bipartite_PhysicalCommunities+FinalLabeling';
    network_data = [network_path network_name '.dlabel.nii'];
    %network_data = ['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20.dlabel.nii'];

    networks = ft_read_cifti_mod(network_data);
    temp_data = networks.data;
    %mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/lpfc_mask.dscalar.nii']);
    mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/group_' mask_type '_mask.dscalar.nii']); 
   
    pfc_idx=mask.data==1;
  
    temp_data(~pfc_idx) = 0;
    temp_data(59413:end) = 0;
    template.data=temp_data(1:59412,1);

    pfc_networks_file = [network_path network_name '_' mask_type '.dscalar.nii'];
    ft_write_cifti_mod(pfc_networks_file,template);

    pfc_networks_label_file = [network_path network_name '_' mask_type '.dlabel.nii'];
    make_label_file(pfc_networks_file,pfc_networks_label_file,color_file);

    pfc_networks_file_alt = [network_path network_name '_' mask_type '_alt.dscalar.nii'];

    make_alt(pfc_networks_file, pfc_networks_file_alt);

    pfc_networks_label_file_alt = [network_path network_name '_' mask_type '_alt.dlabel.nii'];

    make_label_file(pfc_networks_file_alt,pfc_networks_label_file_alt,color_file);

    make_border_file(pfc_networks_label_file, group_surface_L, group_surface_R);

end
    
%% group networks
clear all;
mask_type = 'lpfc';
%color_file='/projects/b1081/NSF_HUBS/resources/yeo_networks_colorfile_renamed.txt';
color_file = '/projects/b1081/NSF_HUBS/resources/kong2019_networks_colorfile_renamed.txt';


group_networks='/projects/b1081/NSF_HUBS/resources/kong2019_group_prior.dscalar.nii';
networks = ft_read_cifti_mod(group_networks);
%networks.data(networks.data == 13) = 11;
pfc_networks=networks;
pfc_mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/group_' mask_type '_mask.dscalar.nii']); 
pfc_idx=pfc_mask.data == 1;
pfc_networks.data(~pfc_idx) = 0;

pfc_networks_file=['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_' mask_type '.dscalar.nii'];
pfc_networks_label_file=['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_' mask_type '.dlabel.nii'];

ft_write_cifti_mod(pfc_networks_file,pfc_networks);

make_label_file(pfc_networks_file,pfc_networks_label_file,color_file);

pfc_networks_file_alt=['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_' mask_type '_alt.dscalar.nii'];
pfc_networks_file_label_alt=['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_' mask_type '_alt.dlabel.nii'];

make_alt(pfc_networks_file, pfc_networks_file_alt);

make_label_file(pfc_networks_file_alt,pfc_networks_file_label_alt,color_file);

%%
clear all;
mask_type = 'lpfc';
color_file='/projects/b1081/NSF_HUBS/resources/yeo_networks_colorfile_renamed.txt';

group_networks='/projects/b1081/NSF_HUBS/resources/yeo17_group_prior.dscalar.nii';
networks = ft_read_cifti_mod(group_networks);
%networks.data(networks.data == 13) = 11;
pfc_networks=networks;
pfc_mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/group_' mask_type '_mask.dscalar.nii']); 
pfc_idx=pfc_mask.data == 1;
pfc_networks.data(~pfc_idx) = 0;

pfc_networks_file=['/projects/b1081/NSF_HUBS/resources/yeo17_group_prior_' mask_type '.dscalar.nii'];
pfc_networks_label_file=['/projects/b1081/NSF_HUBS/resources/yeo17_group_prior_' mask_type '.dlabel.nii'];

ft_write_cifti_mod(pfc_networks_file,pfc_networks);

make_label_file(pfc_networks_file,pfc_networks_label_file,color_file);

pfc_networks_file_alt=['/projects/b1081/NSF_HUBS/resources/yeo17_group_prior_' mask_type '_alt.dscalar.nii'];
pfc_networks_file_label_alt=['/projects/b1081/NSF_HUBS/resources/yeo17_group_prior_' mask_type '_alt.dlabel.nii'];

make_alt(pfc_networks_file, pfc_networks_file_alt);

make_label_file(pfc_networks_file_alt,pfc_networks_file_label_alt,color_file);

%%
%%
clear all;
mask_type = 'lpfc';
color_file='/projects/b1081/NSF_HUBS/resources/kong2019_networks_colorfile_renamed.txt';

group_networks='/projects/b1081/NSF_HUBS/resources/kong2019_group_prior.dscalar.nii';
networks = ft_read_cifti_mod(group_networks);
%networks.data(networks.data == 13) = 11;
pfc_networks=networks;
pfc_mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/group_' mask_type '_mask.dscalar.nii']); 
pfc_idx=pfc_mask.data == 1;
pfc_networks.data(~pfc_idx) = 0;

pfc_networks_file=['/projects/b1081/NSF_HUBS/resources/kong2019_group_prior_' mask_type '.dscalar.nii'];
pfc_networks_label_file=['/projects/b1081/NSF_HUBS/resources/kong2019_group_prior_' mask_type '.dlabel.nii'];

ft_write_cifti_mod(pfc_networks_file,pfc_networks);

make_label_file(pfc_networks_file,pfc_networks_label_file,color_file);

pfc_networks_file_alt=['/projects/b1081/NSF_HUBS/resources/kong2019_group_prior_' mask_type '_alt.dscalar.nii'];
pfc_networks_file_label_alt=['/projects/b1081/NSF_HUBS/resources/kong2019_group_prior_' mask_type '_alt.dlabel.nii'];

make_alt(pfc_networks_file, pfc_networks_file_alt);

make_label_file(pfc_networks_file_alt,pfc_networks_file_label_alt,color_file);
