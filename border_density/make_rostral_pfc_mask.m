%make rostral pfc mask 
clear all;
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dtseries.nii');
subjects = {'HUBS01','HUBS02','HUBS03','HUBS04','HUBS05','HUBS06','HUBS07','HUBS08','HUBS09','HUBS10'};
n=12;
lh_vtx=1:29706;
template.data=zeros(size(template.data));
% define LPFC CO regions 
 for s = 1:length(subjects)
    subject = subjects{s};
    disp(subject);
    dir = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/'];
    lpfc_networks = ft_read_cifti_mod([dir 'Bipartite_PhysicalCommunities+FinalLabeling_lpfc.dlabel.nii']);
    lpfc_networks.data = lpfc_networks.data(lh_vtx); 
    regions = cifti_cluster(lpfc_networks,n,n,50);
    template.data(lh_vtx,1:size(regions,2)) = regions;
    ft_write_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/CO_regions.dtseries.nii'],template);
 end

 %%
%make group rostral pfc mask 
clear all;
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dtseries.nii');
n=12;
lh_vtx=1:29706;
template.data=zeros(size(template.data));
subjects = {'group'};
% define LPFC CO regions 
 for s = 1:length(subjects)
    subject = subjects{s};
    disp(subject);
    lpfc_networks = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc.dlabel.nii']);
    lpfc_networks.data = lpfc_networks.data(lh_vtx); 
    regions = cifti_cluster(lpfc_networks,n,n,50);
    template.data(lh_vtx,1:size(regions,2)) = regions;
    ft_write_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/CO_regions.dtseries.nii'],template);
 end

%%
%manualy pick the rostral CO one [annoying and such a nice element of
%things like the Kong 400

load('/projects/b1081/Scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_large_uint8.mat'); % load dmat
D = distances(1:59412,1:59412);
co_data= readtable('/projects/b1081/NSF_HUBS/resources/CO_regions.txt');
clear distances;
%%
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');
subjects = {'group'};
for s = 1:length(subjects)
    template.data=zeros(size(template.data));
    subject = subjects{s};    
    idx = co_data.Var2(strcmp(co_data.Var1, subject));
    co=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/CO_regions.dtseries.nii']);
    co_idx = find(co.data(:,idx) ==1);
    all_vtx = [];
    for i=1:length(co_idx)
        nearby_vtx = find(D(co_idx(i),:) <= 15);
        all_vtx = [all_vtx nearby_vtx];
    end
    unique_vtx{s} = unique(all_vtx);
    template.data(unique_vtx{s}) =1;
    ft_write_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_lpfc_mask.dscalar.nii'], template);
end

%%
clear all;
subjects = {'HUBS01','HUBS02','HUBS03','HUBS04','HUBS05','HUBS06','HUBS07','HUBS08','HUBS09','HUBS10'};
%subjects = {'group'};
%subjects = {'HUBS04'};
co_data= readtable('/projects/b1081/NSF_HUBS/resources/CO_regions.txt');
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');
rostral_co_lpfc_colorfile = '/projects/b1081/NSF_HUBS/resources/rostral_lpfc_co_colorfile.txt';
surface_l = '/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.L.very_inflated.32k_fs_LR.surf.gii';
surface_r = '/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.R.very_inflated.32k_fs_LR.surf.gii';
for s = 1:length(subjects)
    template.data=zeros(size(template.data));
    subject = subjects{s};    
    idx = co_data.Var2(strcmp(co_data.Var1, subject));
    co=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/CO_regions.dtseries.nii']);
    co_idx = find(co.data(:,idx) ==1);
    template.data(co_idx) =1;
    rostral_co_lpfc = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_co_lpfc.dscalar.nii'];
    rostral_co_lpfc_alt = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_co_lpfc_alt.dscalar.nii'];
    rostral_co_lpfc_alt_label = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_co_lpfc_alt.dlabel.nii'];

    ft_write_cifti_mod(rostral_co_lpfc, template);
    make_alt(rostral_co_lpfc, rostral_co_lpfc_alt);
    make_label_file(rostral_co_lpfc_alt, rostral_co_lpfc_alt_label,rostral_co_lpfc_colorfile);
    make_border_file(rostral_co_lpfc_alt_label,surface_l,surface_r);

end


%%
subjects = {'group'};
surface_l = '/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.L.very_inflated.32k_fs_LR.surf.gii';
surface_r = '/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.R.very_inflated.32k_fs_LR.surf.gii';
rostral_lpfc_colorfile = '/projects/b1081/NSF_HUBS/resources/rostral_lpfc_colorfile.txt';
for s = 1:length(subjects)
    subject = subjects{s};
    rostral_lpfc = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_lpfc_mask.dscalar.nii'];
    rostral_lpfc_label = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_lpfc_mask.dlabel.nii'];
    rostral_lpfc_alt = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_lpfc_mask_alt.dscalar.nii'];
    rostral_lpfc_alt_label = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_lpfc_mask_alt.dlabel.nii'];
    
    make_alt(rostral_lpfc, rostral_lpfc_alt);
    make_label_file(rostral_lpfc_alt, rostral_lpfc_alt_label, rostral_lpfc_colorfile);
    make_label_file(rostral_lpfc_alt, rostral_lpfc_label, rostral_lpfc_colorfile);
    make_border_file(rostral_lpfc_label,surface_l,surface_r);
    make_border_file(rostral_lpfc_alt_label,surface_l,surface_r);
end



%%
sub_num_regions(s,n)=size(cifti_cluster(lpfc_networks,n,n,50),2);

% find anything within 10mm of that thing

% save out mask 