clear all;
% set variables 
template= ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dtseries.nii');
bad_idx_mask = ft_read_cifti_mod('/projects/b1081/Atlases/WashU120_SNR/WashU120_avg_SNR_thresh700.dtseries.nii');
bad_idx = find(bad_idx_mask.data == 1); % identify low SNR vertices
limits = [6 8 10 12 14]; % mm distances
medialwall_idx=find(template.brainstructure ==-1); 

%association_networks = [1 2 7 8 10 11 12]; % lynch association networks 
association_networks = [6 7 8 11 12 13 14 15 16 17]; %yeo assoctiation networks

% load a distance matrix 
load('/projects/b1081/Scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_large_uint8.mat'); % load dmat
D = distances(1:59412,1:59412);
clear distances;

%% Generate Per Subject Border Density

subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
all_networks = [];
for s=1:length(subjects)
    subject = subjects{s};
    disp(subject);
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/yeo/'];
    %network_file = [network_path 'Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii'];       
    network_file = [network_path 'Bipartite_PhysicalCommunities+FinalLabeling_Yeo17.dlabel.nii'];       

    networks = ft_read_cifti_mod(network_file);
    networks.data = networks.data(1:59412);
    all_networks(s,:) = networks.data(1:59412,1);
end

for s = 1:length(subjects)
    subject = subjects{s};
    outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/density/yeo/'];
    if ~isfolder(outFolder), mkdir(outFolder); end

    disp(subject);

    for l=1:length(limits) % generate density per distance [6-14]
        for v=1:size(all_networks,2)
            if(ismember(v,bad_idx))
                sizes(s,v,l) = NaN;
                unique_nets(s,v,l)= NaN;
                unique_nets_norm(s,v,l)=NaN;
            else
                nearby_vtx = [];
                nearby_vtx = find(D(v,:) <= limits(l));
                sizes(s,v,l) = length(nearby_vtx);
                unique_nets(s,v,l) = length(unique(all_networks(s,nearby_vtx)));
                unique_assn_nets(s,v,l) = max(sum(ismember(unique(all_networks(s,nearby_vtx)),association_networks)), 0.000001);
            end
        end
    end

    unique_nets(s,:,l+1) = mean(unique_nets(s,:,1:l),3);
    unique_assn_nets(s,:,l+1) = mean(unique_assn_nets(s,:,1:l),3);

    template.data = squeeze(unique_nets(s,:,:));
    template.time = 1:size(template.data,2);
    template.hdr.dim(6)=size(template.data,2);
    ft_write_cifti_mod([outFolder subject '_network_density.dtseries.nii'],template);

    template.data = squeeze(unique_assn_nets(s,:,:));
    template.time = 1:size(template.data,2);
    template.hdr.dim(6)=size(template.data,2);
    ft_write_cifti_mod([outFolder subject '_association_network_density.dtseries.nii'],template); 
end

mean_unique_nets = squeeze(mean(unique_nets,1)); %average across limits
mean_unique_assn_nets = squeeze(mean(unique_assn_nets,1)); %average across limits

outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/yeo/'];
if ~isfolder(outFolder), mkdir(outFolder); end

template.data = [mean_unique_assn_nets];
template.time = 1:size(template.data,2);
template.hdr.dim(6)=size(template.data,2);
ft_write_cifti_mod([outFolder 'avg_association_network_density.dtseries.nii'],template);

template.data = [mean_unique_nets];
template.time = 1:size(template.data,2);
template.hdr.dim(6)=size(template.data,2);
ft_write_cifti_mod([outFolder 'avg_network_density.dtseries.nii'],template);

%% Calculate group network density

%networks = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/lynch_group_prior50.dlabel.nii');
networks = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/yeo17_group_prior.dlabel.nii');
all_networks = networks.data(1:59412)';

outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/yeo/'];
if ~isfolder(outFolder), mkdir(outFolder); end

for l =1:length(limits)
    for v=1:size(all_networks,2)
        if(ismember(v,bad_idx))
            sizes(v,l) = NaN;
            unique_nets(v,l)= NaN;
            unique_assn_nets(v,l)= NaN;

        else
            nearby_vtx = [];
            nearby_vtx = find(D(v,:) <= limits(l));
            unique_nets(v,l) = length(unique(all_networks(1,nearby_vtx)));
            unique_assn_nets(v,l) = sum(ismember(unique(all_networks(1,nearby_vtx)),association_networks));
        end
    end
end

unique_nets(:,l+1) = mean(unique_nets(:,1:l),2);
unique_assn_nets(:,l+1) = mean(unique_assn_nets(:,1:l),2);

template.data = unique_nets;
template.time = 1:size(template.data,2);
template.hdr.dim(6)=size(template.data,2);
ft_write_cifti_mod([outFolder 'group_network_density.dtseries.nii'],template);

template.data = unique_assn_nets;
template.data(template.data == 0) = 0.00001;
template.time = 1:size(template.data,2);
template.hdr.dim(6)=size(template.data,2);
ft_write_cifti_mod([outFolder 'group_association_network_density.dtseries.nii'],template);

%% Calculate Average Subject - Group Density
clear all;
outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/yeo/'];

template= ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dtseries.nii');

group = ft_read_cifti_mod([outFolder 'group_network_density.dtseries.nii']);
sub_avg = ft_read_cifti_mod([outFolder 'avg_network_density.dtseries.nii']);

template.data = sub_avg.data - group.data;
template.data(template.data == 0) = 0.00001;
template.time = 1:size(template.data,2);
template.hdr.dim(6)=size(template.data,2);
ft_write_cifti_mod([outFolder 'avg-group_network_density.dscalar.nii'],template);

group = ft_read_cifti_mod([outFolder 'group_association_network_density.dtseries.nii']);
sub_avg = ft_read_cifti_mod([outFolder 'avg_association_network_density.dtseries.nii']);

template.data = sub_avg.data - group.data;
template.data(template.data == 0) = 0.00001;
template.time = 1:size(template.data,2);
template.hdr.dim(6)=size(template.data,2);
ft_write_cifti_mod([outFolder 'avg-group_association_network_density.dscalar.nii'],template);

%%
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');

for s = 1:length(subjects)
    subject = subjects{s};
    outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/density/yeo/'];

    data = ft_read_cifti_mod([outFolder subject '_network_density.dtseries.nii']);
    template.data = data.data(:,size(data.data,2));
    ft_write_cifti_mod([outFolder subject '_network_density_avg.dscalar.nii'],template)
    make_alt([outFolder subject '_network_density_avg.dscalar.nii'], [outFolder subject '_network_density_avg_alt.dscalar.nii']);
end

subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');

for s = 1:length(subjects)
    subject = subjects{s};
    outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/density/yeo/'];

    data = ft_read_cifti_mod([outFolder subject '_association_network_density.dtseries.nii']);
    template.data = data.data(:,size(data.data,2));
    ft_write_cifti_mod([outFolder subject '_association_network_density_avg.dscalar.nii'],template)
    make_alt([outFolder subject '_association_network_density_avg.dscalar.nii'], [outFolder subject '_association_network_density_avg_alt.dscalar.nii']);
end

%%
groups = {'group', 'avg', 'avg-group'};
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');
outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/yeo/'];

for g = 1:length(groups)
    group = groups{g};
    data = ft_read_cifti_mod([outFolder group '_network_density.dtseries.nii']);
    template.data = data.data(:,size(data.data,2));
    ft_write_cifti_mod([outFolder group '_network_density_avg.dscalar.nii'],template)
    make_alt([outFolder group '_network_density_avg.dscalar.nii'], [outFolder group '_network_density_avg_alt.dscalar.nii']);
end

%%
groups = {'group', 'avg', 'avg-group'};
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');
outFolder = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/yeo/'];

for g = 1:length(groups)
    group = groups{g};
    data = ft_read_cifti_mod([outFolder group '_association_network_density.dtseries.nii']);
    template.data = data.data(:,size(data.data,2));
    ft_write_cifti_mod([outFolder group '_association_network_density_avg.dscalar.nii'],template)
    make_alt([outFolder group '_association_network_density_avg.dscalar.nii'], [outFolder group '_association_network_density_avg_alt.dscalar.nii']);
end