clear all;
% get their PFC task responses 
tasks = {'tom'};
control = {'spatialwm'};

%tom 
test_subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS07', 'HUBS09', 'HUBS10'}; 
%lang
%test_subjects = {'HUBS03', 'HUBS04', 'HUBS05', 'HUBS07', 'HUBS08', 'HUBS10'}; 

% DMN variants - HUBS06 HUBS08
% HUBS06 = 2 
% HUBS08 = 1 

% LANG variants - HUBS01, HUBS02, HUBS06, HUBS09

subject = 'HUBS01';
task = 'lang';
variant = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/functional_masks/LANG_variant.dscalar.nii']); %edit for ur variant
variant_idx = find(variant.data(:,1) > 0);
target_task_data = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/domain_summaries/sub-' subject '_' task '_zstats_mean.dscalar.nii']);
spatialwm_task_data = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_spatialwm_zstats_mean.dscalar.nii']);
variant_target = mean(target_task_data.data(variant_idx));
variant_spatialwm = mean(spatialwm_task_data.data(variant_idx));
variant_diff = mean(target_task_data.data(variant_idx))-mean(spatialwm_task_data.data(variant_idx));
for s = 1:length(test_subjects)
    test_subject = test_subjects{s};
    target_task_data = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' test_subject '/domain_summaries/sub-' test_subject '_' task '_zstats_mean.dscalar.nii']);
    spatialwm_task_data = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' test_subject '/task_summaries/sub-' test_subject '_spatialwm_zstats_mean.dscalar.nii']);
    sub_variant_response(s) = mean(target_task_data.data(variant_idx));
    sub_variant_spatialwm_response(s) = mean(spatialwm_task_data.data(variant_idx));
    sub_variant_diff(s) = sub_variant_response(s)-sub_variant_spatialwm_response(s);
end

group_mean_target_response=mean(sub_variant_response);
group_mean_spatialwm_response=mean(sub_variant_spatialwm_response);
group_mean_diff = mean(sub_variant_diff);
group_sem_target_response=std(sub_variant_response)/sqrt(length(sub_variant_response));
group_sem_spatialwm_response=std(sub_variant_spatialwm_response)/sqrt(length(sub_variant_response));
group_sem_diff = std(sub_variant_diff)/sqrt(length(sub_variant_response));

fig = figure(1);
set(fig, 'Position', [100, 100, 500, 800]);  % Set figure size

% Data
data = [variant_diff, group_mean_diff];
errors = [0, group_sem_diff];

% Bar plot
b = bar(data, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on

% Error bars
er = errorbar(1:2, data, errors, 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 10);

% Labels and formatting
xticks(1:2)
xticklabels({'Ind', 'Group'})
xtickangle(45)
ylabel('Mean ROI Z value')
ylim([min(data - errors) - abs(0.1*min(data - errors)), max(data + errors) + 0.1*max(data + errors)])
fontsize(fig, 40, "points")
box off

% Save
outdir = '/projects/b1081/NSF_HUBS/images/manuscript/verified/variant_responses';
if ~isfolder(outdir), mkdir(outdir); end
saveas(fig, fullfile(outdir, ['sub-' subject '_LANG.jpg']), 'jpg');



%%
figure()
bar([variant_target group_mean_target_response variant_spatialwm group_mean_spatialwm_response]);
hold on
errorbar([variant_target group_mean_target_response variant_spatialwm group_mean_spatialwm_response],[0 group_sem_target_response 0 group_sem_spatialwm_response], 'LineStyle','none', 'Color','k', 'LineWidth', 1);
xticklabels(categorical({'Ind - TOM', 'Group -TOM', 'Ind - SpatialWM', 'Group - SpatialWM'}))
ylabel('z-statistic')





%% make variants

%%
%make LPFC target regions
clear all;

%HUBS01 2 - but subtract 
%HUBS02 3 dont subtract 
%HUBS06 1 yes subtract 
%HUBS09 3 dont subtract 
subject = 'HUBS06';
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');
template.data=zeros(size(template.data));
group = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/lynch_group_prior20.dlabel.nii');
group_net_idx = find(group.data == 10);
variant = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/LANG_regions.dtseries.nii']);
variant_idx = find(variant.data(:,1) > 0);
clean_variant_idx = setdiff(variant_idx,group_net_idx);
template.data(clean_variant_idx) = 1;
ft_write_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/LANG_variant.dscalar.nii'],template);
