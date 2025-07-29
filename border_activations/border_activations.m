

%% ANALYSIS 1: If you are in the LPFC of FP/DAN/CO, are activated regions closer to the 2 other 2 control nets than non-activated regions?
clear all;
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08','HUBS09', 'HUBS10'};
tasks = {'spatialwm', 'verbalwm', 'viswm', 'audwm', 'vmsit', 'msit', 'visattn', 'audattn'};

nets = [7 8 12];
thresh = 1;
pfc_mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii']); 
pfc_idx=pfc_mask.data == 1;
hi_threshs = [inf thresh];
lo_threshs = [thresh -inf];

%threshs_names= {'95', '90', '85', '80', '75', '70', '65', '60', '60', '55'};
for s= 1:length(subjects)    
    subject = subjects{s};
    disp(subject);
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
    dist_to_net = ft_read_cifti_mod([network_path subject '_distance_to_networks_group_dmat.dtseries.nii']);
    networks_struct = ft_read_cifti_mod([network_path 'Bipartite_PhysicalCommunities+FinalLabeling_lpfc.dlabel.nii']);
    for n =1:length(nets)
        target_net = nets(n);
        border_net = nets;
        border_net(n)=[];
        target_network_idx = find(networks_struct.data(1:59412,1) == target_net);
        for task = 1:length(tasks)
            task_data=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' tasks{task} '_zstats_mean.dscalar.nii']);
            for t =1:length(hi_threshs)
                hi_thresh = hi_threshs(t);
                lo_thresh = lo_threshs(t);
                active_target_network_idx = intersect(find(task_data.data > lo_thresh & task_data.data < hi_thresh),target_network_idx);
                active_dist(s,n,task,t) = mean(min(dist_to_net.data(active_target_network_idx,border_net(1)),dist_to_net.data(active_target_network_idx,border_net(2))));
                %active_dist(s,n,t) = mean(mean([dist_to_net.data(active_target_network_idx,border_net(1)) dist_to_net.data(active_target_network_idx,border_net(2))]));
            end
        end
    end
end

%% plot
outdir = '/projects/b1081/NSF_HUBS/images/manuscript/verified/border_activations_control';
if ~isfolder(outdir), mkdir(outdir); end

for task = 1:length(tasks)
    sub_means = squeeze(mean(active_dist(:,:,task,:),2));
    means = squeeze(nanmean(sub_means));
    stderrs= squeeze(nanstd(sub_means))/sqrt(10);
    
    x = 1:length(means);
    figure(task)
    set(gcf, 'Position', [100, 100, 400, 1000]);  % Set figure size and position
    bar(x, means, 'FaceColor', [0.8 0.8 0.8],'EdgeColor', 'k', 'LineWidth', 1.5);  % Customize color if desired
    hold on;
    errorbar(x, means, stderrs, 'vertical','LineStyle', 'none', 'Color', 'k', 'LineWidth', 2);
    hold off;
    xlim([0.5, length(means)+0.5]);
    disp(xlim);
    xticks(x);
    xticklabels({['Z>' num2str(thresh)], ['Z<' num2str(thresh)]});  % Replace with your cell array of labels if you have one
    ylim([0 9])
    box off;
    [~, p(task), ~, stats]= ttest(sub_means(:,1), sub_means(:,2));
    tStats(task) = stats.tstat;
    cohensd(task) = stats.tstat/sqrt(stats.df+1);
    set(gca, 'FontSize', 40);  % sets font size for current axes
    saveas(gcf, [outdir '/border_activity_' tasks{task} '_' num2str(thresh) '.jpg']);
end

% Create table
T = table(tasks', p', tStats', cohensd',  ...
    'VariableNames', {'Task', 'p', 't', 'Cohens_d'});

% Write table to file
writetable(T, [outdir '/border_activation_stats_' num2str(thresh) '.txt'], 'Delimiter', '\t');

%%


%         active_target_network_idx = intersect(find(task_data.data > active_thresh),target_network_idx);
%         inactive_target_network_idx = intersect(find(task_data.data < active_thresh),target_network_idx);
%         active_size(s,n) = size(active_target_network_idx,1);
%         inactive_size(s,n) = size(inactive_target_network_idx,1);
%         active_dist(s,n) = mean(min(dist_to_net.data(active_target_network_idx,border_net(1)),dist_to_net.data(active_target_network_idx,border_net(2))));
%         inactive_dist(s,n) = mean(min(dist_to_net.data(inactive_target_network_idx,border_net(1)),dist_to_net.data(inactive_target_network_idx,border_net(2))));
%     end
% end

mean_active_dist = mean(active_dist(:));
mean_inactive_dist = mean(inactive_dist(:));
sem_active_dist = std(active_dist(:))/sqrt(10);
sem_inactive_dist = std(inactive_dist(:))/sqrt(10);

figure()
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size (width x height in pixels)

bar([mean_active_dist mean_inactive_dist],'FaceColor', [0.8, 0.8, 0.8])
hold on 
errorbar([mean_active_dist mean_inactive_dist], [sem_active_dist sem_inactive_dist],'LineStyle','none', 'Color','k');
ylabel('Distance to Control Nets')
xticklabels({'Active', 'Non-active'})
set(gca, 'FontSize', 28);

mean_active_size = mean(active_size);
mean_inactive_size = mean(inactive_size);
saveas(gcf, ['/projects/b1081/NSF_HUBS/images/manuscript/border_activations_control/group_nets_FP_border_activity_DANCO_' num2str(active_thresh) '.jpg']);

% stats 
[~, p, ~, stats]= ttest(mean(inactive_dist,2), mean(active_dist,2));
tStats = stats.tstat;
cohensd = stats.tstat/sqrt(stats.df+1);



%% ANALYSIS 2 - For each network compare the distance of top 10% vertices to mean distance from border 

% Q1: What is the distance of a top 10% vtx in FP LPFC to each network
% compared to the mean distance from net to net 

% mean dist of net to net 
clear all;
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08','HUBS09', 'HUBS10'};
target_net = 1;
border_nets = [8 10 2 12 11 7]; % 
target_net_name = 'DNB';
%tasks = {'viswm', 'audwm', 'spatialwm', 'verbalwm', 'vmsit', 'msit', 'visattn', 'audattn','md'};
tasks = {'tom'};
for t=1:length(tasks)
    task = tasks{t};
    for s= 1:length(subjects)    
        subject = subjects{s};   
        disp(subject);
        network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/'];
        dist_to_net = ft_read_cifti_mod([network_path subject '_distance_to_networks.dtseries.nii']);
        networks_struct = ft_read_cifti_mod([network_path 'Bipartite_PhysicalCommunities+FinalLabeling_lpfc.dlabel.nii']);
        target_network_idx = find(networks_struct.data(1:59412,1) == target_net); 
        task_data=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/summary/sub-' subject '_' task '_zstats_mean.dscalar.nii']);
        task_data.data(setdiff(1:numel(task_data.data), target_network_idx)) = -Inf;
        [sorted_vals, sorted_vtx] = sort(task_data.data, 'descend');
        top10vtx = sorted_vtx(1:0.2*size(target_network_idx,1));
        for b = 1:length(border_nets)
            net_dist(s,b) = mean(dist_to_net.data(target_network_idx,border_nets(b)));
            task_dist(s,b) = mean(dist_to_net.data(top10vtx,border_nets(b)));
        end
    end

    ratio_dist = task_dist./net_dist;
    mean_ratio_dist = mean(ratio_dist,1);
    sem_ratio_dist = std(ratio_dist,1)/sqrt(10);
    
    figure(t)
    bar(mean_ratio_dist)
    hold on
    errorbar(mean_ratio_dist,sem_ratio_dist,'LineStyle','none', 'Color','k');
    yline(1)
    %xticklabels({'DAN','CO','DN-B','DN-A','SAL','LANG'})
    %xticklabels({'DAN','CO','FP','DN-A','SAL','LANG'})
    %xticklabels({'DAN','CO','DN-B','DN-A','SAL','FP'})
    xticklabels({'DAN','CO','DN-A','LANG','SAL','FP'})
    
    outdir='/projects/b1081/NSF_HUBS/images/manuscript/border_activations_bar';
    if ~isfolder(outdir), mkdir(outdir); end
    saveas(figure(t),[outdir '/distance_lpfc_' target_net_name '_' task '.jpg'], 'jpg')
    %saveas(figure(t),[outdir '/distance_lpfc_' target_net_name '_' task '.jpg'], 'jpg')
    
    for i=1:length(mean_ratio_dist)
        [~, p(i), ~, stats]= ttest(ratio_dist(:,i),1);
        corrected_pstats(i) = p(i)*length(mean_ratio_dist);
        tStats(i) = stats.tstat;
        cohensd(i) = stats.tstat/sqrt(stats.df+1);
    end
end

%% ANALYSIS 3 - If you are in the top 10% of activated LPFC vertices of cog control, how far are you from the borders?
clear all;
subjects = {'HUBS01','HUBS02','HUBS03','HUBS04','HUBS05','HUBS06','HUBS07','HUBS08','HUBS09','HUBS10'};

for s= 1:length(subjects)    
    subject = subjects{s};
    disp(subject);
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/'];
    dist_to_net = ft_read_cifti_mod([network_path subject '_distance_to_networks.dtseries.nii']);
    networks_struct = ft_read_cifti_mod([network_path 'Bipartite_PhysicalCommunities+FinalLabeling_lpfc.dlabel.nii']);
    %target_network_idx = find(networks_struct.data(1:59412,1) == target_net); 
    task_data=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/summary/sub-' subject '_md_zstats_mean.dscalar.nii']);
    task_data.data(networks_struct.data == 0) = 0;
    [~, sorted_idx] = sort(task_data.data, 'descend');
    top10_vtx = sorted_idx(1:934);
    nets = networks_struct.data(top10_vtx);
    border_dists = [dist_to_net.data(top10_vtx,8) dist_to_net.data(top10_vtx,7) dist_to_net.data(top10_vtx,10)];
    a = min(border_dists);
    disp('hi');
end

%%
clear all;
subjects = {'HUBS01','HUBS02','HUBS03','HUBS04','HUBS05','HUBS06','HUBS07','HUBS08','HUBS09','HUBS10'};
border_net = 12; 
for s = 1:length(subjects)
    subject = subjects{s};
    bordering_region_idx{s} = zeros(59412,1);
    dir = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/'];
    dist_to_net = ft_read_cifti_mod([dir subject '_distance_to_networks.dtseries.nii']);
    target_net = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/FP_regions.dtseries.nii']); 
    for i = 1:size(target_net.data,2)
        region_idx = find(target_net.data(:,i) > 0);
        dist_to_border = min(dist_to_net.data(region_idx,border_net));
        %disp(dist_to_border);
        if dist_to_border <=1
            bordering_region_idx{s}(region_idx) =1;
        end
    end
end

distance_to_border = 0:2:20;
for s= 1:length(subjects)    
    subject = subjects{s};
    disp(subject);
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/'];
    dist_to_net = ft_read_cifti_mod([network_path subject '_distance_to_networks.dtseries.nii']);
    task_data=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/summary/sub-' subject '_md_zstats_mean.dscalar.nii']);
    target_idx = find(bordering_region_idx{s} ==1);
    for d = 1:length(distance_to_border)-1
        idx = target_idx(dist_to_net.data(target_idx,border_net) > distance_to_border(d) & dist_to_net.data(target_idx,border_net) <= distance_to_border(d+1)); %dist to dan
        num_idx(s,d) = length(idx);
        mean_activation(s,d)=mean(task_data.data(idx));
    end
end
%%
figure(1)
mean_activation_dist = nanmean(mean_activation,1);
sem_activation_dist = nanstd(mean_activation,1)/sqrt(10);

bar(mean_activation_dist)
hold on
errorbar(mean_activation_dist,sem_activation_dist);
figure(2)
mean_num_idx = mean(num_idx,1);
bar(mean_num_idx);

%%
 %% ANALYSIS 1: If you are in the LPFC of FP/DAN/CO, are activated regions closer to the 2 other 2 control nets than non-activated regions?
 % activated = z > 1, non-active = z < 1
clear all;
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08','HUBS09', 'HUBS10'};
tasks = {'spatialwm', 'verbalwm', 'viswm', 'audwm', 'vmsit', 'msit', 'visattn', 'audattn'};

target_net = 7;
border_net = [8 10];
%border_net = 1;
active_thresh = 2;
inactive_thresh = 2;
dist_to_net = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/lynch_group_distance_to_networks.dtseries.nii']);
pfc_mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii']); 
pfc_idx=pfc_mask.data == 1;

for s= 1:length(subjects)    
    subject = subjects{s};
    disp(subject);
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
    dist_to_net = ft_read_cifti_mod([network_path subject '_distance_to_networks_group_dmat.dtseries.nii']);
    networks_struct = ft_read_cifti_mod([network_path 'Bipartite_PhysicalCommunities+FinalLabeling_lpfc.dlabel.nii']);
    
    target_network_idx = find(networks_struct.data(1:59412,1) == target_net); 
    task_data=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_md_zstats_mean.dscalar.nii']);
    active_thresh = prctile(task_data.data(pfc_idx),75);
    disp(active_thresh);
    active_target_network_idx = intersect(find(task_data.data > active_thresh),target_network_idx);
    inactive_target_network_idx = intersect(find(task_data.data < active_thresh),target_network_idx);
    active_size(s) = size(active_target_network_idx,1);
    inactive_size(s) = size(inactive_target_network_idx,1);
    active_dist(s) = mean(min(dist_to_net.data(active_target_network_idx,border_net(1)),dist_to_net.data(active_target_network_idx,border_net(2))));
    inactive_dist(s) = mean(min(dist_to_net.data(inactive_target_network_idx,border_net(1)),dist_to_net.data(inactive_target_network_idx,border_net(2))));
    %active_dist(s) = mean(dist_to_net.data(active_target_network_idx,border_net));
    %inactive_dist(s) = mean(dist_to_net.data(inactive_target_network_idx,border_net));
end

mean_active_dist = mean(active_dist);
mean_inactive_dist = mean(inactive_dist);
sem_active_dist = std(active_dist)/sqrt(10);
sem_inactive_dist = std(inactive_dist)/sqrt(10);

figure()
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size (width x height in pixels)

bar([mean_active_dist mean_inactive_dist],'FaceColor', [0.8, 0.8, 0.8])
hold on 
errorbar([mean_active_dist mean_inactive_dist], [sem_active_dist sem_inactive_dist],'LineStyle','none', 'Color','k');
ylabel('Distance to CO or DAN')
xticklabels({'Active', 'Non-active'})
set(gca, 'FontSize', 28);

mean_active_size = mean(active_size);
mean_inactive_size = mean(inactive_size);
saveas(gcf, ['/projects/b1081/NSF_HUBS/images/manuscript/border_activations_control/group_nets_FP_border_activity_DANCO_' num2str(active_thresh) '.jpg']);

% stats 
[~, p, ~, stats]= ttest(inactive_dist, active_dist);
tStats = stats.tstat;
cohensd = stats.tstat/sqrt(stats.df+1);

%% ANALYSIS independant variable = how far from other control nets y axis = mean z-activity 
clear all;
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08','HUBS09', 'HUBS10'};
nets = [7 8 12];
pfc_mask=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii']); 
pfc_idx=pfc_mask.data == 1;
distances_lo = [0 2 4 6 8 10];
distances_hi = [2 4 6 8 10 12];
task = 'verbalwm';
for s= 1:length(subjects)    
    subject = subjects{s};
    disp(subject);
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
    dist_to_net = ft_read_cifti_mod([network_path subject '_distance_to_networks_group_dmat.dtseries.nii']);
    networks_struct = ft_read_cifti_mod([network_path 'Bipartite_PhysicalCommunities+FinalLabeling_lpfc.dlabel.nii']);
    for n =1:length(nets)
        target_net = nets(n);
        border_net = nets;
        border_net(n)=[];
        target_network_idx = find(networks_struct.data(1:59412,1) == target_net); 
        task_data=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_zstats_mean.dscalar.nii']);
        for d = 1:length(distances_hi)
            dist_to_control = min(dist_to_net.data(target_network_idx,border_net(1)), dist_to_net.data(target_network_idx,border_net(2)));
            idx = target_network_idx(find(dist_to_control <= distances_hi(d) & dist_to_control > distances_lo(d)));
            avg_z(s,n,d) = mean(task_data.data(idx));
            sizes(s,n,d) = length(idx);
        end
    end
end

sub_means = squeeze(mean(avg_z,2));
means = squeeze(nanmean(sub_means));
stderrs= squeeze(nanstd(sub_means))/sqrt(10);
distances = {'0-2','2-4', '4-6', '6-8' '8-10', '10-12'};

x = 1:length(means);
figure;
bar(x, means, 'FaceColor', [0.8 0.8 0.8]);  % Customize color if desired
hold on;
errorbar(x, means, stderrs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
hold off;
xlim([0.5, length(means)+0.5]);
xticks(x);
xticklabels(distances);  % Replace with your cell array of labels if you have one
ylabel('Mean Z Value');
xlabel('Distance from Control Border');
box off;
saveas(gcf, ['/projects/b1081/NSF_HUBS/images/manuscript/border_activations_control/activity_x_from_border_' task '.jpg']);
