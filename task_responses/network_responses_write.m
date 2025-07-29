clear all;
% get their PFC task responses 
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
tasks = {'viswm', 'audwm', 'spatialwm', 'verbalwm', 'vmsit', 'msit', 'visattn', 'audattn'};
pfc_mask = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii');
%tasks = {'md', 'epiproj', 'tom', 'lang'};
networks_id = [8 7 12 2 1 10 11];
network_names = {'DAN', 'FP', 'CO', 'DNA', 'DNB', 'LANG', 'SAL'};
%% calculate all responses
T = table(string.empty, string.empty, string.empty, [], ...
    'VariableNames', {'subject', 'task', 'network', 'z'});

for t = 1:length(tasks)
    task=tasks{t};
    all_responses={};
    task_responses=[];
    disp(task);
    for s = 1:length(subjects)
        subject=subjects{s};
        task_data = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_zstats_mean.dscalar.nii']);
        task_data_z = task_data.data(1:59412,size(task_data.data,2));
        network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
        network_file = [network_path 'Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii'];       
        networks = ft_read_cifti_mod(network_file);
        networks.data = networks.data(1:59412);
        networks.data(pfc_mask.data==0) = 0;
        for n = 1:length(networks_id)
            idx = find(networks.data == networks_id(n));
            z_stat = mean(task_data_z(idx));
            stat_row = {subject, task, network_names{n}, z_stat};
            T(end+1, :) = stat_row;
        end
    end
end
output_dir = '/projects/b1081/NSF_HUBS/results';
if ~exist(output_dir, 'dir') 
    mkdir(output_dir) 
end

save('/projects/b1081/NSF_HUBS/results/task_responses.mat', 'T');
clear all;
%% plot responses (i dont actually use this) 
clear all;
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
tasks = {'viswm', 'audwm', 'spatialwm', 'verbalwm', 'vmsit', 'msit', 'visattn', 'audattn'};
%tasks = {'md', 'epiproj', 'tom', 'lang'};

load('/projects/b1081/NSF_HUBS/results/task_responses.mat');
plot_nets = {'FP', 'DAN', 'CO', 'LANG','DN-A', 'DN-B', 'SAL'};
%columns=categorical({'FP', 'DAN', 'CO'});

for n =1:length(plot_nets)
    figure(n); 
    hold on;
    T_use = T(T.network == plot_nets{n},:);
    T_use.task = categorical(T_use.task, tasks, 'Ordinal', true);

    for i = 1:numel(subjects)

        subj_rows = T_use.subject == subjects{i};
        T_subj = T_use(subj_rows, :);
        
        % Sort by task to ensure consistent x-axis
        [~, idx] = sort(T_subj.task);
        plot(T_subj.task(idx), T_subj.z(idx), '-o', 'DisplayName', subjects{i});
        xlabel('Task');
        ylabel('Z');
        legend('Location', 'eastoutside');

        title([plot_nets{n} ' Responses Across Tasks']);
    end
    yline(0, '--k');
    hold off
end

%     network_means=[nanmean(fp_response) nanmean(dan_response) nanmean(co_response)];
%     network_stddevs=[nanstd(fp_response) nanstd(dan_response) nanstd(co_response)];
%     network_stderr = network_stddevs/sqrt(10);
%         
%     fig=figure(t);
%     set(gcf, 'Position', [100, 100, 800, 800]); % Set figure size (width x height in pixels)
%     hb=bar(network_means);
%     
%     hold on
%     er=errorbar(network_means,network_stderr,'LineStyle','none', 'Color','k', 'LineWidth', 1);
%    
%     colors =[
%         254 255 84; % FP   
%         99 214 63; % DAN
%         70 7 147; % CO
%         ];
%     
%     colors=colors/255;
%     fontsize(fig, 24, "points")
%     
%     for k = 1:size(colors,1)
%         hb.FaceColor = 'flat';
%         hb.CData(k,:) = colors(k,:);
%     end
%     %yticks(min(network_means-network_stderr)-0.1:0.5:max(network_means+network_stderr) + 0.3);
%     yticks(-1.5:0.5:1.5);
% 
%     xticklabels(columns);
%     ylabel('mean z-value')
%     title([task ' LPFC responses'])
%     ylim([min(min(network_means - network_stderr) - 0.1,0), max(network_means+network_stderr) + 0.3]);
%     %ylim([-2 2])
%     outdir='/projects/b1081/NSF_HUBS/images/manuscript/networks_lpfc_task_response_baronly';
%     if ~isfolder(outdir), mkdir(outdir); end
%     saveas(figure(t),[outdir '/lpfc_response_z_controlnet' task '.jpg'], 'jpg')
% 
%     all_responses = {fp_response, dan_response, co_response};
%     all_values = [fp_response dan_response co_response];
% 
%     for i = 1:length(all_responses)
%         % Number of data points for each network
%         num_points = length(all_responses{i});
%         % X coordinates for the scatter plot, slightly jittered for better visibility
%         x = repmat(i, num_points, 1) + (rand(num_points, 1) - 0.5) * 0.1;
%         % Y coordinates are the actual data points
%         y = all_responses{i};
%         % Plotting individual data points
%         scatter(x, y, 'MarkerEdgeColor', 'k'); % Unfilled markers with black outlines
%     end
% 
%      ylim([min(all_values) - 0.1, max(all_values) + 0.1]);
%      hold off
% 
%      outdir='/projects/b1081/NSF_HUBS/images/manuscript/networks_lpfc_task_response';
%      if ~isfolder(outdir), mkdir(outdir); end
%      saveas(figure(t),[outdir '/lpfc_response_z_' task '.jpg'], 'jpg')
% 
%     statsT = table( ...
%     strings(0,1), ... % NetworkName
%     zeros(0,1), ...   % PValue
%     zeros(0,1), ...   % CorrectedPValue
%     zeros(0,1), ...   % TStat
%     zeros(0,1), ...   % CohensD
%     'VariableNames', {'NetworkName', 'PValue', 'CorrectedPValue', 'TStat', 'CohensD'});  
% 
%     for i = 1:length(all_responses)
%         [~, p, ~, stats] = ttest(all_responses{i}, 0, 'Tail', 'right'); % One-sample, one-tailed test
%         p_value = round(p,8); % Store the p-value
%         corrected_p_value = round(p*(length(columns)),8); % bonferroni correction
%         tstat = round(stats.tstat,4);
%         cohensd = round(stats.tstat/sqrt(stats.df+1),4);
%         stats_row = {string(columns(i)),p_value,corrected_p_value,tstat,cohensd};
%         statsT = [statsT; stats_row];
%     end
%     
%     writetable(statsT, [outdir '/' task '_stats.txt'], 'Delimiter', '\t');
% end
