clear all;
% get their PFC task responses 
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
tasks = {'epiproj', 'tom', 'lang'};
target_nets = {'DN-A', 'DN-B', 'LANG'};
pfc_mask = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii');

%%
for t = 1:length(tasks)
    task=tasks{t};
    all_responses={};
    task_responses=[];
    disp(task);
    for s = 1:length(subjects)
        subject=subjects{s};
        %task_data = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0_GSR/sub-' subject '/domain_summary/sub-' subject '_' task '_zstats_mean.dscalar.nii']);
        task_data = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/domain_summaries/sub-' subject '_' task '_zstats_mean.dscalar.nii']);
        
        task_data_z = task_data.data(1:59412,size(task_data.data,2));
        network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
        network_file = [network_path 'Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii'];       
        networks = ft_read_cifti_mod(network_file);
        networks.data = networks.data(1:59412);
        networks.data(pfc_mask.data==0) = 0;

        dan_idx= find(networks.data==8);
        fp_idx=find(networks.data==7);
        co_idx=find(networks.data==12);
        dmna_idx=find(networks.data==2);
        dmnb_idx=find(networks.data==1);
        lang_idx=find(networks.data==10);
        sal_idx=find(networks.data==11);
   
        dan_response(s)=mean(task_data_z(dan_idx));
        fp_response(s)=mean(task_data_z(fp_idx));      
        co_response(s)=mean(task_data_z(co_idx));     
        dmna_response(s)=mean(task_data_z(dmna_idx));
        dmnb_response(s)=mean(task_data_z(dmnb_idx));      
        lang_response(s)=mean(task_data_z(lang_idx)); 
        sal_response(s)=mean(task_data_z(sal_idx)); 
    end
    
    columns=categorical({'FP', 'DAN', 'CO', 'DN-A', 'DN-B', 'LANG', 'SAL/PMN'});

    network_means=[nanmean(fp_response) nanmean(dan_response) nanmean(co_response) nanmean(dmna_response) nanmean(dmnb_response) nanmean(lang_response), nanmean(sal_response)];
    network_stddevs=[nanstd(fp_response) nanstd(dan_response) nanstd(co_response) nanstd(dmna_response) nanstd(dmnb_response) nanstd(lang_response) nanstd(sal_response)];
  
    network_stderr = network_stddevs/sqrt(10);
        
    fig=figure(t);
    set(gcf, 'Position', [100, 100, 500, 800]); % Set figure size (width x height in pixels)

    hb=bar(network_means);
    
    hold on
    er=errorbar(network_means,network_stderr,'LineStyle','none', 'Color','k', 'LineWidth', 1);
   
    colors =[
        254 255 84; % FP   
        99 214 63; % DAN
        70 7 147; % CO
        255 255 218;% DMN-A
        234 51 35; % DMN-B
        64 153 153; %LANG
        0 0 0 % SAL/PMN
        ];
    
    colors=colors/255;
    
    fontsize(fig, 28, "points")
    
    for k = 1:size(colors,1)
        hb.FaceColor = 'flat';
        hb.CData(k,:) = colors(k,:);
    end
    %yticks(min(network_means-network_stderr)-0.1:0.5:max(network_means+network_stderr) + 0.3);
    yticks(-1.5:0.5:1.5);

    xticklabels(columns);
    %ylabel('LPFC mean Z value')
    %title([task ' LPFC responses'])
    ylim([min(min(network_means - network_stderr) - 0.1,0), max(network_means+network_stderr) + 0.3]);
    ylim([-2 2])
    outdir='/projects/b1081/NSF_HUBS/images/manuscript/verified/networks_lpfc_task_response_baronly';
    if ~isfolder(outdir), mkdir(outdir); end
    saveas(figure(t),[outdir '/lpfc_response_z_' task '.jpg'], 'jpg')

    %all_responses = {fp_response, dan_response, co_response, dmna_response, dmnb_response, lang_response, sal_response, vis_response, aud_response};
    %all_values = [fp_response dan_response co_response dmna_response dmnb_response lang_response sal_response vis_response, aud_response];

    all_responses = {fp_response, dan_response, co_response, dmna_response, dmnb_response, lang_response, sal_response};
%     all_values = [fp_response dan_response co_response dmna_response dmnb_response lang_response sal_response];
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

    statsT = table( ...
    strings(0,1), ... % NetworkName
    zeros(0,1), ...   % PValue
    zeros(0,1), ...   % CorrectedPValue
    zeros(0,1), ...   % TStat
    zeros(0,1), ...   % CohensD
    'VariableNames', {'NetworkName', 'PValue', 'CorrectedPValue', 'TStat', 'CohensD'});

         % for specialized
         target_idx = find(columns == target_nets{t});
         for i = 1:length(all_responses)
                % Perform paired t-test (one-tailed)
                [~, p, ~, stats] = ttest(all_responses{target_idx}, all_responses{i}, 'Tail', 'right');
                p_value = round(p,8); % Store the p-value
                corrected_p_value = round(p*(length(columns)-1),8); % bonferroni correction
                tstat = round(stats.tstat,4);
                cohensd = round(stats.tstat/sqrt(stats.df+1),4);
                stats_row = {string(columns(i)),p_value,corrected_p_value,tstat,cohensd};
                statsT = [statsT; stats_row];
         end
 writetable(statsT, [outdir '/' task '_stats.txt'], 'Delimiter', '\t');
 end
