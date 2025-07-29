% get dice values for all LPFC networks 
% get LPFC split half netowrks
clear all;

subjects = {'HUBS01', 'HUBS02'};
all_halves = {["HUBS01", "splithalf1"],["HUBS01", "splithalf2"],["HUBS02", "splithalf1"], ["HUBS02", "splithalf2"]};

for i = 1:length(all_halves)
    subject = char(all_halves{i}(1));
    half = char(all_halves{i}(2));
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/' half '/pfm/lynch/'];
    network_file = [network_path 'Bipartite_PhysicalCommunities+FinalLabeling_lpfc.dlabel.nii'];

    networks = ft_read_cifti_mod(network_file);
    networks.data = networks.data(1:59412);
    networks_data(i,:) = networks.data;
    networks.data(networks.data == 13) = 11;
    lpfc_networks{i} = networks.data(networks.data > 0);
    networks_data(i,:) = networks.data;
end


%%

relevant_networks = [1,2,7,8,10,11,12];

for i = 1:length(relevant_networks)
    net=relevant_networks(i);
    for j = 1:size(networks_data,1)
        net1_idx = networks_data(j,:)==net;
        for k =1:size(networks_data,1)
            net2_idx=networks_data(k,:)==net;
            dice_all(i,j,k) = dice(net1_idx,net2_idx);
        end
    end
end


relevant_idx = [5 15 9 14];
mean_dice=squeeze(nanmean(dice_all,1));
sem_dice=squeeze(nanstd(dice_all,1))/sqrt(size(dice_all,1));

cats = {'Sub-1 Half 1 vs. Half 2', 'Sub-2 Half 1 vs. Half 2', 'Sub 1 Half 1 vs. Sub 2 Half 1', 'Sub 1 Half 2 vs. Sub-2 Half 2'};

fig = figure(1);
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]
bar(mean_dice(relevant_idx), 'FaceColor', [0.8, 0.8, 0.8]); % RGB values for light gray

hold on 
xticklabels(cats);
ylabel('Dice Coefficient')

errorbar(mean_dice(relevant_idx), sem_dice(relevant_idx), 'k', 'linestyle', 'none');
fontsize(fig, 18, "points")

colors =[
    234 51 35; % DMN-B
    255 255 218;% DMN-A
    254 255 84; % FP 
    99 214 63; % DAN
    64 153 153; %LANG
    0 0 0; % SAL 
    70 7 147; % CO
    ];


colors=colors/255;

for i = 1:size(dice_all,1)
    num_points = length(relevant_idx);
    % X coordinates for the scatter plot, slightly jittered for better visibility
    %x = repmat(i, num_points, 1) + (rand(num_points, 1) - 0.5) * 0.1;
    x = 1:4;
    x = x + 0.1*(rand - 0.5);
    % Y coordinates are the actual data points
    y = dice_all(i,relevant_idx);
    % Plotting individual data points
    scatter(x, y, 50, colors(i,:), 'filled'); % Unfilled markers with black outlines
end

ylim([0 1]);

%title('LPFC Network Split Half Overlap')
outdir='/projects/b1081/NSF_HUBS/images/manuscript/split_half_networks_lpfc';
if ~exist(outdir,'dir'), mkdir(outdir); end
saveas(figure(1),[outdir '/split_half_dice_bar.jpg'], 'jpg')

%% stats 

dice_relevant = dice_all(:,relevant_idx);
dice_within = mean(dice_relevant(:,[1 2]),2);
dice_between = mean(dice_relevant(:, [3 4]),2);

mean_dice_within = mean(dice_within);
std_dice_within = std(dice_within);

mean_dice_between = mean(dice_between);
std_dice_between = std(dice_between);

[~, p, ~, stats] = ttest(dice_within,dice_between);
p_values = p;
t_statistics = stats.tstat;
cohens_d = stats.tstat/sqrt(stats.df+1);

% Create table with named columns
T = table(mean_dice_within, std_dice_within, ...
          mean_dice_between, std_dice_between, ...
          t_statistics, p_values, cohens_d, ...
          'VariableNames', {'MeanWithin', 'StdWithin', ...
                            'MeanBetween', 'StdBetween', ...
                            'TStatistic', 'PValue', 'CohensD'});

% Write table to text file (tab-delimited for Excel)
writetable(T, [outdir '/split_half_network_dice_stats.txt'], 'Delimiter', '\t');


%% not using

load('/projects/b1081/NSF_HUBS/resources/better_jet_colormap.mat');
final = [mean_dice(5) mean_dice(9); mean_dice(14) mean_dice(15)];
h=heatmap(round(final,2));
colormap(h, better_jet_colormap_diff); % Set colormap to 'parula'
caxis([0 1])
h.XData = {'Sub-1 Half-2', 'Sub-2 Half-2'};
h.YData = {'Sub-1 Half-1', 'Sub-2 Half-1'};
title('Mean Dice Overlap of LPFC networks')


