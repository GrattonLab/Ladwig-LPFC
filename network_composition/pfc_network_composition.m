clear all;
% load in pfc networks 
mask_type = 'lpfc';
template= ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dtseries.nii');
medialwall_idx=find(template.brainstructure ==-1);
surface_area = ft_read_cifti_mod('/projects/b1081/Atlases/32k_ConteAtlas_v2_distribute/Conte69.vertex_areas.dscalar.nii');
surface_area.data(medialwall_idx) = [];
%outdir='/projects/b1081/NSF_HUBS/images/manuscript/lpfc_network_composition';
outdir='/projects/b1081/NSF_HUBS/images/manuscript/verified/lpfc_network_composition_yeo';

if ~exist(outdir,'dir'), mkdir(outdir); end

%% lynch nets
network_ids = [1 2 4 7 8 9 10 11 12 17];
network_labels = {'DMN-B', 'DMN-A','VIS-stream', 'FP', 'DAN', 'Premotor', 'LANG', 'SAL', 'CO'};

colors_idx = [
    234 51 35; %DMN-B
    255 255 218; %DMN-A
    47 117 181; %VIS-stream
    254 255 84; %FP 
    99 214 63; %DAN
    255 130 255; %PRremotor
    64 153 153 %Lang
    0 0 0; %SAL 
    70 7 147; %CO 
    ];

colors_idx=colors_idx/255;

%% yeo nets 

network_ids = [5 6 7 8 11 12 13 14 15 16 17];
network_labels = {'VIS-Stream', 'DAN', 'CO','SAL','PMN', 'FP-A', 'FP-B', 'LANG', 'DN-A', 'DN-B1', 'DN-B2'};

colors_idx = [
    74 155 60; %vis-stream
    0 118 14; %DAN
    196 58 250; %CO
    255 152 213; %SAL 
    119 140 176; %PMN
    230 148 34 %FP-A
    135 50 74; %FP-B 
    12 48 255; %LANG
    0 0 130; %DN-A
    255 255 0; %DMN-B1
    205 62 78; %DMNB-2
    ];

colors_idx=colors_idx/255;

%% group 

%pfc_networks_file=['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc.dscalar.nii'];
pfc_networks_file=['/projects/b1081/NSF_HUBS/resources/yeo17_group_prior_' mask_type '.dscalar.nii'];

networks_struct = ft_read_cifti_mod(pfc_networks_file);
networks =networks_struct.data;
networks(networks==0)=[];

for i =1:length(network_labels)
    group_network_counts(i) = length(find(networks==network_ids(i)));
    group_network_sa(i) = sum(surface_area.data(find(networks==network_ids(i))));
end
    
pct_group_network_counts(:) = 100*group_network_counts/sum(group_network_counts);
pct_group_network_sa(:) = 100*group_network_sa/sum(group_network_sa);

% pct by network 
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'};

for s =1:length(subjects)
    subject=subjects{s};
    disp(subject);

    % load hubs cifti 
    %network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/yeo/'];
    
    %pfc_networks_file = [network_path 'Bipartite_PhysicalCommunities+FinalLabeling_' mask_type '.dlabel.nii'];
    pfc_networks_file = [network_path 'Bipartite_PhysicalCommunities+FinalLabeling_Yeo17_' mask_type '.dlabel.nii'];

    networks_struct = ft_read_cifti_mod(pfc_networks_file);
    networks=networks_struct.data(1:59412);

    networks(networks==0)=[];

    for i =1:length(network_labels)
        all_network_counts(s,i) = length(find(networks==network_ids(i)));
        all_network_sa(s,i) = sum(surface_area.data(find(networks==network_ids(i))));
    end
    
    pct_all_network_counts(s,:) = 100*all_network_counts(s,:)/sum(all_network_counts(s,:));
    pct_all_network_sa(s,:) = 100*all_network_sa(s,:)/sum(all_network_sa(s,:));
end

%%
% Create a table
T = array2table(pct_all_network_sa, ...
    'VariableNames', network_labels(:), ...
    'RowNames', subjects(:));

% Write table to file (tab-delimited text)
writetable(T, [outdir '/surface_area_by_subject.txt'], ...
    'Delimiter', '\t', ...
    'WriteRowNames', true);


%%
mean_pct_all_network_counts=mean(pct_all_network_counts);
std_error_pct_all_network_counts=std(pct_all_network_counts)/sqrt((size(all_network_counts,1)));
stdev_pct_all_network_counts=std(pct_all_network_counts);

mean_pct_all_network_sa=mean(pct_all_network_sa);
std_error_pct_all_network_sa=std(pct_all_network_sa)/sqrt((size(all_network_sa,1)));
stdev_pct_all_network_sa=std(pct_all_network_sa);

%% SA

fig_cat = categorical({'HUBS10', 'HUBS09', 'HUBS08', 'HUBS07', 'HUBS06', 'HUBS05', 'HUBS04', 'HUBS03', 'HUBS02','HUBS01','Group'});
fig_cat = reordercats(fig_cat,cellstr(fig_cat));
pct_fig = [flip(pct_all_network_sa); pct_group_network_sa];

figure(1)
b=barh(fig_cat,pct_fig,'stacked','FaceColor','flat');

for i = 1:size(network_labels,2)
    b(i).CData = colors_idx(i,:);
end
set(gca, 'FontSize', 14);
xlim([0, 100]);
saveas(figure(1),[outdir '/all_' mask_type '_network_composition_sa.jpg'], 'jpg')

%%
figure(4)

mean_pct_all_network_sa=mean(pct_all_network_sa);
std_error_pct_all_network_sa=std(pct_all_network_sa)/sqrt((size(all_network_sa,1)));

hold on
x=linspace(0, 45, 20);
plot(x, x, 'k');
    for i = 1:length(colors_idx)
      errorbar(pct_group_network_sa(i), mean_pct_all_network_sa(i), std_error_pct_all_network_sa(i), 'vertical', 'Color', 'k');
      hold on 
      scatter(pct_group_network_sa(i),mean_pct_all_network_sa(i),[], colors_idx(i,:),"filled","LineWidth",0.5, "MarkerEdgeColor",'k');
      %hold on
      %errorbar(pct_group_network_counts(i), mean_pct_all_network_counts(i), std_error_pct_all_network_counts(i), 'vertical', 'Color', colors_idx(i,:));
    end
xlabel('Group Average');
ylabel('Individuals');
%title('Percent of LPFC Surface Area')
xlim([0, max(pct_group_network_sa+5)]);
ylim([0, max(pct_group_network_sa)+5]);
set(gca, 'FontSize', 14);

saveas(figure(4),[outdir '/' mask_type 'group_vs_individual_sa.jpg'], 'jpg')

%% yeo bar chart 

combined_group = sum(pct_group_network_sa(:,[6 7]));

% yeo only
mean_combined_fp = mean(sum(pct_all_network_sa(:,[6 7]),2));
std_error_combined_fp = std(sum(pct_all_network_sa(:,[6 7]),2))/sqrt(size(all_network_sa,1));
stddev_combined_fp = std(mean(sum(pct_all_network_sa(:,[6 7]),2)));

means = [combined_group mean_combined_fp pct_group_network_sa(6) mean_pct_all_network_sa(6) pct_group_network_sa(7) mean_pct_all_network_sa(7)];
errors = [0 std_error_combined_fp 0 std_error_pct_all_network_sa(6) 0 std_error_pct_all_network_sa(7)];

fig=figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]

b=bar(means);
hold on
columns=categorical({'Combined FP Group', 'Combined FP Individual', 'FP-A Group', 'FP-A Individual', 'FP-B Group', 'FP-B Individual'});
xticklabels(columns);
er = errorbar([2 4 6], means([2 4 6]), errors([2 4 6]), 'k.', 'LineWidth', 1.2);
b.FaceColor = 'flat';
b.CData(1,:) = [0.8 0.8 0.8];
b.CData(2,:) = [0.8 0.8 0.8];
b.CData(3,:) = colors_idx(6,:);
b.CData(4,:) = colors_idx(6,:);
b.CData(5,:) = colors_idx(7,:);
b.CData(6,:) = colors_idx(7,:);

ylabel('Percent of LPFC Surface Area')
fontsize(fig, 20, "points")
ylim([0 max(means)+2])
saveas(gcf,[outdir '/' mask_type '_fp_bar_charts.jpg'], 'jpg')

%% yeo stats just for combined
% Define names and relevant networks
relevant_nets = [6 7];
relevant_net_names = {'FP-A', 'FP-B', 'Combined'};

% Preallocate
n = length(relevant_net_names);
p = zeros(n,1);
tStats = zeros(n,1);
cohensd = zeros(n,1);

% Run t-tests
for i = 1:2
    [~, p(i), ~, stats] = ttest(pct_all_network_sa(:,relevant_nets(i)),pct_group_network_sa(relevant_nets(i)));
    tStats(i) = stats.tstat;
    cohensd(i) = stats.tstat / sqrt(stats.df + 1);
end

% Combined FP-A + FP-B
[~, p(3), ~, stats] = ttest(sum(pct_all_network_sa(:,relevant_nets),2),sum(pct_group_network_sa(relevant_nets)));
tStats(3) = stats.tstat;
cohensd(3) = stats.tstat / sqrt(stats.df + 1);

% Create and save table
T = table(relevant_net_names(:), p, tStats, cohensd, ...
          'VariableNames', {'Network', 'PValue', 'TStatistic', 'CohensD'});
writetable(T, [outdir '/fp_surfacearea_yeo.txt'], 'Delimiter', '\t');



%% stats 
for i = 1:length(network_labels)
    [h(i), p(i), ci, stats]= ttest(pct_all_network_sa(:,i),pct_group_network_sa(i));
    cilo(i) = ci(1);
    cihi(i)= ci(2);
    tStats(i) = stats.tstat;
    cohensd(i) = stats.tstat/sqrt(stats.df+1);
end

% Preallocate
h = zeros(length(network_labels),1);
p = zeros(length(network_labels),1);
cilo = zeros(length(network_labels),1);
cihi = zeros(length(network_labels),1);
tStats = zeros(length(network_labels),1);
cohensd = zeros(length(network_labels),1);
mean_indiv = zeros(length(network_labels),1);
std_indiv = zeros(length(network_labels),1);
group_val = zeros(length(network_labels),1);

% Loop through networks
for i = 1:length(network_labels)
    data_indiv = pct_all_network_sa(:,i);
    data_group = pct_group_network_sa(i);

    [h(i), p(i), ci, stats] = ttest(data_indiv, data_group);
    cilo(i) = ci(1);
    cihi(i) = ci(2);
    tStats(i) = stats.tstat;
    cohensd(i) = stats.tstat / sqrt(stats.df + 1);
    mean_indiv(i) = mean(data_indiv);
    std_indiv(i) = std(data_indiv);
    group_val(i) = data_group;
end

% Build results table
T = table(network_labels(:), mean_indiv(:), std_indiv(:), group_val(:), ...
          tStats(:), p(:), cilo(:), cihi(:), cohensd(:), ...
          'VariableNames', {'Network', 'MeanIndividual', 'StdIndividual', ...
                            'GroupValue', 'TStatistic', 'PValue', ...
                            'CI_Low', 'CI_High', 'CohensD'});

% Write to file
outfile = fullfile(outdir, 'lpfc_network_ttest_stats.txt');
writetable(T, outfile, 'Delimiter', '\t');


