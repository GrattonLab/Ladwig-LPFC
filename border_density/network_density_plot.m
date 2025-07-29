clear all;
network_ids = [1 2 4 7 8 9 10 11 12 17];
network_labels = {'DMN-B', 'DMN-A','VIS-stream', 'FP', 'DAN', 'Premotor', 'LANG', 'SAL', 'CO'};
outdir=['/projects/b1081/NSF_HUBS/images/manuscript/verified/lpfc_association_network_density'];
if ~exist(outdir, 'dir') mkdir(outdir); end

colors_idx = [
    234 51 35; %DMN-B
    255 255 218; %DMN-A
    47 117 181; %VIS-stream
    254 255 84; %FP 
    99 214 63; %DAN
    255 130 255; %Premotor
    64 153 153 %Lang
    0 0 0; %SAL 
    70 7 147; %CO 
    ];

colors_idx=colors_idx/255;

%% Get Per Network Density For the Group
pfc_mask = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii');
networks=ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/lynch_group_prior20.dscalar.nii');
network_density=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/lynch/group_association_network_density_avg.dscalar.nii']);

pfc_density = network_density.data(pfc_mask.data ==1);

networks.data = networks.data(1:59412,1);
networks.data(pfc_mask.data==0) = 0;

for i =1:length(network_labels)
    group_network_density(i) = nanmean(network_density.data(find(networks.data==network_ids(i))));
end

group_lpfc_network_density = nanmean(network_density.data(find(networks.data >0)));

%%
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'};

for s =1:length(subjects)
    subject=subjects{s};

    % load hubs cifti 
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/'];
    networks = ft_read_cifti_mod([network_path 'Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii']);

    networks.data = networks.data(1:59412,1);
    networks.data(pfc_mask.data==0) = 0;

    network_density=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/density/lynch/' subject '_association_network_density_avg.dscalar.nii']);

    for i =1:length(network_labels)
        all_network_density(s,i) = nanmean(network_density.data(find(networks.data==network_ids(i))));
    end

    lpfc_network_density(s) = nanmean(network_density.data(find(networks.data >0)));
    
end
%% group vs individual bar 

mean_lpfc_network_density = nanmean(lpfc_network_density);
std_err_lpfc_network_density = nanstd(lpfc_network_density)/sqrt(length(subjects));
std_lpfc_network_density = nanstd(lpfc_network_density);


%% Make Figure

fig=figure();
set(gcf, 'Position', [100, 100, 1000, 800]); % Set figure size (width x height in pixels)

bar([group_lpfc_network_density mean_lpfc_network_density], 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on
%errorbar([group_lpfc_network_density mean_lpfc_network_density], [0 std_err_lpfc_network_density], ...
   % 'vertical', 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 10);
%errorbar([group_lpfc_network_density mean_lpfc_network_density], [0 std_err_lpfc_network_density], 'vertical', 'Color', 'k', 'LineStyle', 'none');
xticklabels({'Group', 'Individual'});
ylabel('Mean LPFC Density');

set(gca, 'FontSize', 38);

% Number of data points for each network
num_points = length(lpfc_network_density);
% X coordinates for the scatter plot, slightly jittered for better visibility
x = repmat(2, num_points, 1) + (rand(num_points, 1) - 0.5) * 0.1;
% Y coordinates are the actual data points
y = lpfc_network_density;
% Plotting individual data points
scatter(x, y, 100, 'filled', 'MarkerEdgeColor', 'black', 'MarkerFaceColor','#EDB120', 'MarkerEdgeColor', 'k','LineWidth', 2); % Unfilled markers with black outlines
hold off
saveas(gcf,[outdir '/group_vs_individual_association_bar.jpg'], 'jpg')

%% per network graph - i didnt use this
figure(2)
set(gcf, 'Position', [100, 100, 800, 800]); % Set figure size (width x height in pixels)

mean_all_network_density=nanmean(all_network_density);
std_error_all_network_density=nanstd(all_network_density/sqrt((size(all_network_density,1))));

hold on
x=linspace(2, 5, 20);
plot(x, x, 'k');
    for i = 1:length(colors_idx)
      errorbar(group_network_density(i), mean_all_network_density(i), std_error_all_network_density(i), 'vertical', 'Color', 'k');
      hold on 
      scatter(group_network_density(i),mean_all_network_density(i), 60, colors_idx(i,:),"filled","LineWidth",0.5, "MarkerEdgeColor",'k');
    end
xlabel('Group Average');
ylabel('Individuals');
set(gca, 'FontSize', 20);
ylim([2 5])

%title('Percent of LPFC Surface Area')
saveas(figure(2),[outdir '/group_vs_individual_networks.jpg'], 'jpg')

%% stats 

[h, p, ci, stats] = ttest(lpfc_network_density, group_lpfc_network_density);

% Compute t-stat and Cohen's d
tStats = stats.tstat;
cohensd = stats.tstat / sqrt(stats.df + 1);

% Create table
T = table(p, tStats, cohensd, ci(1), ci(2), ...
    'VariableNames', {'p', 't', 'Cohens_d', 'CI_lower', 'CI_upper'});

% Write table to text file
writetable(T, [outdir '/indiv_vs_group_density_stats.txt'], 'Delimiter', '\t');

