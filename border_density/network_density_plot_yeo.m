clear all;

pfc_mask = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii');
networks=ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/yeo17_group_prior.dscalar.nii');
network_density=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/yeo/group_association_network_density_avg.dscalar.nii']);

group_lpfc_network_density = nanmean(network_density.data(pfc_mask.data ==1));
outdir=['/projects/b1081/NSF_HUBS/images/manuscript/lpfc_association_network_density_yeo'];
if ~exist(outdir, 'dir') mkdir(outdir); end
%%

subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'};

for s =1:length(subjects)
    subject=subjects{s};

    % load hubs cifti 
    network_path = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/yeo/'];
    networks = ft_read_cifti_mod([network_path 'Bipartite_PhysicalCommunities+FinalLabeling_Yeo17.dlabel.nii']);

    networks.data = networks.data(1:59412,1);
    networks.data(pfc_mask.data==0) = 0;

    network_density=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/density/yeo/' subject '_association_network_density_avg.dscalar.nii']);

    lpfc_network_density(s) = nanmean(network_density.data(find(networks.data >0)));
    
end

%%

mean_lpfc_network_density = nanmean(lpfc_network_density);
std_err_lpfc_network_density = nanstd(lpfc_network_density)/sqrt(length(subjects));
std_lpfc_network_density = nanstd(lpfc_network_density);


figure(1)
set(gcf, 'Position', [100, 100, 800, 800]); % Set figure size (width x height in pixels)

bar([group_lpfc_network_density mean_lpfc_network_density], 'FaceColor', [0.8, 0.8, 0.8])
hold on
errorbar([group_lpfc_network_density mean_lpfc_network_density], [0 std_err_lpfc_network_density], 'vertical', 'Color', 'k', 'LineStyle', 'none');
xticklabels({'Group', 'Individual'});
ylabel('Mean LPFC Density');
ylim([0 max(lpfc_network_density)+0.1])
set(gca, 'FontSize', 18);

% Number of data points for each network
num_points = length(lpfc_network_density);
% X coordinates for the scatter plot, slightly jittered for better visibility
x = repmat(2, num_points, 1) + (rand(num_points, 1) - 0.5) * 0.1;
% Y coordinates are the actual data points
y = lpfc_network_density;
% Plotting individual data points
scatter(x, y, 100, 'filled', 'MarkerEdgeColor', 'black', 'MarkerFaceColor','#EDB120'); % Unfilled markers with black outlines
hold off
set(gca, 'FontSize', 24)

saveas(gcf,[outdir '/group_vs_individual_association_bar.jpg'], 'jpg')

%% important stats 

[~, p, ~, stats]= ttest(lpfc_network_density,group_lpfc_network_density);
tStats = stats.tstat;
cohensd = stats.tstat/sqrt(stats.df+1);
fid = fopen(fullfile(outdir, 'indiv_vs_group_density_yeo.txt'), 'w');
fprintf(fid, 'p = %.6f\nt = %.6f\nd = %.6f\n', p, tStats, cohensd);
fclose(fid);
disp('hi')

