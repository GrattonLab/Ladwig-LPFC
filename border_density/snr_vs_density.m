clear all;
% load snr file 
snr=ft_read_cifti_mod('/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/sub-group_allsess_task-rest_desc-mode1000_mean_LR_surf_32k_fsLR.dtseries.nii');
% load density 
density = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-group/density/lynch/avg_association_network_density_avg.dscalar.nii');
% pfc only 
group_networks_lpfc = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc.dlabel.nii']);
pfc_mask = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/group_lpfc_mask.dscalar.nii');
bad_idx_mask = ft_read_cifti_mod('/projects/b1081/Atlases/WashU120_SNR/WashU120_avg_SNR_thresh700.dtseries.nii');
bad_idx = find(bad_idx_mask.data == 1);

%%
snr_lpfc = snr.data(pfc_mask.data > 0 & bad_idx_mask.data ~= 1);
density_lpfc = density.data(pfc_mask.data > 0 & bad_idx_mask.data ~= 1);
%bad_idx = density_lpfc == 0;

corr(density_lpfc,snr_lpfc)
%%
x = snr_lpfc(:);
y = density_lpfc(:);

figure('Position', [100, 100, 800, 600])
scatter(x, y, 60, 'filled')
xlabel('Mean BOLD Signal (Mode 1000 Normalized) ', 'FontSize', 24)
ylabel('Association Network Density', 'FontSize', 24)
set(gca, 'FontSize', 24)
hold on

% Fit and plot regression line
p_fit = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(p_fit, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 2)

% Compute correlation and confidence interval
[r_matrix, ~, rlo_matrix, rup_matrix] = corrcoef(x, y);
r_val = r_matrix(2);         % correlation value
rlo = rlo_matrix(2);         % lower CI bound
rup = rup_matrix(2);         % upper CI bound

% Display r and CI on plot
text_pos_x = min(x) + 0.05 * range(x);
text_pos_y = max(y) - 0.03 * range(y);
text(text_pos_x, text_pos_y, ...
    sprintf('r = %.3f, 95%% CI = [%.3f, %.3f]', r_val, rlo, rup), ...
    'FontSize', 24, 'FontWeight', 'bold');

box on
outdir='/projects/b1081/NSF_HUBS/images/manuscript/snr';
% 
if ~isfolder(outdir), mkdir(outdir); end
saveas(gcf,[outdir '/snr_vs_density_scatter.jpg'], 'jpg')
