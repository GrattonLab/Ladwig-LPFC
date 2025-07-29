% do all the rotations - i later did this by batch job so ignore 
clear all;
subjects = {'HUBS01','HUBS02','HUBS03','HUBS04','HUBS05','HUBS06','HUBS07','HUBS08','HUBS09','HUBS10'};
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dtseries.nii');
co_data= readtable('/projects/b1081/NSF_HUBS/resources/CO_regions.txt');
for s = 1:length(subjects)
    template.data=zeros(size(template.data));
    subject = subjects{s};  
    disp(subject);
    idx = co_data.Var2(strcmp(co_data.Var1, subject));
    co=ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/CO_regions.dtseries.nii']);
    co_idx = find(co.data(:,idx) ==1);
    template.data(co_idx) = 1;
    ft_write_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_CO_region.dtseries.nii'], template)
    generate_cifti_rotations(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/rostral_CO_region.dtseries.nii'],10000,['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/'],6);
end

%%
clear all;
association_networks = [1 2 7 8 10 11 12];

subjects = {'HUBS01','HUBS03','HUBS04','HUBS05','HUBS06','HUBS07','HUBS08','HUBS09','HUBS10'};
all_density = [];
 for s = 1:length(subjects)
     subject = subjects{s}; 
     disp(subject);
     networks = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/pfm/lynch/Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii']);       
     rotations = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/functional_masks/rostral_CO_region_ROTATIONS_10000perms.dtseries.nii']);
     density = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/density/lynch/' subject '_association_network_density_avg.dscalar.nii']);
     co_region = ft_read_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subject '/functional_masks/rostral_CO_region.dtseries.nii']);
     co_idx = find(co_region.data ==1);
     co_density(s) = mean(density.data(co_idx,1));
     for i = 1:size(rotations.data,2)
         idx = find(rotations.data(:,i) ==1);
         association_pct(s,i) = sum(ismember(networks.data(idx),association_networks))/length(idx); 
         num_networks(s,i) = numel(unique(networks.data(idx)));
         all_density(s,i) = mean(density.data(idx,1));
     end
     clear density;
 end
%%


%%
% Generate control data (example normal distribution)
valid_spins=all_density(find(num_networks == 1 & association_pct > 0));

fig = figure();
set(gcf, 'Position', [100, 100, 800, 600]); % Set figure size (width x height in pixels)
co_color= [93 37 160];
co_color = co_color/255;
control_values = valid_spins;

real_values = co_density;

% Plot histogram for control values
%histogram(b, 'FaceColor', 'k', 'EdgeColor', 'k', 'Normalization', 'pdf'); 
%histogram(b, 'FaceColor', 'k', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinMethod', 'auto');
histogram(valid_spins, 'FaceColor', 'k', 'EdgeColor', 'k', 'Normalization', 'probability', 'NumBins', 10);

hold on;

% Overlay real values as red dots
%scatter(real_values, y_max(2) * 0.95 * ones(size(real_values)), 100, 'r', 'filled');

real_val = mean(co_density); % Replace with your desired X value
se = std(co_density)/sqrt(size(co_density,2));
p = mean(valid_spins >= real_val);

hold on; % Keep the current plot
fill([real_val-se real_val+se real_val+se real_val-se], ...
     [0 0 max(ylim) max(ylim)], ...
     co_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none'); % semi-transparent band
% Optional: display p-value
text(real_val, max(ylim)*0.85, ...
     {'Rostral CO', ['p = ' num2str(round(p,3))]}, ...
     'Color', co_color, 'HorizontalAlignment', 'center', 'FontSize', 12);
hold off; % Release the hold

% Labels and title
xlabel('Mean Association Network Density');
ylabel('Proportion of Spin Test Samples');
hold off;
fontsize(fig, 22, "points")

outdir='/projects/b1081/NSF_HUBS/images/manuscript/rostral_co_rotation_analysis';
if ~isfolder(outdir), mkdir(outdir); end
saveas(gcf,[outdir '/co_rotation_vs_distribution.jpg'], 'jpg')
