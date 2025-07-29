clear all;

subjects = {'HUBS01', 'HUBS02'};
all_halves = {["HUBS01", "half_1"],["HUBS01", "half_2"],["HUBS02", "half_1"], ["HUBS02", "half_2"]};
tasks = {'audattn', 'visattn', 'vmsit', 'msit', 'audwm', 'viswm', 'verbalwm', 'spatialwm'};

task_half_data = [];
template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_cortexonly.dscalar.nii');
for t = 1:length(tasks)
        task=tasks{t};
    for i = 1:length(all_halves)
        subject = char(all_halves{i}(1));
        half = char(all_halves{i}(2));
        task_file = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_zstats.dtseries.nii'];
        task_data =ft_read_cifti_mod(task_file);
        if(strcmp('half_1', half))
            task_half_data(i,:) = task_data.data(1:59412,size(task_data.data,2)-2);
            template.data = task_data.data(1:59412,size(task_data.data,2)-2);
            ft_write_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_' half '_zstats_mean.dscalar.nii'],template);
            make_alt(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_' half '_zstats_mean.dscalar.nii'],['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_' half '_zstats_mean_alt.dscalar.nii']);
        else
            task_half_data(i,:) = task_data.data(1:59412,size(task_data.data,2)-1);
            template.data = task_data.data(1:59412,size(task_data.data,2)-1);
            ft_write_cifti_mod(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_' half '_zstats_mean.dscalar.nii'],template);
            make_alt(['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_' half '_zstats_mean.dscalar.nii'],['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' subject '/task_summaries/sub-' subject '_' task '_' half '_zstats_mean_alt.dscalar.nii']);     
         end
    end
    all_corrs(t,:,:) = corr(task_half_data');
end
%%

cats = {'Sub-1 Half 1 vs. Half 2', 'Sub-2 Half 1 vs. Half 2', 'Sub 1 Half 1 vs. Sub 2 Half 1', 'Sub 1 Half 2 vs. Sub-2 Half 2'};
relevant_idx = [5 15 9 14];

mean_corr=squeeze(nanmean(all_corrs,1));
sem_corr=squeeze(nanstd(all_corrs,1))/sqrt(size(all_corrs,1));

fig = figure(1);
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]
bar(mean_corr(relevant_idx), 'FaceColor', [0.8, 0.8, 0.8]); % RGB values for light gray

hold on 
xticklabels(cats);
ylabel('Correlation')

errorbar(mean_corr(relevant_idx), sem_corr(relevant_idx), 'k', 'linestyle', 'none');
fontsize(fig, 18, "points")

task_colors = [255 149 202; 
    255 66 161; 
    212 24 118;
    151 14 83;
    238 34 12;
    181 23 0;
    255 100 78;
    255 150 141;
    ];  % red, green, blue

task_colors = task_colors/255;


for i = 1:size(all_corrs,1)
    num_points = length(relevant_idx);
    % X coordinates for the scatter plot, slightly jittered for better visibility
    %x = repmat(i, num_points, 1) + (rand(num_points, 1) - 0.5) * 0.1;
    x = 1:4;
    x = x + 0.15*(rand - 0.5);
    % Y coordinates are the actual data points
    y = all_corrs(i,relevant_idx);
    % Plotting individual data points
    scatter(x, y, 50,task_colors(i,:), 'filled'); % Unfilled markers with black outlines
end
%legend(tasks,'Location','southeast')
%title('Task Activity Split Half Correlation')
outdir='/projects/b1081/NSF_HUBS/images/manuscript/split_half_tasks';
if ~exist(outdir,'dir'), mkdir(outdir); end

saveas(figure(1),[outdir '/split_half_corr_bar.jpg'], 'jpg')

%% stats
% Compute stats
corr_relevant = all_corrs(:, relevant_idx);
corr_within = mean(corr_relevant(:, [1 2]), 2);
corr_between = mean(corr_relevant(:, [3 4]), 2);

mean_corr_within = mean(corr_within);
std_corr_within = std(corr_within);

mean_corr_between = mean(corr_between);
std_corr_between = std(corr_between);

[~, p, ~, stats] = ttest(corr_within, corr_between);
p_values = p;
t_statistics = stats.tstat;
cohens_d = stats.tstat / sqrt(stats.df + 1);

% Create table
T = table(mean_corr_within, std_corr_within, ...
          mean_corr_between, std_corr_between, ...
          t_statistics, p_values, cohens_d, ...
          'VariableNames', {'MeanWithin', 'StdWithin', ...
                            'MeanBetween', 'StdBetween', ...
                            'TStatistic', 'PValue', 'CohensD'});

% Write to text file
writetable(T, [outdir '/split_half_task_stats.txt'], 'Delimiter', '\t');


