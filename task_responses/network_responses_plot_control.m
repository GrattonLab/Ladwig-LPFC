clear all;
load('/projects/b1081/NSF_HUBS/results/task_responses.mat');
% Make sure your networks and tasks are in correct order
tasks = unique(T.task, 'stable');  % 8 tasks
networks = {'FP','DAN','CO','DNA','DNB','LANG','SAL'};  % 7 networks
network_names = {'FP','DAN','CO','DN-A','DN-B','LANG','SAL/PMN'};  % 7 networks

% Initialize
mean_z = zeros(length(tasks), length(networks));    % 8x7
stderr_z = zeros(length(tasks), length(networks));

% Compute means and SEM
for i = 1:length(networks)
    for j = 1:length(tasks)
        idx = strcmp(T.network, networks{i}) & strcmp(T.task, tasks{j});
        z_vals = T.z(idx);
        mean_z(j,i) = mean(z_vals);
        stderr_z(j,i) = std(z_vals) / sqrt(length(z_vals));
    end
end

% Colors for each network (7 total)
colors = [
    254 255 84;   % FP   
    99 214 63;    % DAN
    70 7 147;     % CO
    255 255 218;  % DN-A
    234 51 35;    % DN-B
    64 153 153;   % LANG
    0 0 0         % SAL/PMN
] / 255;

% Create grouped bar plot by network
figure('Position', [100, 100, 1000, 900]);  % [left, bottom, width, height]

hold on;
nbars = length(tasks);      % 8 tasks per group
ngroups = length(networks); % 7 networks

groupwidth = min(0.8, nbars/(nbars + 1.5));


for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    h = bar(x, mean_z(i,:), groupwidth / nbars, 'FaceColor', 'flat', ...
        'CData', colors, 'EdgeColor', 'k', 'LineWidth', 0.1);  % <--- added EdgeColor and LineWidth
    errorbar(x, mean_z(i,:), stderr_z(i,:), 'k', 'linestyle', 'none', 'linewidth', 1);
end

xticks(1:ngroups);
xlim([0.3, ngroups + 0.7]);  % Tightens padding on both sides
xticklabels(network_names);
ylim([min(mean_z(:))-max(stderr_z(:)) max(mean_z(:))+max(stderr_z(:))]);
ylabel('Mean LPFC Network Z Value');
%legend(tasks, 'Location', 'eastoutside', 'FontSize', 16);
set(gca, 'FontSize', 26, 'LineWidth', 1.5);
box on;
outdir='/projects/b1081/NSF_HUBS/images/manuscript/verified/all_tasks_network_bar';
% 
if ~isfolder(outdir), mkdir(outdir); end
saveas(gcf,[outdir '/all_control_tasks.jpg'], 'jpg')


%%
clear all;
load('/projects/b1081/NSF_HUBS/results/task_responses.mat');
outdir='/projects/b1081/NSF_HUBS/images/manuscript/all_tasks_network_bar';

% Standardize network labels in T.network
T.network = replace(T.network, 'DNA', 'DN-A');
T.network = replace(T.network, 'DNB', 'DN-B');
T.network = replace(T.network, 'SAL', 'SAL/PMN');

% Define mapping from short task names to long names
task_map = containers.Map(...
    {'audattn', 'visattn', 'viswm', 'audwm', ...
     'msit', 'vmsit', 'spatialwm', 'verbalwm'}, ...
    {'Auditory Attention', 'Visual Attention', 'Visual N-back', 'Auditory N-back', ...
     'MSIT', 'vMSIT', 'Spatial WM', 'Verbal WM'});

tasks = unique(T.task);
networks = unique(T.network);

results = table();  % initialize
row = 1;

% Step 1: Run one-sided t-tests for each task × network
for i = 1:numel(tasks)
    for j = 1:numel(networks)
        idx = T.task == tasks(i) & T.network == networks(j);
        zvals = T.z(idx);
        [~, p, ~, stats] = ttest(zvals, 0, 'Tail', 'right');
        
        % Compute Cohen's d
        d = mean(zvals) / std(zvals);

        % Store raw results
        results.Task(row)      = tasks(i);
        results.TaskLabel(row) = {task_map(char(tasks(i)))};
        results.Network(row)   = networks(j);
        results.MeanZ(row)     = mean(zvals);
        results.Tstat(row)     = stats.tstat;
        results.Pval(row)      = p;
        results.CohensD(row)   = d;
        row = row + 1;
    end
end

% Step 2: Bonferroni correction within each task
results.Bonf_Pval = nan(height(results), 1);
for i = 1:numel(tasks)
    idx = results.Task == tasks(i);
    corrected = results.Pval(idx) * sum(idx);
    corrected(corrected > 1) = 1;
    results.Bonf_Pval(idx) = corrected;
end

% Step 3: Round all numeric results
results.MeanZ     = round(results.MeanZ, 2);
results.Tstat     = round(results.Tstat, 2);
results.Pval      = round(results.Pval, 4);
results.Bonf_Pval = round(results.Bonf_Pval, 4);
results.CohensD   = round(results.CohensD, 2);


% Step 4: Sort by custom task and network order
task_order = {'viswm', 'audwm', 'spatialwm', 'verbalwm', ...
              'msit', 'vmsit', 'visattn', 'audattn'};
results.Task = categorical(results.Task, task_order, 'Ordinal', true);

network_order = {'FP', 'DAN', 'CO', 'DN-A', 'DN-B', 'LANG', 'SAL/PMN'};
results.Network = categorical(results.Network, network_order, 'Ordinal', true);

results = sortrows(results, {'Task', 'Network'});

% Step 5: Write to file
writetable(results, [outdir '/all_control_tasks_stats.txt'], ...
    'Delimiter', '\t', 'FileType', 'text');


%%
clear all;
load('/projects/b1081/NSF_HUBS/results/task_responses.mat');

[p, tbl, stats] = anovan(T.z, ...
    {T.task, T.network, T.subject}, ...
    'model', 'interaction', ...
    'random', 3, ...          % subject is random
    'varnames', {'Task', 'Network', 'Subject'});
