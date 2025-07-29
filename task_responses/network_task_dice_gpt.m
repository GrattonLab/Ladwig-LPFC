%% Setup
clear; clc;

subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', ...
            'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
tasks = {'epiproj', 'tom', 'lang'};
network_labels = [1, 2, 10]; % DMN-B, DMN-A, LANG
pct_thresh = [75 80 85];
nS = numel(subjects); nT = numel(tasks); nP = numel(pct_thresh);

% Load group-level networks
group_networks_lpfc = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/lynch_group_prior20_lpfc.dlabel.nii');
group_networks = ft_read_cifti_mod('/projects/b1081/NSF_HUBS/resources/lynch_group_prior20.dlabel.nii');
group_data = group_networks.data;

% Masks
lpfc_mask = group_networks_lpfc.data(1:59412) > 0;
group_net_lpfc = group_data(lpfc_mask);
group_net_nonlpfc = group_data(~lpfc_mask);

%% Preload data
disp('Preloading subject data...');
task_data_all = cell(nS, nT);
ind_networks_all = cell(nS, 1);

for s = 1:nS
    subj = subjects{s};

    % Load all task data
    for t = 1:nT
        fname = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0_GSR/sub-' subj ...
                 '/domain_summary/sub-' subj '_' tasks{t} '_zstats_mean.dscalar.nii'];
        task_cifti = ft_read_cifti_mod(fname);
        task_data_all{s,t} = task_cifti.data(1:59412, end); % Use last column
    end

    % Load resting-state networks
    network_file = ['/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-' subj ...
                    '/pfm/lynch/Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii'];
    ind_cifti = ft_read_cifti_mod(network_file);
    ind_networks_all{s} = ind_cifti.data;
end

%% Main Analysis Loop
disp('Beginning main loop...');
for p = 1:nP
    for s = 1:nS
        for t = 1:nT
            task_data = task_data_all{s,t};
            task_lpfc = task_data(lpfc_mask);
            task_nonlpfc = task_data(~lpfc_mask);

            % Thresholding
            z_lpfc = prctile(task_lpfc, pct_thresh(p));
            z_nonlpfc = prctile(task_nonlpfc, pct_thresh(p));
            roi_lpfc = task_lpfc > z_lpfc;
            roi_nonlpfc = task_nonlpfc > z_nonlpfc;

            % Individual network overlaps
            for s2 = 1:nS
                subj_net = ind_networks_all{s2};
                net_lpfc = subj_net(lpfc_mask);
                net_nonlpfc = subj_net(~lpfc_mask);

                overlap_dmna_ind(p,s,t,s2) = dice(roi_lpfc, net_lpfc == 2);
                overlap_dmnb_ind(p,s,t,s2) = dice(roi_lpfc, net_lpfc == 1);
                overlap_lang_ind(p,s,t,s2) = dice(roi_lpfc, net_lpfc == 10);

                overlap_dmna_ind_nonlpfc(p,s,t,s2) = dice(roi_nonlpfc, net_nonlpfc == 2);
                overlap_dmnb_ind_nonlpfc(p,s,t,s2) = dice(roi_nonlpfc, net_nonlpfc == 1);
                overlap_lang_ind_nonlpfc(p,s,t,s2) = dice(roi_nonlpfc, net_nonlpfc == 10);
            end

            % Group overlaps
            overlap_dmna_group(p,s,t) = dice(roi_lpfc, group_net_lpfc == 2);
            overlap_dmnb_group(p,s,t) = dice(roi_lpfc, group_net_lpfc == 1);
            overlap_lang_group(p,s,t) = dice(roi_lpfc, group_net_lpfc == 10);

            overlap_dmna_group_nonlpfc(p,s,t) = dice(roi_nonlpfc, group_net_nonlpfc == 2);
            overlap_dmnb_group_nonlpfc(p,s,t) = dice(roi_nonlpfc, group_net_nonlpfc == 1);
            overlap_lang_group_nonlpfc(p,s,t) = dice(roi_nonlpfc, group_net_nonlpfc == 10);
        end
    end
end

%% Compute Mean Overlaps by Task

% Manual mapping: task 1 = dmna (2), task 2 = dmnb (1), task 3 = lang (10)
mean_overlap_group = {
    mean(overlap_dmna_group(:,:,1), 1);
    mean(overlap_dmnb_group(:,:,2), 1);
    mean(overlap_lang_group(:,:,3), 1);
};

mean_overlap_group_nonlpfc = {
    mean(overlap_dmna_group_nonlpfc(:,:,1), 1);
    mean(overlap_dmnb_group_nonlpfc(:,:,2), 1);
    mean(overlap_lang_group_nonlpfc(:,:,3), 1);
};

mean_overlap = {
    squeeze(mean(overlap_dmna_ind(:,:,1,:), 1));
    squeeze(mean(overlap_dmnb_ind(:,:,2,:), 1));
    squeeze(mean(overlap_lang_ind(:,:,3,:), 1));
};

mean_overlap_nonlpfc = {
    squeeze(mean(overlap_dmna_ind_nonlpfc(:,:,1,:), 1));
    squeeze(mean(overlap_dmnb_ind_nonlpfc(:,:,2,:), 1));
    squeeze(mean(overlap_lang_ind_nonlpfc(:,:,3,:), 1));
};

%%
%% Per-Task Plots and Stats
colors = [255 255 218; 234 51 35; 64 153 153]/255; % DMN-A, DMN-B, LANG
outdir = '/projects/b1081/NSF_HUBS/images/manuscript/verified/network_task_overlap';
if ~isfolder(outdir), mkdir(outdir); end

for t = 1:nT
    task = tasks{t};
    
    % Self vs others
    self = diag(mean_overlap{t});
    nonself = (sum(mean_overlap{t}, 2) - self) / (nS - 1);

    self_nonlpfc = diag(mean_overlap_nonlpfc{t});
    nonself_nonlpfc = (sum(mean_overlap_nonlpfc{t}, 2) - self_nonlpfc) / (nS - 1);

    % Group overlap means
    mean_self(t) = mean(self);
    stderr_self(t) = std(self) / sqrt(nS);
    mean_non(t) = mean(nonself);
    stderr_non(t) = std(nonself) / sqrt(nS);

    mean_group(t) = mean(mean_overlap_group{t});
    stderr_group(t) = std(mean_overlap_group{t}) / sqrt(nS);

    mean_self_nonlpfc(t) = mean(self_nonlpfc);
    stderr_self_nonlpfc(t) = std(self_nonlpfc) / sqrt(nS);
    mean_non_nonlpfc(t) = mean(nonself_nonlpfc);
    stderr_non_nonlpfc(t) = std(nonself_nonlpfc) / sqrt(nS);
    mean_group_nonlpfc(t) = mean(mean_overlap_group_nonlpfc{t});
    stderr_group_nonlpfc(t) = std(mean_overlap_group_nonlpfc{t}) / sqrt(nS);

    % Improvement metric: self - others
    improvement{t} = arrayfun(@(col) mean(mean_overlap{t}(col, col) - mean_overlap{t}(setdiff(1:nS, col), col)), 1:nS);
    improvement_nonlpfc{t} = arrayfun(@(col) mean(mean_overlap_nonlpfc{t}(col, col) - mean_overlap_nonlpfc{t}(setdiff(1:nS, col), col)), 1:nS);

    % T-tests
    [~, p_vs_group(t), ~, stats_vs_group(t)] = ttest(self, mean_overlap_group{t}', 'Tail', 'right');
    [~, p_vs_others(t), ~, stats_vs_others(t)] = ttest(self, nonself, 'Tail', 'right');
    [~, p_improvement_lpfc_vs_nonlpfc(t), ~, stats_improvement(t)] = ttest(improvement{t}, improvement_nonlpfc{t}, 'Tail', 'right');

    % Table of stats
    statsT = table( ...
        ["Group"; "Others"; "LPFC vs NonLPFC"], ...
        [p_vs_group(t); p_vs_others(t); p_improvement_lpfc_vs_nonlpfc(t)], ...
        [p_vs_group(t)*2; p_vs_others(t)*2; p_improvement_lpfc_vs_nonlpfc(t)], ...
        [stats_vs_group(t).tstat; stats_vs_others(t).tstat; stats_improvement(t).tstat], ...
        [stats_vs_group(t).tstat/sqrt(stats_vs_group(t).df); ...
         stats_vs_others(t).tstat/sqrt(stats_vs_others(t).df); ...
         stats_improvement(t).tstat/sqrt(stats_improvement(t).df)], ...
        [stats_vs_group(t).df; stats_vs_others(t).df; stats_improvement(t).df], ...
        'VariableNames', {'Comparison', 'PValue', 'CorrectedPValue', 'TStat', 'CohensD', 'df'});
    
    % Round numeric columns
    numericCols = varfun(@isnumeric, statsT, 'OutputFormat', 'uniform');
    statsT(:, numericCols) = varfun(@(x) round(x, 3, 'significant'), statsT(:, numericCols));

    % Plot
    fig = figure('Position', [100, 100, 800, 600]);
    b = bar([mean_self(t), mean_group(t), mean_non(t)], ...
            'FaceColor', colors(t,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    errorbar([mean_self(t), mean_group(t), mean_non(t)], ...
             [stderr_self(t), stderr_group(t), stderr_non(t)], ...
             'k', 'LineStyle', 'none', 'LineWidth', 2);
    ylim([0 mean_self(t) + stderr_self(t) + 0.2 * mean_self(t)]);

    ylabel('Dice Overlap');
    xticklabels({'Individual', 'Group', 'Other'});
    set(gca, 'FontSize', 36);

    % Save outputs
    saveas(fig, fullfile(outdir, [task '.jpg']));
    writetable(statsT, fullfile(outdir, [task '_stats.txt']), 'Delimiter', '\t');
end

%% the overall stuff
overall_mean_improvement = mean(cell2mat(improvement(:)),1);
mean_improvement = mean(cell2mat(improvement(:)),2);
stderr_improvement = std(cell2mat(improvement(:)),[],2)/sqrt(10);
mean_overall_mean_improvement = mean(overall_mean_improvement);
sem_overall_mean_improvement = std(overall_mean_improvement)/sqrt(size(overall_mean_improvement,2));

overall_mean_improvement_nonlpfc = mean(cell2mat(improvement_nonlpfc(:)));
mean_improvement_nonlpfc = mean(cell2mat(improvement_nonlpfc(:)),2);
stderr_improvement_nonlpfc = std(cell2mat(improvement_nonlpfc(:)),[],2)/sqrt(10);
mean_overall_mean_improvement_nonlpfc = mean(overall_mean_improvement_nonlpfc);
sem_overall_mean_improvement_nonlpfc = std(overall_mean_improvement_nonlpfc)/sqrt(size(overall_mean_improvement_nonlpfc,2));

[~, p_overall_lpfc_vs_non, ~, stats_overall] = ttest(overall_mean_improvement, overall_mean_improvement_nonlpfc,'Tail','right');

fig = figure('Position', [100, 100, 800, 600]);
hb=bar([mean_overall_mean_improvement mean_overall_mean_improvement_nonlpfc mean_improvement(1) mean_improvement_nonlpfc(1) mean_improvement(2) mean_improvement_nonlpfc(2) mean_improvement(3) mean_improvement_nonlpfc(3)],'EdgeColor', 'k', 'LineWidth', 1.5);
hold on
err = errorbar([mean_overall_mean_improvement mean_overall_mean_improvement_nonlpfc mean_improvement(1) mean_improvement_nonlpfc(1) mean_improvement(2) mean_improvement_nonlpfc(2) mean_improvement(3) mean_improvement_nonlpfc(3)], [sem_overall_mean_improvement sem_overall_mean_improvement_nonlpfc stderr_improvement(1) stderr_improvement_nonlpfc(1) stderr_improvement(2) stderr_improvement_nonlpfc(2) stderr_improvement(3) stderr_improvement_nonlpfc(3)], 'vertical','LineStyle', 'none', 'Color', 'k', 'LineWidth', 2);
colors =[
    128 128 128;
    128 128 128;
    255 255 218;% DMN-A
    255 255 218;% DMN-A
    234 51 35; % DMN-B
    234 51 35; % DMN-B
    64 153 153; %LANG
    64 153 153; %LANG
    ];
colors = colors/255;
ylabel('Dice: Individual - Other');
xticklabels({'LPFC', 'Non-LPFC', 'LPFC','Non-LPFC', 'LPFC', 'Non-LPFC', 'LPFC', 'Non-LPFC'})
for k = 1:size(colors,1)
    hb.FaceColor = 'flat';
    hb.CData(k,:) = colors(k,:);
end

statsT = table( ...
strings(0,1), ... % Comparison
zeros(0,1), ...   % PValue
zeros(0,1), ...   % CorrectedPValue
zeros(0,1), ...   % TStat
zeros(0,1), ...   % CohensD
zeros(0,1), ... % df 
'VariableNames', {'Comparison', 'PValue', 'CorrectedPValue', 'TStat', 'CohensD', 'df'});

stats_row_group = {'Overall LPFC vs Not',p_overall_lpfc_vs_non,p_overall_lpfc_vs_non,stats_overall.tstat,stats_overall.tstat/sqrt(stats_overall.df), stats_overall.df};
statsT = [statsT; stats_row_group];
numVars = varfun(@isnumeric, statsT, 'OutputFormat', 'uniform');

% Apply rounding only to numeric variables
statsT(:, numVars) = varfun(@(x) round(x, 3, 'significant'), statsT(:, numVars));    
writetable(statsT, [outdir '/overall_stats.txt'], 'Delimiter', '\t');
set(gca, 'FontSize', 24);

saveas(gcf,[outdir '/lpfc_vs_nonlpfc_dice_improvement.jpg'], 'jpg')
