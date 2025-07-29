clear all;
subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'}; 
%subjects = {'HUBS04'};
for s = 1:length(subjects)
    subject = subjects{s};
    disp(subject);
    datalist = ['/projects/b1081/NSF_HUBS/Ladwig_LPFC/datalists/' subject '_task_datalist.txt']; % change this to datalist for your project
    dataInfo = readtable(datalist);
    dataInfo = dataInfo(strcmpi(dataInfo.task, 'epiproj'), :);
    %dataInfo = dataInfo(~strcmp(dataInfo.task, 'arithmetic'), :);

    numdatas=size(dataInfo.sub,1); %number of datasets to analyses (subs X sessions)
    outdir = ['/projects/b1081/NSF_HUBS/Ladwig_LPFC/task_processing/timing/' subject];
    system(['mkdir -p ' outdir]);
    
    for i=1:numdatas
        %disp(i);
        run_nums = str2double(regexp(dataInfo.runs{i},',','split'))';
     % get runs, converting to numerical array (other orientiation since that's what's expected
        overall_run_nums = str2double(regexp(dataInfo.overall_runs{i},',','split'))'; % get runs, converting to numerical array (other orientiation since that's what's expected   
        run_nums=run_nums(1:end-1);
        overall_run_nums=overall_run_nums(1:end-1);
        
        for r = 1:length(run_nums)
            if(length(run_nums)==1)
                outfile=[outdir '/sub-' dataInfo.sub{i}, '_ses-' num2str(dataInfo.sess(i)) '_task-' dataInfo.task{i}];
            else
                outfile=[outdir '/sub-' dataInfo.sub{i}, '_ses-' num2str(dataInfo.sess(i)) '_task-' dataInfo.task{i} '_run-0' num2str(run_nums(r))];
            end
            if(strcmp(dataInfo.task{i},'epiproj') || strcmp(dataInfo.task{i},'audviswm') || strcmp(dataInfo.task{i},'audvisattn') || strcmp(dataInfo.task{i},'tomfalse') || strcmp(dataInfo.task{i},'tompain') || strcmp(dataInfo.task{i}, 'langlocaud') || strcmp(dataInfo.task{i}, 'langlocvis'))
                func = str2func(dataInfo.task{i});
                func(subject, overall_run_nums(r), outfile)       
            else 
                get_md_timing(subject, dataInfo.task{i}, overall_run_nums(r), outfile);
            end
        end     
    end
end

function tompain(subject, run, outfile)
    % Load the relevant .mat file and extract the table
    T = load(['/projects/b1081/NSF_HUBS/task_behavior/' subject '/tompain/run' num2str(run) '.mat']);
    T = T.allData;

    % Extract and clean labels and onset times
    labels = replace(lower(string(T(:,1))), '_', '');
    onset = cell2mat(T(:,8));

    % Remove 'fix' rows
    valid = labels ~= "fix";
    labels = labels(valid);
    onset = onset(valid);

    % Get unique condition names
    conds = unique(labels);

    % Map labels for file output
    labelMap = containers.Map({'Emo', 'Phys'}, {'emo', 'phy'});

    % Write each conditions onsets to its own .1D file
    for i = 1:numel(conds)
        idx = labels == conds(i);
        cleanLabel = conds(i);

        if isKey(labelMap, cleanLabel)
            fileLabel = labelMap(cleanLabel);
        else
            fileLabel = cleanLabel;  % fallback: use as-is
        end

        dlmwrite([outfile '_' char(fileLabel) '.1D'], round(onset(idx), 1)', 'delimiter', ' ');
    end
end

function tomfalse(subject, run, outfile)
    % Load the relevant .mat file and extract the table
    T = load(['/projects/b1081/NSF_HUBS/task_behavior/' subject '/tomfalse/run' num2str(run) '.mat']);
    T = T.allData;

    % Extract and clean labels and onset times
    labels = replace(lower(string(T(:,1))), '_', '');
    onset = cell2mat(T(:,10));

    % Remove 'fix' rows
    valid = labels ~= "fix";
    labels = labels(valid);
    onset = onset(valid);

    % Get unique condition names
    conds = unique(labels);

    % Map labels for file output
    labelMap = containers.Map({'belief', 'photo'}, {'bel', 'pic'});

    % Write each conditions onsets to its own .1D file
    for i = 1:numel(conds)
        idx = labels == conds(i);
        cleanLabel = conds(i);

        if isKey(labelMap, cleanLabel)
            fileLabel = labelMap(cleanLabel);
        else
            fileLabel = cleanLabel;  % fallback: use as-is
        end

        dlmwrite([outfile '_' char(fileLabel) '.1D'], round(onset(idx), 1)', 'delimiter', ' ');
    end
end

%% epiproj new 
function epiproj(subject, run, outfile)
    T = load(['/projects/b1081/NSF_HUBS/task_behavior/' subject '/epiproj/run' num2str(run) '.mat']);
    T = T.allData;
    labels = replace(lower(string(T(:,1))), '_', '');
    onset = cell2mat(T(:,8));

    % Mask out 'fix' rows
    valid = labels ~= "fix";
    labels = labels(valid);
    onset = onset(valid)+5;

    % Get unique condition names
    conds = unique(labels);

    % Write each condition's onsets to file
    for i = 1:numel(conds)
        idx = labels == conds(i);
        dlmwrite([outfile '_' char(conds(i)) '.1D'], round(onset(idx), 1)', 'delimiter', ' ');
    end
end    

%% audviswm
function audviswm(subject, run, outfile)
%AUDVISWM Extracts block onset times (first appearance only) by condition
%
% Inputs:
%   matfile - path to .mat file with 266x10 cell array T
%   outfile - base name for output files (e.g., 'sub-HUBS02_audviswm')
    T = load(['/projects/b1081/NSF_HUBS/task_behavior/' subject '/audviswm/run' num2str(run) '.mat']);

    % Load data
    T = T.allData;

    % Extract relevant columns
    labels = string(T(:,1));
    onset  = cell2mat(T(:,9));

    label_map = containers.Map( ...
        {'Vis Female Control', 'Vis Fem', ...
         'Vis Male Control',   'Vis Male', ...
         'Aud Dog Control',    'Aud Dog', ...
         'Aud Cat Control',    'Aud Cat', ...
         'Cue'}, ...
        {'femvispas', 'femvisact', ...
         'malvispas', 'malvisact', ...
         'dogaudpas', 'dogaudact', ...
         'cataudpas', 'cataudact', ...
         'cue'} ...
    );

    % Loop over mapped labels
    for k = keys(label_map)
        orig_label = k{1};
        short_label = label_map(orig_label);

        % Get all matching rows
        idx = find(labels == orig_label);

        % Special case: Cue  keep all
        if strcmpi(short_label, 'cue')
            onsets_out = round(onset(idx), 1);
        else
            if isempty(idx)
                continue;  % skip if not found
            end
            onsets_out = round(onset(idx(1)), 1);  % only first
        end

        % Write to file
        dlmwrite([outfile '_' short_label '.1D'], onsets_out', 'delimiter', ' ');
    end
end

%% audvisattn new 
function audvisattn(subject, run, outfile)
%AUDVISATTN Extract block onsets for each condition, using 'cue' as block marker
% - All conditions: grab rows that follow a 'cue'
% - 'cue': grab all cue rows

    % Load .mat file
    X = load(['/projects/b1081/NSF_HUBS/task_behavior/' subject '/audvisattn/run' num2str(run) '.mat']);
    T = X.allData;

    % Extract condition labels and onset times
    raw_labels = string(T(:,1));
    onset = cell2mat(T(:,11));
    onset_passive = cell2mat(T(:,12));
    norm_labels = replace(lower(raw_labels), ' ', '');

    % Identify all unique conditions
    all_conditions = unique(norm_labels);
    all_conditions(all_conditions == "fix") = [];

    % Find all cue rows
    cue_idx = find(norm_labels == 'cue');

    for i = 1:numel(all_conditions)
        cond = char(all_conditions(i));

        if strcmp(cond, 'cue')
            % All cue rows
            idx = cue_idx;
        else
            % Rows following a cue AND matching this condition
            idx = cue_idx + 1;
            idx = idx(idx <= numel(norm_labels) & norm_labels(idx) == cond);
        end

        % Write file if any onsets exist
        if ~isempty(idx)
            if(strcmp(cond,'passive')) % we changed how we were tracking this condition after HUBS01 so need to do this
                onsets_out = round(onset_passive(idx-1),1);
                dlmwrite([outfile '_' cond '.1D'], onsets_out', 'delimiter', ' ');
            else
                onsets_out = round(onset(idx), 1);
                dlmwrite([outfile '_' cond '.1D'], onsets_out', 'delimiter', ' ');
            end
        end
    end
end

%% langlocaud
function langlocaud(subject, run, outfile)
%LANGLOCAUD Extract onset times for 'int' and 'degr' stimuli from timing file
%
% Inputs:
%   table_file - path to timing .txt or .csv file
%   outfile    - base name for output files (e.g., 'sub-HUBS02_langlocaud')

    % Read table
    T = readtable(['/projects/b1081/NSF_HUBS/task_behavior/' subject '/langlocaud/run' num2str(run) '.txt']);

    % Extract timing and stimulus code columns
    onset = T{:,3};     % Block onset time
    stim_code = T{:,4}; % Stimulus code (0=rest, 1=int, 2=degr)

    % Extract and round onset times
    onset_int  = round(onset(stim_code == 1), 0);
    onset_degr = round(onset(stim_code == 2), 0);

    % Write to .1D files
    dlmwrite([outfile '_int.1D'],  onset_int',  'delimiter', ' ');
    dlmwrite([outfile '_deg.1D'], onset_degr', 'delimiter', ' ');
end
%% langlocvis 
function langlocvis(subject, run, outfile)
%LANGLOCVIS Extract onset times of every 3-trial set for 'S' and 'N' stimuli
% Uses onset from column 7, offset so first event starts at 14.0s

    X = ['/projects/b1081/NSF_HUBS/task_behavior/' subject '/langlocvis/run' num2str(run) '.mat'];
    % Load .mat file
    T = load(X).allData;

    % Extract columns
    stim_type = string(T(:,3));
    list_id   = cell2mat(T(:,2));
    onset     = cell2mat(T(:,7));

    % Offset so first onset = 14.0
    onset = onset + (14.0 - onset(1));

    % Get onsets for S where list_id is 1,4,7,...
    is_S = stim_type == "S";
    onset_sen = round(onset(is_S & mod(list_id - 1, 3) == 0), 1);

    % Get onsets for N where list_id is 1,4,7,...
    is_N = stim_type == "N";
    onset_non = round(onset(is_N & mod(list_id - 1, 3) == 0), 1);

    % Write to .1D files
    dlmwrite([outfile '_sen.1D'], onset_sen', 'delimiter', ' ');
    dlmwrite([outfile '_non.1D'], onset_non', 'delimiter', ' ');
end

%%
function get_md_timing(subject, task, run, outfile)
%GET_TASK_TIMING Extracts and writes block start times for 'easy' and 'hard'
% Includes a block if the label changed OR it's been trials_per_block rows
%
% Inputs:
%   matfile   - path to .mat file with cell array T
%   task_name - string name of the task (e.g., 'spatialwm')
%   outfile   - base name for output files (e.g., 'sub-HUBS01_spatialwm')

    Tinfo=readtable('/projects/b1081/NSF_HUBS/Ladwig_LPFC/task_processing/resources/task_reference.txt');
    blocktype_col = Tinfo.blocktype_col(strcmpi(Tinfo.task,task));
    timing_col = Tinfo.timing_col(strcmpi(Tinfo.task, task));
    trials_per_block = Tinfo.trials_per_block(strcmpi(Tinfo.task,task));

    X = ['/projects/b1081/NSF_HUBS/task_behavior/' subject '/' task '/run' num2str(run) '.mat'];
    % Load .mat file
    S = load(X);
    T = S.allData;

    rawVals = T(:, blocktype_col);  % Get column 2
    vals = cellfun(@(x) string(x), rawVals, 'UniformOutput', false);  % Try to convert each to string
    vals(cellfun(@isempty, rawVals)) = {""};  % Replace empty entries
    block_type = lower(string(vals));  % Final string array    

    %block_type = lower(string(T(:,blocktype_col)));
    start_time = cell2mat(T(:, timing_col));
    n = size(T, 1);

    easy_rows = [];
    hard_rows = [];
    last_easy_idx = -Inf;
    last_hard_idx = -Inf;

    for i = 1:n
        curr = block_type(i);
        prev = ""; if i > 1, prev = block_type(i-1); end

        % Easy logic
        if curr == "easy" || curr == "1"
            if ~strcmp(curr, prev) || (i - last_easy_idx) >= trials_per_block
                easy_rows(end+1) = i;
                last_easy_idx = i;
            end
        end

        % Hard logic
        if curr == "hard" || curr == "2"
            if ~strcmp(curr, prev) || (i - last_hard_idx) >= trials_per_block
                hard_rows(end+1) = i;
                last_hard_idx = i;
            end
        end
    end

    % Round and save
    easy_array = round(start_time(easy_rows), 1);
    hard_array = round(start_time(hard_rows), 1);

    dlmwrite([outfile '_easy.1D'], easy_array', 'delimiter', ' ');
    dlmwrite([outfile '_hard.1D'], hard_array', 'delimiter', ' ');
end
