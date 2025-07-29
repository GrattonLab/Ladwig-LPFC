clear all;
base_input_dir = '/projects/b1081/NSF_HUBS/Nifti/derivatives/postFCproc_CIFTI_23.2.0/';
base_output_dir = '/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/';
taskname = 'task-rest_desc-mode1000_mean_LR_surf_32k_fsLR.dtseries.nii';

for subj_num = 1:10
    subj_id = sprintf('HUBS%02d', subj_num);
    subj_dir = fullfile(base_input_dir, ['sub-' subj_id]);
    ses_dirs = dir(fullfile(subj_dir, 'ses-*'));
    
    all_data = [];
    n_sess = 0;
    
    for s = 1:length(ses_dirs)
        ses_name = ses_dirs(s).name;
        cifti_path = fullfile(subj_dir, ses_name, 'cifti_timeseries_normalwall', ...
            sprintf('sub-%s_%s_%s', subj_id, ses_name, taskname));
        
        if exist(cifti_path, 'file')
            cifti = ft_read_cifti_mod(cifti_path);
            if isempty(all_data)
                all_data = double(cifti.data);
            else
                all_data = all_data + double(cifti.data);
            end
            n_sess = n_sess + 1;
        else
            warning('Missing file: %s', cifti_path);
        end
    end
    
    if n_sess > 0
        avg_data = all_data / n_sess;
        cifti.data = avg_data;
        
        out_dir = fullfile(base_output_dir, ['sub-' subj_id]);
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        
        out_fname = sprintf('sub-%s_allsess_task-rest_desc-mode1000_mean_LR_surf_32k_fsLR.dtseries.nii', subj_id);
        out_path = fullfile(out_dir, out_fname);
        ft_write_cifti_mod(out_path, cifti);
    else
        warning('No valid sessions found for %s', subj_id);
    end
end

%%
clear all;
input_group_dir = '/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/';
output_group_dir = fullfile(input_group_dir, 'sub-group');
output_group_fname = 'sub-group_allsess_task-rest_desc-mode1000_mean_LR_surf_32k_fsLR.dtseries.nii';

all_group_data = [];
n_subjects = 0;

for subj_num = 1:10
    subj_id = sprintf('HUBS%02d', subj_num);
    subj_file = fullfile(input_group_dir, ['sub-' subj_id], ...
        sprintf('sub-%s_allsess_task-rest_desc-mode1000_mean_LR_surf_32k_fsLR.dtseries.nii', subj_id));
    
    if exist(subj_file, 'file')
        cifti = ft_read_cifti_mod(subj_file);
        if isempty(all_group_data)
            all_group_data = double(cifti.data);
        else
            all_group_data = all_group_data + double(cifti.data);
        end
        n_subjects = n_subjects + 1;
    else
        warning('Missing subject summary: %s', subj_file);
    end
end

if n_subjects > 0
    group_avg = all_group_data / n_subjects;
    cifti.data = group_avg;
    
    if ~exist(output_group_dir, 'dir')
        mkdir(output_group_dir);
    end
    
    out_path = fullfile(output_group_dir, output_group_fname);
    ft_write_cifti_mod(out_path, cifti);
else
    warning('No subject summaries found  no group average created.');
end
%%
make_alt('/projects/b1081/NSF_HUBS/Nifti/')
