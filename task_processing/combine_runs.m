clear all;

subjects = {'HUBS01', 'HUBS02', 'HUBS03', 'HUBS04', 'HUBS05', 'HUBS06', 'HUBS07', 'HUBS08', 'HUBS09', 'HUBS10'};
%subjects = {'HUBS02'};
dscalar_template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_full.dscalar.nii');
dtseries_template = ft_read_cifti_mod('/projects/b1081/Atlases/cifti_template_full.dtseries.nii');
%subjects = {'HUBS04'};
for s = 1:length(subjects)
sub = subjects{s};
disp(sub);
contrasts=readtable("/projects/b1081/NSF_HUBS/Ladwig_LPFC/task_processing/resources/all_contrasts.csv", "Delimiter", ",");
%contrasts = contrasts(find(strcmp(contrasts.task, 'langlocaud')),:);

dataInfo=readtable(['/projects/b1081/NSF_HUBS/Ladwig_LPFC/datalists/' sub '_task_datalist.txt']);

dataInfoExt = innerjoin(dataInfo,contrasts);

input_dir=['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' sub '/task_glm_outputs'];
output_dir=['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' sub '/task_summaries'];
system(['mkdir -p ' output_dir]);
badruns = readtable('/projects/b1081/NSF_HUBS/Ladwig_LPFC/datalists/all_bad_runs.txt');

for i = 1:size(contrasts,1)
    contrast=contrasts.name{i};
    disp(contrast);
    if(contains(contrast,'epiproj'))
        contrast_name = 'epiproj';
    else
        contrast_name = contrast;
    end
    all_ses=find(strcmp(dataInfoExt.name, contrasts.name{i}));
    run=1;
    z_data =[];
    t_data = [];
    b_data = [];
    for j = 1:size(all_ses)
        row = all_ses(j);
        idx=dataInfoExt.idx(all_ses(j));
        run_nums = str2double(regexp(dataInfoExt.runs{row},',','split'))'; % get runs, converting to numerical array (other orientiation since that's what's expected
        run_nums=run_nums(1:end-1);
        overall_run_nums = str2double(regexp(dataInfoExt.overall_runs{row},',','split'))'; % get runs, converting to numerical array (other orientiation since that's what's expected
        overall_run_nums=overall_run_nums(1:end-1);
        for r = 1:length(run_nums)
            if(length(run_nums)==1)   
                file = sprintf('sub-%s_ses-%d_task-%s',dataInfoExt.sub{row},dataInfoExt.sess(row),dataInfoExt.task{row});
            else
                file = sprintf('sub-%s_ses-%d_task-%s_run-0%d',dataInfoExt.sub{row},dataInfoExt.sess(row),dataInfoExt.task{row},run_nums(r));
            end

            % check if bad run
            bad_run = any(string(badruns.subject_motion) == sub & string(badruns.contrast_name)==contrast_name & badruns.overall_run == overall_run_nums(r));
            if ~bad_run
                z_cifti = ft_read_cifti_mod([input_dir '/' dataInfoExt.task{row} '/stats/' file '_zstats.dtseries.nii']);
%                 t_cifti = ft_read_cifti_mod([input_dir '/' dataInfoExt.task{row} '/stats/' file '_stats.dtseries.nii']);
%                 b_cifti = ft_read_cifti_mod([input_dir '/' dataInfoExt.task{row} '/stats/' file '_stats.dtseries.nii']);
                z_data(:,run) = z_cifti.data(:,idx);
%                 t_data(:,run) = t_cifti.data(:,idx);
%                 b_data(:,run) = b_cifti.data(:,idx-1);

                run=run+1;
            else 
                disp(['bad run ' sub ' ' contrast ' ' num2str(run)]);
            end

        end
    end
    num_runs=run-1;  
    max_split_1 = [1 4 6 7 9 12 14 15];
    max_split_2 = [2 3 5 8 10 11 13 16];

    if(ismember(num_runs,max_split_1))
        split_1=max_split_1(1:ceil(num_runs/2));
        split_2=max_split_2(1:floor(num_runs/2));
    else
        split_1=max_split_1(1:floor(num_runs/2));
        split_2=max_split_2(1:ceil(num_runs/2));  
    end
    z_data(:,run) = mean(z_data(:,split_1),2);
    z_data(:,run+1) = mean(z_data(:,split_2),2);
    z_data(:,run+2) = mean(z_data(:,[1:num_runs]),2);

%     t_data(:,run) = mean(t_data(:,split_1),2);
%     t_data(:,run+1) = mean(t_data(:,split_2),2);
%     t_data(:,run+2) = mean(t_data(:,[1:num_runs]),2);
% 
%     b_data(:,run) = mean(b_data(:,split_1),2);
%     b_data(:,run+1) = mean(b_data(:,split_2),2);
%     b_data(:,run+2) = mean(b_data(:,[1:num_runs]),2);    

    dtseries_template.data=z_data;
    dtseries_template.time=1:size(z_data,2);
    dtseries_template.hdr.dim(6)=size(z_data,2);
    ft_write_cifti_mod([output_dir '/sub-' sub '_' contrast '_zstats.dtseries.nii'],dtseries_template);
    dscalar_template.data = z_data(:,size(z_data,2));
    ft_write_cifti_mod([output_dir '/sub-' sub '_' contrast '_zstats_mean.dscalar.nii'],dscalar_template);
    if ~exist([output_dir '/alt'], 'dir'), mkdir([output_dir '/alt']); end
    make_alt([output_dir '/sub-' sub '_' contrast '_zstats_mean.dscalar.nii'],[output_dir '/alt/sub-' sub '_' contrast '_zstats_mean_alt.dscalar.nii'])

%     dtseries_template.data=t_data;
%     dtseries_template.time=1:size(t_data,2);
%     dtseries_template.hdr.dim(6)=size(t_data,2);
%     ft_write_cifti_mod([output_dir '/sub-' sub '_' contrast '_tstats.dtseries.nii'],dtseries_template);
%     dscalar_template.data = t_data(:,size(t_data,2));
%     ft_write_cifti_mod([output_dir '/sub-' sub '_' contrast '_tstats_mean.dscalar.nii'],dscalar_template);
%     make_alt([output_dir '/sub-' sub '_' contrast '_tstats_mean.dscalar.nii'],[output_dir '/sub-' sub '_' contrast '_tstats_mean_alt.dscalar.nii'])
%     
%     dtseries_template.data=b_data;
%     dtseries_template.time=1:size(b_data,2);
%     dtseries_template.hdr.dim(6)=size(b_data,2);
%     ft_write_cifti_mod([output_dir '/sub-' sub '_' contrast '_betas.dtseries.nii'],dtseries_template);
%     dscalar_template.data = b_data(:,size(b_data,2));
%     ft_write_cifti_mod([output_dir '/sub-' sub '_' contrast '_betas_mean.dscalar.nii'],dscalar_template);
%     make_alt([output_dir '/sub-' sub '_' contrast '_betas_mean.dscalar.nii'],[output_dir '/sub-' sub '_' contrast '_betas_mean_alt.dscalar.nii'])

end
end
