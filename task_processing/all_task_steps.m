function all_task_steps(subject)

    addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED/'));
    
    datalist = ['/projects/b1081/NSF_HUBS/Ladwig_LPFC/datalists/' subject '_task_datalist.txt']; % change this to datalist for your project
    
    dataInfo_all = readtable(datalist);

    % filter to just ones you want
    %dataInfo = dataInfo_all;
    dataInfo = dataInfo_all(find(strcmp(dataInfo_all.task, 'langlocaud')),:);
    
    numdatas=size(dataInfo.sub,1); %number of datasets to analyses (subs X sessions)
    stim_dir = ['/projects/b1081/NSF_HUBS/Ladwig_LPFC/task_processing/timing/' subject];
   
    for i=1:numdatas
        data_dir=['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskProc_CIFTI_23.2.0/sub-' dataInfo.sub{i} '/ses-' num2str(dataInfo.sess(i)) '/cifti_timeseries_normalwall'];
        scale_dir=['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' dataInfo.sub{i} '/task_glm_outputs/' dataInfo.task{i} '/scale']; 
        stats_dir=['/projects/b1081/NSF_HUBS/Nifti/derivatives/TaskStats_CIFTI_23.2.0/sub-' dataInfo.sub{i} '/task_glm_outputs/' dataInfo.task{i} '/stats']; 
    
        system(['mkdir -p ' scale_dir]);
        system(['mkdir -p ' stats_dir]);
    
        run_nums{i} = str2double(regexp(dataInfo.runs{i},',','split'))'; % get runs, converting to numerical array (other orientiation since that's what's expected
        run_nums{i}=run_nums{i}(1:end-1);
    
        for r = 1:length(run_nums{i})
            if(length(run_nums{i})==1)
                file = sprintf('sub-%s_ses-%d_task-%s_LR_surf_subcort_222_32k_fsLR_smooth1',dataInfo.sub{i},dataInfo.sess(i),dataInfo.task{i});
            else
                file = sprintf('sub-%s_ses-%d_task-%s_run-0%d_LR_surf_subcort_222_32k_fsLR_smooth1',dataInfo.sub{i},dataInfo.sess(i),dataInfo.task{i},run_nums{i}(r));
            end
            disp(file);
    
            % make a fake nifti 
            fake_nifti(data_dir, data_dir, file);
    
            % run task proc 
            task_proc(data_dir, scale_dir, file);
    
            % run task script
            task_glm(scale_dir, stats_dir, file, stim_dir, dataInfo.task{i});
    
            % convert to cifti 
            fake_nifti_to_cifti(stats_dir, stats_dir, file);
        end
    end
end

function fake_nifti(input_dir, output_dir, file)
        workbenchdir = '/projects/b1081/Scripts/workbench_2.1.0/bin_rh_linux64/';
        cifti_in = [input_dir '/' file '.dtseries.nii'];
        nifti_out = [output_dir '/' file '.nii'];
        system(['module load connectome_workbench/1.5.0 && wb_command -cifti-convert -to-nifti ' cifti_in ' ' nifti_out]);
        disp('done fake nifti');
end

function task_proc(input_dir, output_dir, file)
        system(['sh task_proc.sh ' input_dir ' ' file ' ' output_dir]);
end

function task_glm(input_dir, output_dir, file, stim_dir, task)
    taskref=readtable('/projects/b1081/NSF_HUBS/Ladwig_LPFC/task_processing/resources/task_reference.txt');
    file=extractBefore(file, '_LR');

    if(strcmp(task,'epiproj') || strcmp(task,'audviswm') || strcmp(task,'audvisattn'))
        system(['sh task_script_' task '.sh ' input_dir ' ' file ' ' output_dir ' ' stim_dir]);
    else
        condition1 = taskref(strcmp(taskref.task,task),:).condition1{1};
        condition2 = taskref(strcmp(taskref.task,task),:).condition2{1};
        block = taskref(strcmp(taskref.task,task),:).block;
        system(['sh task_script_gen.sh ' input_dir ' ' file ' ' output_dir ' ' stim_dir ' ' task ' ' num2str(block) ' ' condition1 ' ' condition2]);
    end
end

function fake_nifti_to_cifti(input_dir, output_dir, file)
    file=extractBefore(file, '_LR');
    stats = {'stats', 'zstats'};
    for i =1:length(stats)
        stats_file=[file '_' stats{i}];
        workbenchdir = '/projects/b1081/Scripts/workbench2/bin_linux64/';
        cifti_template = '/projects/b1081/Atlases/cifti_template_full.dtseries.nii';
        cifti=ft_read_cifti_mod(cifti_template);
        cifti_out = [output_dir '/' stats_file '.dtseries.nii'];
    
        nifti_in = [input_dir '/' stats_file '.nii'];
        nifti=load_nii(nifti_in);
        num_images=size(nifti.img, length(size(nifti.img)));
        cifti.time=1:num_images;
        cifti.hdr.dim(6)=num_images;
        cifti.data = cifti.data(:,1:num_images);
    
        ft_write_cifti_mod([output_dir '/template_cifti.dtseries.nii'], cifti);
        cifti_template_new=[output_dir '/template_cifti.dtseries.nii'];
        system(['module load connectome_workbench/1.5.0 && wb_command -cifti-convert -from-nifti ' nifti_in ' ' cifti_template_new ' ' cifti_out]);
        delete([output_dir '/template_cifti.dtseries.nii']);
        disp('done fake nifti to cifti');
    end
end
