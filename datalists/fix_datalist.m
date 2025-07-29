function fix_datalist(subject, type)
    disp(subject);
    path='/projects/b1081/NSF_HUBS/datalists/';
    datalist = [path subject '_' type '_datalist.txt']; % change this to datalist for your project
    
    dataInfo = readtable(datalist);
    dataInfo = sortrows(dataInfo, dataInfo.sess);
    tasks = unique(dataInfo.task);
    sessions = unique(dataInfo.sess);
    
    %%
    to_delete = [];
    for i=1:length(tasks)
        run_counter=1;
        for j=1:length(sessions)
            overall_runs=[];
            run_idx=find(strcmp(dataInfo.task, tasks{i}) & dataInfo.sess==sessions(j));
            if(length(run_idx)==0)
                continue;
            end
    
            run_nums=[];
            for k = 1:length(run_idx)
                run_nums= [run_nums num2str(k) ','];
                overall_runs= [overall_runs num2str(run_counter) ','];
                run_counter=run_counter+1;
            end
            to_delete = [to_delete; run_idx(1:length(run_idx)-1)];
            dataInfo.new_runs{run_idx(length(run_idx))} = run_nums;
            dataInfo.overall_runs{run_idx(length(run_idx))} = overall_runs;
            %dataInfo.new_runs{run_idx(length(run_idx))} = run_nums(1:end-1);
            %dataInfo.overall_runs{run_idx(length(run_idx))} = overall_runs(1:end-1);
        end
    end
    dataInfo(to_delete,:) = [];
    dataInfo.runs=dataInfo.new_runs;
    dataInfo.new_runs=[];
    system(['rm ' datalist]);
    writetable(dataInfo,[path subject '_' type '_datalist.txt'],'WriteRowNames',true);
    exit;
end