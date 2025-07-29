README.txt 

Surface Processing (already done)
----------------------------------------
1. generate task datalist for subject (/GrattonLab-General-Repo/SurfacePipeline/gen_datalist.sh, extend_datalist.m)
2. put tasks on the surface (post_fc_processing_batch_GrattonLab_HUBS_task.m)


Task Analysis 
--------------------------------------
1. Generate timing files if not already done (generate_timing_files.m) - uses the behavioral data from each subject

2. kickoff_all_task_steps.sh  
- runs the GLM on a per run basis 
- outputs a dtseries per run, including all contrasts 
 
3. combine runs.m
- selects only the contrasts we care about and creates a dtseries per task per subject (including all runs) 
- includes a single images per run + split half 1 + split 2 + overall
- also outputs a dscalar with just the average map per subject/task

4. combine tasks.m
- combines related tasks into cognitive domains of interest (multiple demand, theory of mind, episodic projection, language)
- outputs a dtseries with single images per run (of all relevant tasks) + split halves + an overall 
- also outputs just a dscalar with the average per sub/domain




