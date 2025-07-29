#!/bin/bash  

#SBATCH -A p31306
#SBATCH -n 30
#SBATCH -N 1  ## number of nodes
#SBATCH -t 20:00:00
#SBATCH --mem=100G
#SBATCH -p normal   
#SBATCH --job-name=kickoff_array_task

## set your working directory  
cd /projects/b1081/NSF_HUBS/Ladwig_LPFC/task_processing/

subject=$1

module load matlab/r2022a

matlab -nosplash -nodesktop -singleCompThread -r "all_task_steps('${subject}'); exit"
