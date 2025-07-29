#!/bin/bash  

#SBATCH -A p31306
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem=0
#SBATCH -p normal   
#SBATCH --job-name=kickoff_array_rotations
#SBATCH --array=1-10

cd /projects/b1081/NSF_HUBS/Ladwig_LPFC/border_density

module load matlab/r2022a

SUBJECT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /projects/b1081/NSF_HUBS/Ladwig_LPFC/resources/subjects.txt)

outputdir="/projects/b1081/NSF_HUBS/Nifti/derivatives/RestStats_CIFTI_2320/sub-${SUBJECT}/functional_masks"
inputcifti="${outputdir}/rostral_CO_region.dtseries.nii"
echo "generate_cifti_rotations('${inputcifti}',10000,'${outputdir}',10); exit"

matlab -nosplash -nodesktop -singleCompThread -r "addpath(genpath('/projects/b1081/Scripts/CIFTI_RELATED')); generate_cifti_rotations('${inputcifti}',10000,'${outputdir}',10); exit"
