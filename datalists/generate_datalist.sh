#!/bin/bash

header="sub,sess,task,TR,dropFr,topDir,dataFolder,confoundsFolder,FDtype,runs"

subject=$1
type=$2
tr="1.355";
frames="5";
topDir="/projects/b1081/NSF_HUBS/Nifti/derivatives/";
dataFolder="preproc_fmriprep-23.2.0a2_all/";
confoundsFolder="preproc_fmriprep-23.2.0a2_all/";
FDtype="fFD";

output_dir="/projects/b1081/NSF_HUBS/datalists"
output_file="${subject}_${type}_datalist.txt"
# Directory to search for files
directory="/projects/b1081/NSF_HUBS/Nifti/derivatives/preproc_fmriprep-23.2.0a2_all/fmriprep/sub-${subject}/"

mkdir -p "${output_dir}"
echo $header >> "${output_dir}/${output_file}"

# Find all files ending in "*res-2_desc-preproc_bold.nii.gz"
if [[ "$type" == "rest" ]]
then 
    files=$(find "$directory" -type f -name "*res-2_desc-preproc_bold.nii.gz" -name "*rest*")
else 
    files=$(find "$directory" -type f -name "*res-2_desc-preproc_bold.nii.gz" ! -name "*rest*")
fi 


# Loop through each file
for file in $files; do
    # Extract the word after "task-" and before "_"
    task=$(echo "$file" | grep -oP "(?<=task-).*?(?=_)")
    
    # Extract the number after "ses-" and before "task"
    session=$(echo "$file" | grep -oP "(?<=_ses-).*?(?=_task)")

    # Extract the number after "run-" and before "_"
    run=$(echo "$file" | grep -oP "(?<=run-).*?(?=_)")

    if [ -z "$run" ]; then
        run=1
        echo "$subject,$session,$task,$tr,$frames,$topDir,$dataFolder,$confoundsFolder,$FDtype,$run">> "${output_dir}/${output_file}"
    elif [ "$run" -gt 0 ]; then
        echo "$subject,$session,$task,$tr,$frames,$topDir,$dataFolder,$confoundsFolder,$FDtype,$run">> "${output_dir}/${output_file}"
    fi

done

module load matlab/r2021b
matlab -nosplash -nodesktop -singleCompThread -r "fix_datalist('$subject', '$type')"
