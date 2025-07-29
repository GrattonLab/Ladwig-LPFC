#!/bin/tcsh -xef

input_dir=$1
file=$2
output_dir=$3
stim_dir=$4
cue_time=3;
block_time=40;
model_type='BLOCK'

module load afni
mkdir -p $output_dir
cd $output_dir

# do our regression 
3dDeconvolve                                                              \
    -force_TR 1.355                                                       \
    -input ${input_dir}/${file}_LR_surf_subcort_222_32k_fsLR_smooth1_scale.nii          \
    -polort A -float                                                   \
    -num_stimts 9                                                      \
    -stim_times 1 $stim_dir/${file}_cue.1D "${model_type}(${cue_time},1)"        \
    -stim_label 1 cue                                               	                    \
    -stim_times 2 $stim_dir/${file}_femvisact.1D "${model_type}(${block_time},1)"        \
    -stim_label 2 femvisact                                                                 \
    -stim_times 3 $stim_dir/${file}_malvisact.1D "${model_type}(${block_time},1)"        \
    -stim_label 3 malvisact                                               		            \
    -stim_times 4 $stim_dir/${file}_dogaudact.1D "${model_type}(${block_time},1)"        \
    -stim_label 4 dogaudact                                                                 \
    -stim_times 5 $stim_dir/${file}_cataudact.1D "${model_type}(${block_time},1)"        \
    -stim_label 5 cataudact                                                                 \
    -stim_times 6 $stim_dir/${file}_femvispas.1D "${model_type}(${block_time},1)"        \
    -stim_label 6 femvispas                                                                 \
    -stim_times 7 $stim_dir/${file}_malvispas.1D "${model_type}(${block_time},1)"        \
    -stim_label 7 malvispas                                                                 \
    -stim_times 8 $stim_dir/${file}_dogaudpas.1D "${model_type}(${block_time},1)"        \
    -stim_label 8 dogaudpas                                                                 \
    -stim_times 9 $stim_dir/${file}_cataudpas.1D "${model_type}(${block_time},1)"        \
    -stim_label 9 cataudpas                                                                 \
    -gltsym "SYM: femvisact +malvisact -dogaudact -cataudact"                               \
    -glt_label 1 vis_act-aud_act                                                            \
    -gltsym "SYM: femvisact +malvisact -femvispas -malvispas"                               \
    -glt_label 2 vis_act-vis_pas                                                            \
    -gltsym "SYM: dogaudact +cataudact -dogaudpas -cataudpas"                               \
    -glt_label 3 aud_act-aud_pas                                                            \
    -jobs 8                                                                                 \
    -fout -tout                                                \
    -x1D ${output_dir}/${file}_X.xmat.1D                    \
    -xjpeg ${output_dir}/${file}_X.jpg                                     \
    -fitts ${output_dir}/${file}_fitts.nii                     \
    -errts ${output_dir}/${file}_errts.nii                  \
    -bucket ${output_dir}/${file}_stats.nii                     \
    -overwrite    \

df=$(3dAttribute BRICK_STATAUX ${output_dir}/${file}_stats.nii'[35]' | awk '{print $NF}')
echo $df

3dcalc \
    -a ${output_dir}/${file}_stats.nii \
    -expr "fitt_t2z (a,${df})" \
    -prefix ${output_dir}/${file}_zstats.nii \
    -float \
    -overwrite
