#!/bin/tcsh -xef

input_dir=$1
file=$2
output_dir=$3
stim_dir=$4
task='epiproj'
block_time=10;
model_type='BLOCK'

module load afni
mkdir -p $output_dir
cd $output_dir

# do our regression 
3dDeconvolve                                                              \
    -force_TR 1.355                                                       \
    -input ${input_dir}/${file}_LR_surf_subcort_222_32k_fsLR_smooth1_scale.nii   \
    -polort A -float                                                   \
    -num_stimts 6                                                      \
    -stim_times 1 $stim_dir/${file}_futureself.1D "${model_type}(${block_time},1)"         \
    -stim_label 1 futureself                                                               \
    -stim_times 2 $stim_dir/${file}_presentself.1D "${model_type}(${block_time},1)"        \
    -stim_label 2 presentself                                              		            \
    -stim_times 3 $stim_dir/${file}_pastself.1D "${model_type}(${block_time},1)"            \
    -stim_label 3 pastself                                                            \
    -stim_times 4 $stim_dir/${file}_futurenonself.1D "${model_type}(${block_time},1)"        \
    -stim_label 4 futurenonself                                                                 \
    -stim_times 5 $stim_dir/${file}_presentnonself.1D "${model_type}(${block_time},1)"        \
    -stim_label 5 presentnonself                                                                 \
    -stim_times 6 $stim_dir/${file}_pastnonself.1D "${model_type}(${block_time},1)"        \
    -stim_label 6 pastnonself                                                              \
    -gltsym "SYM: pastself -presentself"                                                   \
    -glt_label 1 past_present_self                                                         \
    -gltsym "SYM: futureself -presentself"                                                 \
    -glt_label 2 future_present_self                                                       \
    -gltsym "SYM: futureself +pastself +presentself -futurenonself -pastnonself -presentnonself"  \
    -glt_label 3 self-nonself                                                            \
    -jobs 8                                                                                 \
    -fout -tout                                                \
    -x1D ${output_dir}/${file}_X.xmat.1D                    \
    -xjpeg ${output_dir}/${file}_X.jpg                         \
    -fitts ${output_dir}/${file}_fitts.nii                     \
    -errts ${output_dir}/${file}_errts.nii                  \
    -bucket ${output_dir}/${file}_stats.nii                 \
    -overwrite      \

df=$(3dAttribute BRICK_STATAUX ${output_dir}/${file}_stats.nii'[26]' | awk '{print $NF}')
echo $df

3dcalc \
    -a ${output_dir}/${file}_stats.nii \
    -expr "fitt_t2z (a,${df})" \
    -prefix ${output_dir}/${file}_zstats.nii \
    -float \
    -overwrite 