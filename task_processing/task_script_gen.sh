#!/bin/tcsh -xef

input_dir=$1
file=$2
output_dir=$3
stim_dir=$4
task=$5
block=$6
condition1=$7
condition2=$8

module load afni

mkdir -p $output_dir
cd $output_dir

3dDeconvolve                                                              \
    -force_TR 1.355                                                       \
    -input ${input_dir}/${file}_LR_surf_subcort_222_32k_fsLR_smooth1_scale.nii -polort A -float  \
    -num_stimts 2                                                      \
    -stim_times 1 ${stim_dir}/${file}_${condition1}.1D "BLOCK(${block},1)"        \
    -stim_label 1 ${condition1}                                                            \
    -stim_times 2 ${stim_dir}/${file}_${condition2}.1D "BLOCK(${block},1)"        \
    -stim_label 2 ${condition2}                                       \
    -gltsym "SYM: ${condition1} -${condition2}"                     \
    -glt_label 1 ${condition1}-${condition2}                          \
    -jobs 8                                                            \
    -fout -tout                                                        \
    -x1D ${output_dir}/${file}_X.xmat.1D                    \
    -xjpeg ${output_dir}/${file}_X.jpg                                     \
    -fitts ${output_dir}/${file}_fitts.nii                     \
    -errts ${output_dir}/${file}_errts.nii                  \
    -bucket ${output_dir}/${file}_stats.nii             \
    -overwrite                                       \

df=$(3dAttribute BRICK_STATAUX ${output_dir}/${file}_stats.nii'[8]' | awk '{print $NF}')
echo $df

3dcalc \
    -a ${output_dir}/${file}_stats.nii \
    -expr "fitt_t2z (a,${df})" \
    -prefix ${output_dir}/${file}_zstats.nii \
    -float \
    -overwrite

