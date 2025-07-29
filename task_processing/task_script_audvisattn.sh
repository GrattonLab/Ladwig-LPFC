#!/bin/tcsh -xef

input_dir=$1
file=$2
output_dir=$3
stim_dir=$4
cue_time=2.6
block_time=26

module load afni
mkdir -p $output_dir
cd $output_dir

# do our regression 
3dDeconvolve                                                              \
    -force_TR 1.355                                                       \
    -input ${input_dir}/${file}_LR_surf_subcort_222_32k_fsLR_smooth1_scale.nii   \
    -polort A -float                                                   \
    -num_stimts 6                                                      \
    -stim_times 1 $stim_dir/${file}_watchright.1D "BLOCK(${block_time},1)"        \
    -stim_label 1 watchright                                               	   \
    -stim_times 2 $stim_dir/${file}_watchleft.1D "BLOCK(${block_time},1)"        \
    -stim_label 2 watchleft                                                           \
    -stim_times 3 $stim_dir/${file}_listenright.1D "BLOCK(${block_time},1)"        \
    -stim_label 3 listenright                                               		   \
    -stim_times 4 $stim_dir/${file}_listenleft.1D "BLOCK(${block_time},1)"        \
    -stim_label 4 listenleft                                                            \
    -stim_times 5 $stim_dir/${file}_passive.1D "BLOCK(${block_time},1)"        \
    -stim_label 5 passive                                                           \
    -stim_times 6 $stim_dir/${file}_cue.1D "BLOCK(${cue_time},1)"        \
    -stim_label 6 cue                                                           \
    -gltsym "SYM: watchright +watchleft -listenright -listenleft"                     \
    -glt_label 1 watch-listen                                               \
    -gltsym "SYM: watchright +watchleft -passive"                        \
    -glt_label 2 watch-passive                   \
    -gltsym "SYM: listenright +listenleft -passive"                        \
    -glt_label 3 listen-passive                   \
    -jobs 8                                                            \
    -fout -tout                                                     \
    -x1D ${output_dir}/${file}_X.xmat.1D                    \
    -xjpeg ${output_dir}/${file}_X.jpg                                     \
    -fitts ${output_dir}/${file}_fitts.nii                     \
    -errts ${output_dir}/${file}_errts.nii                  \
    -bucket ${output_dir}/${file}_stats.nii                 \
    -overwrite \

df=$(3dAttribute BRICK_STATAUX ${output_dir}/${file}_stats.nii'[23]' | awk '{print $NF}')
echo $df

3dcalc \
    -a ${output_dir}/${file}_stats.nii \
    -expr "fitt_t2z (a,${df})" \
    -prefix ${output_dir}/${file}_zstats.nii \
    -float \
    -overwrite
