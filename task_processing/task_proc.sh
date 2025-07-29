
input_dir=$1
input_file=$2
output_dir=$3

module load afni

mkdir -p $output_dir
cd $output_dir

3dTstat -prefix ${input_file}_mean.nii ${input_dir}/${input_file}.nii
                                           
3dcalc                                                                 \
    -a ${input_dir}/${input_file}.nii                                              \
    -b ${input_file}_mean.nii                                                \
    -expr 'min(200, a/b*100)*step(a)*step(b)'                      \
    -prefix ${input_file}_scale.nii                                 \
    -overwrite

rm ${output_dir}/${input_file}_mean.nii 
