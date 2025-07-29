function generate_cifti_rotations(input_cifti, num_iterations, output_dir, num_pools)

% This function (adapted from EMG 2020 PNAS) takes in a CIFTI as
% input, generates num_iterations rotations of vertices in the file, and 
% saves out all rotated maps into a singular output CIFTI.
%
%   input_cifti = '/path/to/cifti.dtseries.nii' 
%       **this should be a single-map CIFTI, otherwise script will use map #1
%   num_iterations = number of times to perform random rotation of input
%   output_dir = folder to output the rotated maps (e.g., 'Rotations1000.dtseries.nii')
%   num_pools (OPTIONAL) = default=6; number of parallel processors to use
%
%   OUTPUTS are:
%   1) one .dtseries.nii file in output_dir, with num_iterations maps
%   2) a .mat with (num_vertices x num_iterations) matrix with vertex
%   rotation indices
%
%   NOTE - medial wall vertices are typically excluded in our CIFTIs.
%   Rotated maps will reflect this as 0s in certain locations.
%

% set up variables
%set parallel pools to 6 if not specified
if nargin==3
    num_pools=6;
end

%read input to be rotated
if ischar(input_cifti)
    cifti_template = ft_read_cifti_mod(input_cifti);
    [~, fname, ~] = fileparts(input_cifti);
else
    error('Input must be path to CIFTI filename')
end

% format coordinates to rotate
sphere = gifti('/projects/b1081/Scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.coord.gii');
sphereLcoords = sphere.vertices;
sphere = gifti('/projects/b1081/Scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.coord.gii');
sphereRcoords = sphere.vertices;
allspherecoords = [sphereLcoords ; sphereRcoords];

ncortverts = nnz(cifti_template.brainstructure(cifti_template.brainstructure>0)<3);
nsurfverts = nnz(cifti_template.brainstructure<3);
nsurfvertsL = nsurfverts ./ 2;
nsurfvertsR = nsurfverts ./ 2;

ciftidata_orig = cifti_template.data(1:ncortverts, 1);
cifti_template.data((ncortverts+1):end,:) = [];
cifti_template.pos(cifti_template.brainstructure>2,:) = [];
cifti_template.brainstructure(cifti_template.brainstructure>2) = [];
cifti_template.brainstructurelabel(3:end) = [];

origspherecoords = allspherecoords(cifti_template.brainstructure>0 & cifti_template.brainstructure<3,:);
origspherehems = cifti_template.brainstructure(cifti_template.brainstructure>0 & cifti_template.brainstructure<3);
brainstructure_surfacespaceL = cifti_template.brainstructure(1:nsurfvertsL);
brainstructure_surfacespaceR = cifti_template.brainstructure((nsurfvertsL+1):(nsurfvertsL + nsurfvertsR));

% start parallel workers
if isempty(gcp('nocreate'))
    pool = parpool(num_pools);
end

% generate rotated maps
remapped = zeros(size(cifti_template.data,1), num_iterations);
rot_maps = zeros(size(cifti_template.data,1), num_iterations);

parfor iter = 1:num_iterations
    
    disp(iter)
    test = zeros(length(brainstructure_surfacespaceL),3);
    
    ok = 0;
    while ok==0
        
        xrot = rand * 2*pi;
        yrot = rand * 2*pi;
        zrot = rand * 2*pi;
        
        %i think this is implementing a minimum rotation amount?
        if (min([xrot-0, 2*pi - xrot]) + min([yrot-0, 2*pi - yrot]) + min([zrot-0, 2*pi - zrot])) > (pi/8)
            ok = 1;
        end
    end
    
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    xrotcoords = rotmat_x * origspherecoords';
    xyrotcoords = rotmat_y * xrotcoords;
    xyzrotcoords = rotmat_z * xyrotcoords;
    
    yrot = -yrot;
    zrot = -zrot;
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    xrotcoords = rotmat_x * origspherecoords';
    xyrotcoords = rotmat_y * xrotcoords;
    xyzrotcoords_mirror = rotmat_z * xyrotcoords;
    
    this_remapped = zeros(size(cifti_template.data,1),1);
    
    for ind = 1:ncortverts
        
        if origspherehems(ind)==1
            test(:,1) = xyzrotcoords(1,ind); test(:,2) = xyzrotcoords(2,ind); test(:,3) = xyzrotcoords(3,ind);
            [~, rot_ind] = min(sum(abs(sphereLcoords - test),2));
            if brainstructure_surfacespaceL(rot_ind) > 0
                rot_cifti_ind = rot_ind - nnz(brainstructure_surfacespaceL(1:rot_ind)<0);
                this_remapped(ind) = rot_cifti_ind;
            end
        else
            test(:,1) = xyzrotcoords_mirror(1,ind); test(:,2) = xyzrotcoords_mirror(2,ind); test(:,3) = xyzrotcoords_mirror(3,ind);
            [~, rot_ind] = min(sum(abs(sphereRcoords - test),2));
            if brainstructure_surfacespaceR(rot_ind) > 0
                rot_cifti_ind = rot_ind - nnz(brainstructure_surfacespaceR(1:rot_ind)<0) + nnz(brainstructure_surfacespaceL>0);
                this_remapped(ind) = rot_cifti_ind;
            end
        end
    end

    remapped(:,iter) = this_remapped;
    non_medialwall_inds = this_remapped ~= 0; %should still be fine if there are 0s in original datafile
    this_rot_map = zeros(size(cifti_template.data,1),1);
    this_rot_map(non_medialwall_inds) = ciftidata_orig(this_remapped(non_medialwall_inds));
    rot_maps(:,iter) = this_rot_map;
end

out_filename = strtok(fname, '.');
out = cifti_template;
% out.data = remapped;
% out.dimord = 'pos_time';
% ft_write_cifti_mod([output_dir '/Rotations_' num2str(num_iterations) 'iter_indices.dtseries.nii'], out);
save([output_dir '/' out_filename '_ROTATIONS_' num2str(num_iterations) 'perms_indices.mat'], 'remapped')
out.data = rot_maps;
ft_write_cifti_mod([output_dir '/' out_filename '_ROTATIONS_' num2str(num_iterations) 'perms.dtseries.nii'], out);

poolobj = gcp('nocreate');
delete(poolobj)

end

