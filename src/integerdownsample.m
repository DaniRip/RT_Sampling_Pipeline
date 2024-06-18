function [d_mat, final_num_voxels, final_voxel_coords, runtime] = integerdownsample(d_mat, voxel_indices, n)
% INTEGERDOWNSAMPLE A quick function to return every n voxels
%
% Updated by Danielle Ripsman June 16, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

    timerVal = tic;

    final_voxel_coords = voxel_indices(1:n:end);
    d_mat=d_mat(final_voxel_coords, :);
    final_num_voxels = numel(final_voxel_coords);

    runtime = toc(timerVal);

end
