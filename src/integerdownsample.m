function [d_mat, final_voxel_indices, runtime] = integerdownsample(d_mat, n)
% INTEGERDOWNSAMPLE A quick function to return every n voxels
%
% Updated by Danielle Ripsman June 16, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

    timerVal = tic;

    voxel_indices = 1:(size(d_mat,1));
    final_voxel_indices = voxel_indices(1:n:end);
    d_mat=d_mat(final_voxel_coords, :);

    runtime = toc(timerVal);

end
