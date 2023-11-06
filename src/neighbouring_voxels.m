function [voxels_idx] = neighbouring_voxels(queryIdx, coordinates, kdtree)
    queryVoxel = coordinates(queryIdx, :);
    [voxels_idx, ~] = knnsearch(kdtree, queryVoxel, 'K', 21);
end


