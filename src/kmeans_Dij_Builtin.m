function [centr, newVoxelIndices, newVoxelCoords, clusters, runtime] = kmeans_Dij_Builtin(coordinates, d_mat, k)
% KMEANS_DIJ Simple k-means function with random start - clusters
% on Dij of voxels
%
% Inputs:
%   > coordinates: All voxel coordinates
%   > d_mat_orig: Original Dij
%   > k: number of clusters
%
% Updated by Danielle Ripsman June 19, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

    %format long
    timerVal = tic;

    totalVoxels = size(d_mat,1);
    nBeamlets = size(d_mat,2);

    % Initialize our k-means with k random voxels
    newVoxelCoords = zeros(k,3);
    centr = zeros(k,nBeamlets);

    clusters = kmeans(d_mat,k);

    for kk = 1:k
        newVoxelCoords(kk,:) = mean(coordinates(clusters==kk,:),1);
        centr(kk,:) = mean(d_mat(clusters==kk,:),1);
    end
    
    % Remove all unassigned clustered
    unrealThings = isnan(sum(newVoxelCoords,2));
    newVoxelCoords(unrealThings,:) = [];
    centr(unrealThings,:) = [];
    
    actualClusters = size(centr,1);
    newVoxelIndices = 1:actualClusters;

    fprintf('There are actually k=%d clusters\n', actualClusters);

    runtime = toc(timerVal);
end
