function [d_mat, final_voxel_coords, runtime] = surfaceSampler_Gatik_v3(d_mat, coordinates, nLayer, k)
%%
   % Produces clusters of voxels with layers working inwards. Each
   % cluster(layer) consists of list of indices of struct_voxels.
   % Note: If you notice an isolated boundary point in an alphaShape, try
   % increasing the alpha value.
%%
% Inputs:
%   > d_mat: Any structure's Dij matrix
%   > coordinates: The x-y-z coordinates over every voxel in d_mat
%   > nLayer: The sampling rate of an every n-layer sampling
%   > k: Total number of desired voxels for a reverse sampling approach
%   Note: either nLayer or k may be empty, but not both!
%
% Written by Gatik Gola Nov. 2023
% Updated by Dani June 18, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


    timerVal = tic;    

    % Build out the layers:
    points = coordinates;
    layer = 0;
    indices = 1:size(points, 1);
    indices=indices';
    alpha_value = 3.0;%1.5;%2.0; % Pick max voxel dimension from (x, y, z)
    while ~isempty(points)
        layer = layer + 1;
        x = points(:, 1);
        y = points(:, 2);
        z = points(:, 3);
        shp = alphaShape(x, y, z, alpha_value);
        boundary_facets = boundaryFacets(shp);
        boundary_points_indices = unique(boundary_facets(:));
        if isempty(boundary_points_indices)
            % Put the remaining points in the last layer
            points = [];
            indices = [];
        else 
            layer_voxels(layer).Vox = indices(boundary_points_indices);
            points(boundary_points_indices, :) = [];
            indices(boundary_points_indices)=[];
        end
    end

    nTotal=0;
    final_voxel_coords=[];
    % Sample every n layers:
    if ~isempty(nLayer)
        for ii = 1:nLayer:layer
            final_voxel_coords =[final_voxel_coords; layer_voxels(ii).Vox];
            nTotal = nTotal+numel(layer_voxels(ii).Vox);
        end
    else
        % If we're targetting a specific number of layers, 
        isFull = false;
        % Calculate the size of each vector using arrayfun
        layerSizes = arrayfun(@(x) numel(x.Vox), layer_voxels);
        usedIndices = zeros(1,numel(layerSizes));
        while ~isFull
            currentFit = layerSizes<k-nTotal;
            currentFit(usedIndices==1)=0;
            if sum(currentFit)==0
                isFull=true;
            else 
                [~,currInd] = find(currentFit);
                final_voxel_coords =[final_voxel_coords; layer_voxels(currInd(1)).Vox];
                usedIndices(currInd(1))=1;
                nTotal = nTotal+layerSizes(currInd(1));
            end
        end
    end

    final_voxel_coords=sort(final_voxel_coords);
    d_mat = d_mat(final_voxel_coords,:);

    runtime = toc(timerVal);
end
    