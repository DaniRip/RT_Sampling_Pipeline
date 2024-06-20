function [d_mat, newVoxelIndices, centr, clusters, runtime] = kmeans_coordinates(coordinates, d_mat_orig, k)
% KMEANS_COORDINATES Simple k-means function with random start - clusters
% on coordinates of voxels
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

    totalVoxels = size(coordinates,1);
    nBeamlets = size(d_mat_orig,2);

    % Initialize our k-means with k random voxels
    randomClusters = randperm(totalVoxels, k);
    centr = coordinates(randomClusters,:);

    clusters = zeros(totalVoxels,1);
    clusters(randomClusters) = 1:k; % assign the assigned clusters

    d_mat = zeros(k,nBeamlets);

    stillChanging = true;
    iter = 0;
    
    while stillChanging 
        stillChanging = false;
        maxMinDist = 0;
        for ii = 1:totalVoxels
            % Switched to non-normed for speed, should be equivalent (can't sqrt 2 nums and have the bigger one be smaller!)
            nonNormedDist = sum((coordinates(ii,:) - centr).^2, 2);%normDist = sqrt(sum((coordinates(ii,:) - centr).^2, 2));
            [curMin,kk] = min(nonNormedDist);
            if clusters(ii)~=kk
                if curMin>maxMinDist; maxMinDist = curMin; end
                clusters(ii)= kk;
                stillChanging=true;
            end
        end
        if~stillChanging % We can already end, nothing changed
            break;
        end
        for kk = 1:k
            centr(kk,:) = mean(coordinates(clusters==kk,:),1);
        end
        iter = iter+1;
        if mod(iter, 10)==0
            fprintf('On iteration %d, current min %f\n', iter, maxMinDist);
        end
    end

    fprintf('Algorithm took %d iterations\n', iter);

    for kk = 1:k
        d_mat(kk,:) = mean(d_mat_orig(clusters==kk,:),1);
    end

    % Remove all unassigned clustered
    unrealThings = isnan(sum(d_mat,2));
    d_mat(unrealThings,:) = [];
    centr(unrealThings,:) = [];
    
    actualClusters = size(d_mat,1);
    newVoxelIndices = 1:actualClusters;

    fprintf('There are actually k=%d clusters\n', actualClusters);

    runtime = toc(timerVal);
end
