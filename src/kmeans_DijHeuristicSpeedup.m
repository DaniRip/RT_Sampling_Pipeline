function [centr, newVoxelIndices, newVoxelCoords, clusters, runtime] = kmeans_DijHeuristicSpeedup(coordinates, d_mat, k)
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

    disp('Starting k-mD Heuristic')

    %format long
    timerVal = tic;

    totalVoxels = size(d_mat,1);
    nBeamlets = size(d_mat,2);

    % Initialize our k-means with k random voxels
    randomClusters = randperm(totalVoxels, k);
    centr = d_mat(randomClusters,:);

    clusters = zeros(totalVoxels,1);
    clusters(randomClusters) = 1:k; % assign the assigned clusters

    fullInds = 1:k;

    newVoxelCoords = zeros(k,3);

    stillChanging = true;
    iter = 0;
    
    while stillChanging && iter < 75
        stillChanging = false;
        maxMinDist = 0;
        rowSumsCentr = sum(centr,2);% Trying something out as a speed boost
        rowSumsDmat = sum(d_mat,2);
        for ii = 1:totalVoxels
            [~,inds] = mink(abs(rowSumsDmat(ii) - rowSumsCentr),20); % 20 Smallest deviations in totality
            % Switched to non-normed for speed, should be equivalent (can't sqrt 2 nums and have the bigger one be smaller!)
            nonNormedDist = sum((d_mat(ii,:) - centr(inds,:)).^2,2);%sqrt(sum((d_mat(ii,:) - centr).^2, 2));
            [curMin,kk] = min(nonNormedDist);
            if clusters(ii)~=inds(kk)
                if curMin>maxMinDist; maxMinDist = curMin; end
                clusters(ii)= inds(kk);
                stillChanging=true;
            end
        end
        if~stillChanging % We can already end, nothing changed
            break;
        end
        for kk = 1:k
            centr(kk,:) = mean(d_mat(clusters==kk,:),1);
        end
        iter = iter+1;
        if mod(iter, 5)==0
            fprintf('On iteration %d, current min %f\n', iter, maxMinDist);
        end
    end

    fprintf('Algorithm took %d iterations\n', iter);
    if (stillChanging); disp('Stopped due to timeout'); end


    for kk = 1:k
        newVoxelCoords(kk,:) = mean(coordinates(clusters==kk,:),1);
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
