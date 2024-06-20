function [centr, newVoxelIndices, newVoxelCoords, clusters, runtime] = kmeans_DijNeighbour(coordinates, d_mat, k)
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

    % Initialize our voxels randomly
    randPerm = randperm(totalVoxels);
    clusters = mod(randPerm, k) + 1;

    centr = zeros(k,nBeamlets);
    for kk = 1:k
        centr(kk,:) = mean(d_mat(clusters==kk,:),1);
    end

    % Going to assume things are still roughly cubic, so everything has
    % about 26 nearest neighbours. Except the boundaries, so there might be
    % a few extra checks here.

    NNmatrix = zeros(totalVoxels,26);
    for ii = 1:totalVoxels
        % The lowest value will always be 0, i.e., current voxel, so drop it
        [~,inds] = mink(sum((coordinates(ii,:) - coordinates).^2, 2),27);
        NNmatrix(ii,:) = inds(2:end)';
    end

    newVoxelCoords = zeros(k,3);

    stillChanging = true;
    iter = 0;
    
    while stillChanging 
        stillChanging = false;
        maxMinDist = 0;
        for ii = 1:totalVoxels
            relCentres = unique(clusters(NNmatrix(ii,:)))';
            % Switched to non-normed for speed, should be equivalent (can't sqrt 2 nums and have the bigger one be smaller!)
            nonNormedDist = sum((d_mat(ii,:) - centr(relCentres,:)).^2,2);%sqrt(sum((d_mat(ii,:) - centr).^2, 2));
            [curMin,kk] = min(nonNormedDist);
            if clusters(ii)~=relCentres(kk)
                if curMin>maxMinDist; maxMinDist = curMin; end
                clusters(ii)= relCentres(kk);
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
        if mod(iter, 10)==0
            fprintf('On iteration %d, current min %f\n', iter, maxMinDist);
        end
    end

    fprintf('Algorithm took %d iterations\n', iter);


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
