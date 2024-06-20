function [d_mat, final_voxel_coords, clusters, runtime] = fountainclustering(d_mat, k, type)
% FOUNTAINCLUSTERING Cluster by dominant beamlet, then find relevant
% beamlets from that cluster.
%
% Inputs:
%   > d_mat: Any structure's Dij matrix
%   > n: Number of samples desired
%   > type: either - max or med. 
%       - max = take the max values from each cluster (up to n voxels)
%       - med = take the median voxels from each cluster (up to n voxels)
%
% Updated by Danielle Ripsman June 17, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

    timerVal = tic;

    if ~exist('k','var') || isempty(k); k = size(d_mat,2); end
    if ~exist('type','var') || isempty(type); type = 'max'; end

    final_voxel_coords = zeros(k,1);

    % Clustering the voxels based on their maximum beamlets:
    [~,clusters] = max(d_mat,[],2);
    
    % Add voxels to each cluster_assignment (target and OAR)
    unique_clusters = unique(clusters); % set of all unique beamlet indices with highest doses
    totalClusters = numel(unique_clusters);

    % histc is outdated, so we'll do this more confusing version:
    counts = histcounts(clusters, [unique_clusters; max(clusters)+1]');
    %histc(clusters, unique_clusters);

    % Calculate the number of items to pull from each cluster - try to
    % assign things systematically, only account for items that there are
    % enough of in the cluster, as there isn't enough, remove from pool
    includeCounts = zeros(totalClusters,1);
    idx = 1:totalClusters;
    iter = 1;
    totalCounted = 0;
    while totalCounted < k
        includeCounts(idx(1:min(k-totalCounted,length(idx))))=iter;
        totalCounted = totalCounted + min(k-totalCounted,length(idx));
        idx(counts(idx) == iter)=[];
        iter=iter+1;
    end

    % Find and sort the voxels we are keeping
    voxInd = 1;
    for j = 1:totalClusters
        currCluster = find(clusters==unique_clusters(j));
        [~, vSeq] = sort(d_mat(currCluster,unique_clusters(j)),'descend');
        if strcmp(type, 'max')
            final_voxel_coords(voxInd:voxInd+includeCounts(j)-1)=currCluster(vSeq(1:includeCounts(j)));
        elseif strcmp(type, 'med')
            % Take the midpoint, then subtract # of include values from it.
            % Both floored down so we bias backwards once, forwards once.
            medPoint = max(1,floor(numel(currCluster)/2)-floor(includeCounts(j)/2));
            final_voxel_coords(voxInd:voxInd+includeCounts(j)-1)=currCluster(vSeq(medPoint:medPoint+includeCounts(j)-1));
        end
        voxInd=voxInd+includeCounts(j);
    end

    final_voxel_coords=sort(final_voxel_coords);
    d_mat = d_mat(final_voxel_coords,:);

    runtime = toc(timerVal);
end
