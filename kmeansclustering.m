function [ cluster, centr ] = kmeansclustering(mat_file, k )

%kMeans Clusters data points into k clusters.
%   Input args: k: number of clusters; 
%   Output args: cluster: 1-by-n array with values of 0,...,k-1
%   representing in which cluster the corresponding point lies in
%   centr: m-by-k matrix of the m-dimensional centroids of the k clusters

input = load(mat_file);
P = input.guiInput;
num_beamlets = size(input.Dij,2);
target_dose = input.targetDose;
beam_width = input.beamWidth;
d_beam_indicies = input.beamIndicies;

P = transpose(P);
numP = size(P,2); % number of points
dimP = size(P,1); % dimension of points

% choose k unique random indices between 1 and size(P,2) (number of points)
randIdx = randperm(numP,k);
% initial centroids
centr = P(:,randIdx);

% init cluster array
cluster = zeros(1,numP);

% init previous cluster array clusterPrev (for stopping criterion)
clusterPrev = cluster;

% for reference: count the iterations
iterations = 0;

% init stopping criterion
stop = false; % if stopping criterion met, it changes to true

while stop == false
    
    % for each data point 
    for idxP = 1:numP
        % init distance array dist
        dist = zeros(1,k);
        % compute distance to each centroid
        for idxC=1:k
            dist(idxC) = norm(P(:,idxP)-centr(:,idxC));
        end
        % find index of closest centroid (= find the cluster)
        [~, clusterP] = min(dist);
        cluster(idxP) = clusterP;
    end
    
    % Recompute centroids using current cluster memberships:
        
    % init centroid array centr
    centr = zeros(dimP,k);
    % for every cluster compute new centroid
    for idxC = 1:k
        % find the points in cluster number idxC and compute row-wise mean
        centr(:,idxC) = mean(P(:,cluster==idxC),2);
    end
    
    % Checking for stopping criterion: Clusters do not chnage anymore
    if clusterPrev==cluster
        stop = true;
    end
    % update previous cluster clusterPrev
    clusterPrev = cluster;
    
    iterations = iterations + 1;
    
end

fprintf('kMeans.m used %d iterations of changing centroids.\n',iterations);
end