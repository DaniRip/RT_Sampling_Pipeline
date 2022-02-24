function [ cluster, centr ] = kmeansclustering(mat_file, k )

%kMeans Clusters data points into k clusters.
%   Input args: k: number of clusters; 
%   Output args: cluster: 1-by-n array with values of 0,...,k-1
%   representing in which cluster the corresponding point lies in
%   centr: m-by-k matrix of the m-dimensional centroids of the k clusters

%input = load(mat_file);
P = mat_file;
%num_beamlets = size(input.Dij,2);
%target_dose = input.targetDose;
%beam_width = input.beamWidth;
%d_beam_indicies = input.beamIndicies;

P = transpose(P);
numP = size(P,2);
dimP = size(P,1);

% choose k unique random indices between 1 and size(P,2) (number of points)
randIdx = randperm(numP,k);
centr = P(:,randIdx);

% init cluster array
cluster = zeros(1,numP);

% init previous cluster array clusterPrev (for stopping criterion)
clusterPrev = cluster;

iterations = 0;
stop = false;

while stop == false
    
    
    for idxP = 1:numP
        
        dist = zeros(1,k);
        
        for idxC=1:k
            dist(idxC) = norm(P(:,idxP)-centr(:,idxC));
        end
        
        [~, clusterP] = min(dist);
        cluster(idxP) = clusterP;
    end
    
    centr = zeros(dimP,k);
    for idxC = 1:k
        centr(:,idxC) = mean(P(:,cluster==idxC),2);
    end
    
    if clusterPrev==cluster
        stop = true;
    end
    clusterPrev = cluster;
    
    iterations = iterations + 1;
    
end

fprintf('kMeans.m used %d iterations of changing centroids.\n',iterations);
end