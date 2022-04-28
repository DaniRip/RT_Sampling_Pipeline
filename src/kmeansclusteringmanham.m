function [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, runtime] = kmeansclusteringmanham(mat_file, k, sample_interval)
    %kMeans Clusters data points into k clusters.
    %   Input args: k: number of clusters;
    %   Output args: cluster: 1-by-n array with values of 0,...,k-1
    %   representing in which cluster the corresponding point lies in
    %   centr: m-by-k matrix of the m-dimensional centroids of the k clusters

    timerVal = tic;
    [d_target, d_OAR, ~, ~, num_beamlets, target_dose, ~] = integerdownsample(mat_file, sample_interval);
    state = false;
    count = 0;
    
    % these currently aren't downsampled!
    target_indicies = input.voxelIndicies.target;
    OAR_indicies = input.voxelIndicies.OAR;
    indicies = [target_indicies;OAR_indicies];
    
    while count < 2
        if state==false
            P = d_target;
        else
            P = d_OAR;
        end

        P = transpose(P);
        numP = size(P,2); %voxels
        dimP = size(P,1); %beamlets

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
            % for each voxel, find if its closer to a neghbour's cluster
            for idxP = 1:numP
                
                dist = zeros(1,2);
                
                self_cluster = cluster(indxP);
                self_dist = norm(P(:,idxP)-centr(:,self_cluster));
                % change this to radius of nearby indicies
                neighbours = [idxP-1 idxP+1];
                for count = 1:length(neighbours)
                    dist(count) = norm(P(:,idxP)-centr(:,cluster(neighbours(count))));
                    if dist(count) < self_dist
                        self_cluster = cluster(neighbours(count));
                    end
                end
                cluster(idxP) = self_cluster;
            end
            
            % calculate centroids of each cluster
            centr = zeros(dimP,k);
            for idxC = 1:k
                centr(:,idxC) = mean(P(:,cluster==idxC),2);
            end

            if clusterPrev==cluster
                stop = true;
                cluster = transpose(cluster);
                centr = transpose(centr);
                if state==false
                    d_target = centr;
                    state = true;
                else
                    d_OAR = centr;
                end
            end
            clusterPrev = cluster;

            iterations = iterations + 1;
        end
        count = count + 1;
        %fprintf('Used %d iterations of changing centroids.\n',iterations);
    end
    num_target_voxels = size(d_target,1);
    num_OAR_voxels = size(d_OAR,1);
    runtime = toc(timerVal);
end
