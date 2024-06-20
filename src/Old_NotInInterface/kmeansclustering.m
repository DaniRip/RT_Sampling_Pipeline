function [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, runtime] = kmeansclustering(mat_file, k, sample_interval)
    %kMeans Clusters data points into k clusters.
    %   Input args: k: number of clusters; 
    %   Output args: cluster: 1-by-n array with values of 0,...,k-1
    %   representing in which cluster the corresponding point lies in
    %   centr: m-by-k matrix of the m-dimensional centroids of the k clusters

    timerVal = tic;
    [d_target, d_OAR, ~, ~, num_beamlets, target_dose, ~] = integerdownsample(mat_file, sample_interval);
    state = false;
    count = 0;
    while count < 2
        if state==false
            P = d_target;
        else
            P = d_OAR;
        end
        
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