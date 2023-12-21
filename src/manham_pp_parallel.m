function [target_voxels_cluster, OAR_voxels_cluster, target_cluster_centr, OAR_cluster_centr, runtime] = manham_pp_parallel(mat_file, k)
    %%
       % Produces clusters of voxels separated using kmeans mahnam algorithm. 
       % Clusters are of type hashmap ["centroid": [struct_voxels]
    %%

    timerVal = tic;    
    [d_target, d_OAR, num_target_voxels, num_OAR_voxels, voxel_coord_target, voxel_coord_OAR, ~, ~, ~] = integerdownsample(mat_file, 1);
    assert(k <= num_OAR_voxels && k <= num_target_voxels, 'Number of voxels need to be greater than the number of clusters')
    
    state = false;
    count = 0; % for target and OAR 

    % [{cluster: [list_of_voxels]}]
    target_voxels_cluster = containers.Map;
    OAR_voxels_cluster = containers.Map;
    % [{cluster: centroid}]
    target_cluster_centr = containers.Map;
    OAR_cluster_centr = containers.Map;
    
    parpool('local');
    pool_obj = gcp;
    
    while count < 2
        if state == false
            P = d_target;
            coordinates = voxel_coord_target;
        else
            P = d_OAR;
            coordinates = voxel_coord_OAR;
        end

        P = transpose(P);
        % broadcast variable for parfor loop so it doesn't get reinitialized 
        % for every parallel process
        pool_constant = parallel.pool.Constant(P);
        [~, numP] = size(P); % [beamlets, voxels]

        % The k-means++ initialization (with P)
        centr = P(:,1+round(rand*(numP-1)));
        cluster = ones(1,numP);
        for i = 2:k
            D = P-centr(:,cluster);
            D = cumsum(sqrt(dot(D,D,1)));
            if D(end) == 0, centr(:,i:k) = P(:,ones(1,k-i+1)); return; end
            centr(:,i) = P(:,find(rand < D/D(end),1)); % rand = (0,1) is 
            % checking whether a randomly generated number between 0 and 1 is 
            % less than the corresponding probability for each data point.
            % Higher the ratio, higher the chances that the rand
            % condition for that index will hold true.
            [~,cluster] = max(bsxfun(@minus,2*(centr'*P),dot(centr,centr,1).'));
        end

        cluster_prev = 0;

        kdtree = KDTreeSearcher(coordinates); % for nearest neighbour search
        iterations = 0;
        stop = false;

        %% Precalculate neighbouring voxels
        nei_voxels_idx = cell(1, numP);
        parfor idxP = 1:numP
            nei_voxels_idx{idxP} = neighbouring_voxels(idxP, coordinates, kdtree);
        end
        
        while stop == false
            for idxP = 1:numP
                cluster_idx = cluster(idxP); % current cluster
                %% Check with neighbouring cluster                
                for nei_voxel_idx = nei_voxels_idx{idxP}
                    cluster_centr = centr(:, cluster_idx); % current cluster centroid
                    nei_cluster_idx = cluster(nei_voxel_idx);
                    nei_cluster_centr = centr(:, nei_cluster_idx);
                    d_nei_diff = P(:, idxP) - nei_cluster_centr;
                    d_nei = dot(d_nei_diff, d_nei_diff); % dot-product is equivalent to squared 2-norm in R^n
                    d_centr_diff = P(:, idxP) - cluster_centr;
                    d_centr = dot(d_centr_diff, d_centr_diff);
                    if d_nei < d_centr
                        cluster_idx = nei_cluster_idx;
                    end
                end
                cluster(idxP) = cluster_idx;  % assign voxel to the updated cluster
            end
            %% stopping criterion
            if cluster_prev == cluster
                stop = true;
                if state == false
                    state = true;
                end
            else
                %% update cluster centres
                parfor i = 1:k
                    centr(:, i) = mean(pool_constant.Value(:, cluster == i), 2);
                end
            end
            cluster_prev = cluster;
            iterations = iterations + 1;
        end
        %% Add voxels/centroid to each cluster (target and OAR)
        for i = 1:k
            if count == 0
                target_voxels_cluster(num2str(i)) = find(cluster == i);
                target_cluster_centr(num2str(i)) = centr(:, i);
            else
                OAR_voxels_cluster(num2str(i)) = find(cluster == i);
                OAR_cluster_centr(num2str(i)) = centr(:, i);
            end
        end
        count = count + 1;
        %fprintf('Used %d iterations of changing centroids.\n',iterations);
    end
    runtime = toc(timerVal);
    delete(pool_obj);
end
