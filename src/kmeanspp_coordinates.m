function [target_voxels_cluster, OAR_voxels_cluster, target_cluster_centr, OAR_cluster_centr, runtime] = kmeanspp_coordinates(mat_file, k)
    %%
       % Produces clusters of voxels separated using kmeans algorithm. 
       % Clusters are of type hashmap ["coordinates_centroid": [struct_voxels]
    %%

    format long
    timerVal = tic;
    [~, ~, num_target_voxels, num_OAR_voxels, voxel_coord_target, voxel_coord_OAR, ~, ~, ~] = integerdownsample(mat_file, 1);
    assert(k <= num_OAR_voxels && k <= num_target_voxels, 'Number of voxels need to be greater than the number of clusters')

    state = false;
    count = 0; % for target and OAR 

    % [{cluster_centroid: [list_of_voxels]}]
    target_voxels_cluster = containers.Map;
    OAR_voxels_cluster = containers.Map;
    % [{cluster: centroid}]
    target_cluster_centr = containers.Map;
    OAR_cluster_centr = containers.Map;

    parpool('local');
    pool_obj = gcp;

    while count < 2
        if state == false
            coordinates = voxel_coord_target;
        else
            coordinates = voxel_coord_OAR;
        end

        coordinates = transpose(coordinates);
        pool_constant = parallel.pool.Constant(coordinates); % broadcast variable for parfor loop

        % The k-means++ initialization (with coordinates)
        centr = coordinates(:,1+round(rand*(size(coordinates,2)-1)));
        cluster_assignment = ones(1,size(coordinates,2));
        for i = 2:k
            D = coordinates-centr(:,cluster_assignment);
            D = cumsum(sqrt(dot(D,D,1)));
            if D(end) == 0, centr(:,i:k) = coordinates(:,ones(1,k-i+1)); return; end
            centr(:,i) = coordinates(:,find(rand < D/D(end),1)); % rand = (0,1) is 
            % checking whether a randomly generated number between 0 and 1 is 
            % less than the corresponding probability for each data point.
            % Higher the probability, higher the chances that the rand
            % condition will hold true. A new cluster_assignment gets chosen
            [~,cluster_assignment] = max(bsxfun(@minus,2*(centr'*coordinates),dot(centr,centr,1).'));
        end
        
        cluster_prev = 0;
        iterations = 0;
        stop = false;
        while stop == false
            %% stopping criterion
            if cluster_prev == cluster_assignment
                stop = true;
                if state == false
                    state = true;
                end
            else
                %% Update cluster centers
                parfor i = 1:k
                    centr(:, i) = mean(pool_constant.Value(:, cluster_assignment == i), 2);
                end
                [~,cluster_assignment] = max(bsxfun(@minus,2*(centr'*coordinates),dot(centr,centr,1).'),[],1);
            end
            cluster_prev = cluster_assignment;
            iterations = iterations + 1;
        end
        %% Add voxels/centroids to each cluster_assignment (target and OAR)
        for i = 1:k
            if count == 0
                target_voxels_cluster(num2str(i)) = find(cluster_assignment == i);
                target_cluster_centr(num2str(i)) = centr(:, i);
            else
                OAR_voxels_cluster(num2str(i)) = find(cluster_assignment == i);
                OAR_cluster_centr(num2str(i)) = centr(:, i);
            end
        end
        count = count + 1;
    end
    runtime = toc(timerVal);
    delete(pool_obj);
end
