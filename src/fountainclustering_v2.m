function [target_voxels_cluster, OAR_voxels_cluster, target_cluster_centr, OAR_cluster_centr, runtime] = fountainclustering_v2(mat_file)
    %%
       % Produces clusters of voxels grouped wrt dominant beamlets, bemalet
       % with the maximum intensity.
    %%
    timerVal = tic;
    [d_target, d_OAR, ~, ~, ~, ~, ~, ~, ~] = integerdownsample(mat_file, 1);
    state = false;
    count = 0;
    
    % [{cluster_beamlet: [list_of_voxels]}]
    target_voxels_cluster = containers.Map;
    OAR_voxels_cluster = containers.Map;
    % [{cluster_beamlet: centroid}]
    target_cluster_centr = containers.Map;
    OAR_cluster_centr = containers.Map;

    % go through each row, cluster each voxel based off its highest 
    while count < 2
        if state == false
            P = d_target;
        else
            P = d_OAR;
        end

        numP = size(P,2); % beamlets
        dimP = size(P,1); % voxels
        
        % init cluster array
        cluster = zeros(1,dimP); % each idx (voxel) in cluster belongs to cluster(idx)
        centr = cell(1, numP); % each value corresponds to {highest dose, voxel index}
        for i = 1:numel(centr)
            centr{i} = {0, 0};
        end
        
        for idxP = 1:dimP
            [~,idx] = max(P(idxP,:)); % find dominant beamlet
            cluster(idxP) = idx;
            x = isequal(centr{idx}, {0, 0});
            if x
                centr{idx} = {P(idxP, idx), idxP};
            else
                existing_val = centr{idx};
                bool = P(idxP, idx) > existing_val{1};
                centr{idx} = {bool*P(idxP, idx) + (1-bool)*existing_val{1},...
                    bool*idx + (1-bool)*existing_val{2}};
            end
        end
        %% Add voxels to each cluster_assignment (target and OAR)
        unique_clusters = unique(cluster); % set of all unique beamlet indices with highest doses
        for unique_cluster = unique_clusters
            if count == 0
                target_voxels_cluster(num2str(unique_cluster)) = find(cluster == unique_cluster);
                val = centr{unique_cluster};
                target_cluster_centr(num2str(unique_cluster)) = val{2};
            else
                OAR_voxels_cluster(num2str(unique_cluster)) = find(cluster == unique_cluster);
                val = centr{unique_cluster};
                OAR_cluster_centr(num2str(unique_cluster)) = val{2};
            end
        end
        state = true;
        count = count + 1;
    end
    runtime = toc(timerVal);
end
