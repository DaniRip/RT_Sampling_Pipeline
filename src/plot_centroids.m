function plot_centroids(target_or_OAR, cluster, centr, mat_file)
    %%
        % Plots centroids of clusters. Expects a hashmap cluster with values
        % as list of struct_voxels, centr
    %%
    [~, ~, ~, ~, voxel_coord_target, voxel_coord_OAR, ~, ~, ~] = integerdownsample(mat_file, 1);
    if target_or_OAR == "target"
        coordinates = voxel_coord_target;
    else
        coordinates = voxel_coord_OAR;
    end
    
    % Prepare for plotting the centroids
    centroids_x_beam1 = []; % beamlet idx <= 760
    centroids_x_beam2 = []; % beamlet idx > 760
    centroids_y_beam1 = [];
    centroids_y_beam2 = [];
    centroids_z_beam1 = [];
    centroids_z_beam2 = [];
    
    % fountain clustering
    if ~isempty(centr)
        for key = keys(centr)
            idx = centr(key{1});
            % centroid indices
            centr_x = coordinates(idx, 1);
            centr_y = coordinates(idx, 2);
            centr_z = coordinates(idx, 3);
            % Append centroid coordinates
            if str2double(key) <= 760
                centroids_x_beam1 = [centroids_x_beam1; centr_x];
                centroids_y_beam1 = [centroids_y_beam1; centr_y];
                centroids_z_beam1 = [centroids_z_beam1; centr_z];
            else
                centroids_x_beam2 = [centroids_x_beam2; centr_x];
                centroids_y_beam2 = [centroids_y_beam2; centr_y];
                centroids_z_beam2 = [centroids_z_beam2; centr_z];
            end
        end
    % any k-means clustering
    else
        for key = keys(cluster)
            idx = cluster(key{1});
            % voxels indices
            voxels_x = coordinates(idx, 1);
            voxels_y = coordinates(idx, 2);
            voxels_z = coordinates(idx, 3);
            % centroid indices
            centr_x = mean(voxels_x);
            centr_y = mean(voxels_y);
            centr_z = mean(voxels_z);
            
            % Append centroid coordinates
            centroids_x_beam1 = [centroids_x_beam1; centr_x];
            centroids_y_beam1 = [centroids_y_beam1; centr_y];
            centroids_z_beam1 = [centroids_z_beam1; centr_z];
        end
    end
        
    % Plot only the centroids
    scatter3(centroids_x_beam1, centroids_y_beam1, centroids_z_beam1, 100, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black');
    hold on
    scatter3(centroids_x_beam2, centroids_y_beam2, centroids_z_beam2, 100, 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
    hold off
    % zmax = 0;
    % zlim([-50, zmax]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Plot of Centroids');
end
