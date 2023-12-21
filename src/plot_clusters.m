function plot_clusters(target_or_OAR, cluster, mat_file)
    %%
        % Plots clusters of voxels. Expects a hashmap cluster with values
        % as list of struct_voxels
    %%
    [~, ~, ~, ~, voxel_coord_target, voxel_coord_OAR, ~, ~, ~] = integerdownsample(mat_file, 1);
    if target_or_OAR == "target"
        coordinates = voxel_coord_target;
    else
        coordinates = voxel_coord_OAR;
    end
    colors = lines(length(keys(cluster)));

    for k = 1:length(keys(cluster))
        key = keys(cluster);
        idx = cluster(key{k});

        % Voxel indices
        voxels_x = coordinates(idx, 1);
        voxels_y = coordinates(idx, 2);
        voxels_z = coordinates(idx, 3);

        scatter3(voxels_x, voxels_y, voxels_z, 50, colors(k, :), 'filled');
        hold on
    end

    hold off
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    % zmax = 0;
    % zlim([-50, zmax]);
    title('Plot of Clusters');
end