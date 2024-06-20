function plot_clusters(coordinates)
    %%
        % Plots clusters of voxels. Expects a hashmap cluster with values
        % as list of struct_voxels
        % A prettier plotter than the other
    %%

    numClusters = size(coordinates,1);
    colors = lines(numClusters);

    for k = 1:numClusters
        % Voxel indices
        voxels_x = coordinates(k, 1);
        voxels_y = coordinates(k, 2);
        voxels_z = coordinates(k, 3);

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