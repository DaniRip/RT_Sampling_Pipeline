function plot_clusters_updated(coordinates, clusters, k)
    %%
        % Plots clusters of voxels. Expects an array of cluster assignments
        % A prettier plotter than the other
    %%
    colors = lines(k);

    for i = 1: k
        ind = clusters==i;
        voxels_x = coordinates(ind, 1);
        voxels_y = coordinates(ind, 2);
        voxels_z = coordinates(ind, 3);
        scatter3(voxels_x, voxels_y, voxels_z, 50, colors(i, :), 'filled');
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