function [target_cluster, OAR_cluster, runtime] = surfaceSampler_Gatik_v3(mat_file)
    %%
       % Produces clusters of voxels with layers working inwards. Each
       % cluster(layer) consists of list of indices of struct_voxels.
       % Note: If you notice an isolated boundary point in an alphaShape, try
       % increasing the alpha value.
    %%
    timerVal = tic;    
    [~, ~, ~, ~, voxel_coord_target, voxel_coord_OAR, ~, ~, ~] = integerdownsample(mat_file, 1);

    state = false;
    count = 0; % for target and OAR

    while count < 2
        if state == false
            coordinates = voxel_coord_target;
        else
            coordinates = voxel_coord_OAR;
        end
        points = coordinates;
        layer = 1;
        cur_length = size(points, 1);
        prev_length = 0;
        while prev_length ~= cur_length
            prev_length = cur_length;
            x = points(:, 1);
            y = points(:, 2);
            z = points(:, 3);
            alpha_value = 2.0; % Pick max voxel dimension from (x, y, z)
            shp = alphaShape(x, y, z, alpha_value);
            boundary_facets = boundaryFacets(shp);
            boundary_points_indices = unique(boundary_facets(:));
            layer_voxels(layer).Vox = boundary_points_indices;
            points(boundary_points_indices, :) = [];
            layer = layer + 1;
            cur_length = size(points, 1);
            disp(size(points));
            if size(points, 1) == 48
                break;
            end
        end
    
        % % Plot the alpha shape and points on the surface
        % figure;
        % plot(shp);
        % hold on;
        % points_on_boundary = points(boundary_points_indices, :);
        % plot3(points_on_boundary(:, 1), points_on_boundary(:, 2), points_on_boundary(:, 3), 'r.', 'MarkerSize', 10);
        % zmin = -35;
        % zmax = -12;
        % zlim([zmin, zmax]);
        % title('Alpha Shape and Points on Surface');
        % xlabel('X-axis');
        % ylabel('Y-axis');
        % zlabel('Z-axis');
        % hold off;
    
        target_cluster = containers.Map;
        OAR_cluster = containers.Map;
        for layer = 1:length(layer_voxels)
            layer_str = num2str(layer);
            if count == 0
                target_cluster(layer_str) = layer_voxels(layer).Vox;
            else
                OAR_cluster(layer_str) = layer_voxels(layer).Vox;
            end
        end
        state = true;
        count = count + 1;
    end
    runtime = toc(timerVal);
end
    