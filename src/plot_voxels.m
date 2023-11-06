function plot_voxels(target_or_OAR, mat_file)
    [~, ~, ~, ~, voxel_coord_target, voxel_coord_OAR, ~, ~, ~] = integerdownsample(mat_file, 1);
    if target_or_OAR == "target"
        coordinates = voxel_coord_target;
    else
        coordinates = voxel_coord_OAR;
    end
    X = coordinates(:, 1);
    Y = coordinates(:, 2);
    Z = coordinates(:, 3);

    plot3(X, Y, Z, 'o');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Plot of Voxels');
end