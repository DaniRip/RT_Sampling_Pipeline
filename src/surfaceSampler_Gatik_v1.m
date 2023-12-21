function [cluster] = surfaceSampler_Gatik(mat_file)
    timerVal = tic;    
    [d_target, d_OAR, ~, ~, voxel_coord_target, voxel_coord_OAR, ~, ~, ~] = integerdownsample(mat_file, 1);

    state = false;
    count = 0; % for target and OAR 

    if state == false
        P = d_target;
        coordinates = voxel_coord_target;
    else
        P = d_OAR;
        coordinates = voxel_coord_OAR;
    end
    P = transpose(P);

    %% Assigning voxels to grid cells
    X_coord = coordinates(:, 1);
    Y_coord = coordinates(:, 2);
    Z_coord = coordinates(:, 3);
    % TO DO: Find the best nbins for efficiency and accuracy of algorithm
    % in general cases
    nbins_x = 57;
    nbins_y = 175;
    nbins_z = 70;
    % bin allocation of data points for each axis
    [~, ~, bin_x] = histcounts(X_coord, nbins_x);
    [~, ~, bin_y] = histcounts(Y_coord, nbins_y);
    [~, ~, bin_z] = histcounts(Z_coord, nbins_z);

    mask = zeros(nbins_x, nbins_y, nbins_z); % discretized grid
    lin_ind = sub2ind([nbins_x nbins_y nbins_z], bin_x, bin_y, bin_z); % Convert bin coordinates to linear indices
    %% Maintain hashmap of struct voxels idx in each target grid cell
    % hashmap: linear_idx : [struct_voxel_idx]
    unique_lin_inds = unique(lin_ind);
    grid_cell_voxels = struct('vox', cell(1, max(unique_lin_inds)));
    for i = 1:numel(unique_lin_inds)
        unique_lin_ind = unique_lin_inds(i);
        [x, y, z] = ind2sub([nbins_x nbins_y nbins_z], unique_lin_ind); % (bin_x, bin_y, bin_z)
        mask(x, y, z) = 1;
        idx = find(bin_x == x & bin_y == y & bin_z == z);
        grid_cell_voxels(unique_lin_ind).vox = idx;
    end
    
    %% Fill the empty spots in mask array with artificial binary ones
    for i = 1:nbins_z
        count = 0;
        while count < 2
            row_matrix = mask(:, :, i);
            if count == 1 row_matrix = transpose(row_matrix); end
            [row_ind, col_ind] = find(row_matrix == 1);
            unique_row_ind = unique(row_ind);
            first_occurrence = arrayfun(@(r) find(row_ind == r, 1, 'first'), unique_row_ind);
            last_occurrence = arrayfun(@(r) find(row_ind == r, 1, 'last'), unique_row_ind);
            col_lower_bound = col_ind(first_occurrence);
            col_upper_bound = col_ind(last_occurrence);
            for j = 1:length(col_lower_bound)
                row_matrix(unique_row_ind(j), col_lower_bound(j):col_upper_bound(j)) = 1;
            end
            if count == 0
                mask(:, :, i) = row_matrix;
            else
                mask(:, :, i) = transpose(row_matrix);
            end
            count = count+1;
        end
    end
    for i = 1:nbins_x
        count = 0;
        while count < 2
            row_matrix = squeeze(mask(i, :, :));
            if count == 1 row_matrix = transpose(row_matrix); end
            [row_ind, col_ind] = find(row_matrix == 1);
            unique_row_ind = unique(row_ind);
            first_occurrence = arrayfun(@(r) find(row_ind == r, 1, 'first'), unique_row_ind);
            last_occurrence = arrayfun(@(r) find(row_ind == r, 1, 'last'), unique_row_ind);
            col_lower_bound = col_ind(first_occurrence);
            col_upper_bound = col_ind(last_occurrence);
            for j = 1:length(col_lower_bound)
                row_matrix(unique_row_ind(j), col_lower_bound(j):col_upper_bound(j)) = 1;
            end
            if count == 0
                mask(i, :, :) = squeeze(row_matrix);
            else
                mask(i, :, :) = transpose(squeeze(row_matrix));
            end
            count = count+1;
        end
    end
    for i = 1:nbins_y
        count = 0;
        while count < 2
            row_matrix = squeeze(mask(:, i, :));
            if count == 1 row_matrix = transpose(row_matrix); end
            [row_ind, col_ind] = find(row_matrix == 1);
            unique_row_ind = unique(row_ind);
            first_occurrence = arrayfun(@(r) find(row_ind == r, 1, 'first'), unique_row_ind);
            last_occurrence = arrayfun(@(r) find(row_ind == r, 1, 'last'), unique_row_ind);
            col_lower_bound = col_ind(first_occurrence);
            col_upper_bound = col_ind(last_occurrence);
            for j = 1:length(col_lower_bound)
                row_matrix(unique_row_ind(j), col_lower_bound(j):col_upper_bound(j)) = 1;
            end
            if count == 0
                mask(:, i, :) = squeeze(row_matrix);
            else
                mask(:, i, :) = transpose(squeeze(row_matrix));
            end
            count = count+1;
        end
    end
    %%
    cur_layer = 0;
    mapLength = 15;
    cmap = colormap(jet(4*mapLength));
    cmap = cmap(1:4:4*mapLength,:);
    cmap = flip(cmap,1);

    surf_points = getSurfacePoints(mask);
    surf_x = surf_points(:, 1);
    surf_y = surf_points(:, 2);
    surf_z = surf_points(:, 3);

    figure;
    p = patch(isosurface(mask)); % The second argument is the isovalue (0.5 for binary masks)
    p.FaceColor = 'blue'; % Set face color
    p.EdgeColor = 'none'; % Set edge color
    daspect([1,1,1]); % Adjust aspect ratio
    view(3); % Set the view to 3D
    axis tight; % Tighten the axis limits
    grid on; % Display grid

    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');

    title('3D Binary Mask Visualization with Isosurface');

    log = ismember([surf_x surf_y surf_z], [bin_x bin_y bin_z], 'rows');
    cluster = containers.Map;
    indices = find(log == 1);
    disp(size(indices));
    disp(indices(1: 100));
    cluster("1") = indices;
    

    % %% Assigning voxels to grid cells
    % X_coord = coordinates(:, 1);
    % Y_coord = coordinates(:, 2);
    % Z_coord = coordinates(:, 3);
    % % TO DO: Use binary search technique to find nbins for general tumour cases
    % % The idea is to find the maximum pixelation such that there are no
    % % empty bins (bins with no voxels) because that will affect the logic
    % % of surface point detection
    % nbins_x = 175;
    % nbins_y = 57;
    % nbins_z = 70;
    % % bin allocation of data points for each axis
    % [N_x, ~, bin_x] = histcounts(X_coord, nbins_x);
    % [N_y, ~, bin_y] = histcounts(Y_coord, nbins_y);
    % [N_z, ~, bin_z] = histcounts(Z_coord, nbins_z);
    % disp(N_x);
    % disp(N_y);
    % disp(N_z);
    % 
    % mask = zeros(nbins_z, nbins_x, nbins_y); % discretized grid
    % lin_ind = sub2ind([nbins_z nbins_x nbins_y], bin_z, bin_x, bin_y); % Convert bin coordinates to linear indices
    % %% Maintain hashmap of struct voxels idx in each target grid cell
    % unique_lin_inds = unique(lin_ind);
    % disp('unique_lin_inds')
    % disp(size(unique_lin_inds));
    % grid_cell_voxels = struct('vox', cell(1, max(unique_lin_inds)));
    % for i = 1:numel(unique_lin_inds)
    %     unique_lin_ind = unique_lin_inds(i);
    %     [z, x, y] = ind2sub([nbins_z nbins_x nbins_y], unique_lin_ind); % (bin_z, bin_x, bin_y)
    %     mask(z, x, y) = 1;
    %     idx = find(bin_z == z & bin_x == x & bin_y == y);
    %     grid_cell_voxels(unique_lin_ind).vox = idx;
    % end
    % 
    % %%
    % cur_layer = 0;
    % mapLength = 15;
    % cmap = colormap(jet(4*mapLength));
    % cmap = cmap(1:4:4*mapLength,:);
    % cmap = flip(cmap,1);
    % 
    % surf_points = getSurfacePoints(mask);
    % z_val = surf_points(:, 1);
    % x_val = surf_points(:, 2);
    % y_val = surf_points(:, 3);
    % 
    % idx = find(bin_z' == surf_points(:,1) & bin_x' == surf_points(:,2) & bin_y' == surf_points(:,3));
    % disp(size(idx));
    % disp(idx(1: 100));
    % figure;
    % p = patch(isosurface(mask)); % The second argument is the isovalue (0.5 for binary masks)
    % p.FaceColor = 'blue'; % Set face color
    % p.EdgeColor = 'none'; % Set edge color
    % daspect([1,1,1]); % Adjust aspect ratio
    % view(3); % Set the view to 3D
    % axis tight; % Tighten the axis limits
    % grid on; % Display grid
    % 
    % xlabel('X-axis');
    % ylabel('Y-axis');
    % zlabel('Z-axis');
    % 
    % title('3D Binary Mask Visualization with Isosurface');
    % %% Peel back layer by layer until there is nothing left
    % while sum(sum(sum(mask))) > 0
    %     cur_layer = cur_layer + 1;
    %     surf_points = getSurfacePoints(mask);
    %     disp(size(surf_points(:, 1)));
    %     disp(size(surf_points(:, 2)));
    %     disp(size(surf_points(:, 3)));
    %     idx = find(bin_z' == surf_points(:,1) & bin_x' == surf_points(:,2) & bin_y' == surf_points(:,3));
    %     outer_layers(cur_layer).vox = idx;
    %     % Plot every layer
    % %     hold on; p = patch(isosurface(mask));
    % %     p.EdgeColor ='none'; p.FaceAlpha = .2; p.FaceColor =cmap(mod(curLayer-1,mapLength)+1,:);
    %     mask(surf_points(:,1), surf_points(:,2), surf_points(:,3)) = 0;
    % end
    % %%
    % sample_voxels = [];
    % 
    % if cur_layer == 0
    %     error('no structure input');
    % elseif cur_layer == 1
    %     sample_voxels = outer_layers(cur_layer).vox;
    % else
    %     chosen_layers = [1,cur_layer];
    %     if cur_layer > 2 && num_layers > 2
    %         num_layers = min(cur_layer,num_layers); % in case less layers than requested exist
    %         all_layers = 2:cur_layer-1;
    %         index = 0:(cur_layer-2)/(num_layers-1):cur_layer-2;
    %         index = round(index(2:end));
    %         index = index(1:num_layers-2);
    %         chosen_layers = [chosen_layers,all_layers(index)];        
    %     end
    % 
    %     chosen_layers = sort(chosen_layers);
    % 
    %     % Go through layers backwards, since most inner was last
    %     axis;
    %     for ii = 1:length(chosen_layers)
    %        cur_ind = chosen_layers(end-ii+1); 
    %        if ii > 1 % Always starts at 1 for now
    %            for jj = cur_ind+1:chosen_layers(end-ii+2)-1
    %                mask(outer_layers(jj).vox)=1;
    %            end
    %        end
    %        mask(outer_layers(cur_ind).vox)=1;
    %        % Will plot every layer but the center if its too small but nice!
    %        hold on; p = patch(isosurface(mask));
    %        p.EdgeColor ='none'; p.FaceAlpha = .2; p.FaceColor =cmap(mod(ii-1,mapLength)+1,:);
    % %        blankMask = zeros(size(mask));
    % %        blankMask(OuterLayers(curInd).Vox) = 1;
    % %        [i,j,k] = ind2sub(size(mask),find(blankMask));
    % %        hold on; s = scatter3(i,j,k);
    % %        s.MarkerFaceColor = 'none'; s.MarkerEdgeColor = cmap(mod(curLayer-1,mapLength)+1,:); s.MarkerEdgeAlpha = .2;
    %        sample_voxels = [sample_voxels; outer_layers(cur_ind).vox];
    %     end
    % 
    % end
    % 
    % % Will print every layer, but not as nicely
    % mask = repmat(logical(0), scanSize);
    % mask(sample_voxels) = 1;
    % figure; hold on; p = patch(isosurface(mask));
    % p.EdgeColor ='none'; p.FaceAlpha = .2; p.FaceColor =cmap(mod(ii-1,mapLength)+1,:);
    % 
    % display(strcat([num2str(length(sample_voxels)), ' out of ', num2str(totalVoxels), ' selected (', num2str(totalVoxels/length(sample_voxels)), ' sampling)']));
    % sample_voxels = unique(sample_voxels);
end