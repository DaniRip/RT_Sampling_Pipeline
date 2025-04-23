function [d_mat, sampled_voxels_indices, runtime] = probingDoseOptimization(d_mat, target_dose, max_iter, k)
% Probing-based Dose Optimization with Subset Sampling
% Uses projected gradient descent to find the most important voxels
% 
% Inputs:
%   d_mat           - dose-influence matrix
%   target_dose     - target dose vector (n x 1)
%   max_iter        - maximum number of iterations
%   coordinates     - voxel coordinates matrix (n x 3) 
%   k               - number of sampled voxels required
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-% 
    timerVal = tic;
    % Dimensions
    [n, nb] = size(d_mat);

    % Step 1: Probing (Projected Gradient Descent)
    x = zeros(nb, 1);
    [~,D] = eig(2*(d_mat')*d_mat);
    [d,~] = sort(diag(D));
    min_eig = d(1);
    max_eig = d(end);
    eta = 2/(min_eig + max_eig);
    tol = 1e-4;
    grad_f = 2 * d_mat' * (d_mat * x - target_dose);
    iter = 0;
    for i = 1:max_iter
        % Compute gradient
        grad_f = 2 * d_mat' * (d_mat * x - target_dose);
        % Gradient update with projection onto non-negative orthant
        x = x - eta * grad_f;
        x = max(x, 0);  % Element-wise projection
        iter = iter + 1;
    end
    x_probe = x;

    % Step 2: Compute Importance Scores
    r = d_mat * x_probe - target_dose;  % Residual
    importance_scores = zeros(n, 1);
    
    for v = 1:n
        % Compute gradient norm for each voxel
        row_v = d_mat(v, :)';
        grad_v = 2 * row_v * r(v);
        importance_scores(v) = norm(grad_v);
    end

    % Step 3: Sampling
    % Normalize importance scores
    q = importance_scores / sum(importance_scores);
    
    % Custom implementation of weighted sampling without replacement
    sampled_voxels_indices = sample_weighted_scores(n, k, q);
    % sampled_voxels_coord = coordinates(sampled_voxels_indices, :);

    d_mat = d_mat(sampled_voxels_indices,:);

    runtime = toc(timerVal);
    % Note: Visualize the sampled voxels using plot_voxels function:
    % plot_voxels(sampled_voxels_coord)
end

function sampled_indices = sample_weighted_scores(n, k, weights)
    % Sample the top k weight scores
    % Inputs:
    %   n       - total number of items
    %   m       - number of items to sample
    %   weights - probability weights for sampling
    
    % Validate inputs
    assert(length(weights) == n, 'Weights must match number of items');
    assert(k <= n, 'Cannot sample more items than available');
   
    [~, original_indices] = sort(weights, 'descend');
    
    sampled_indices = zeros(k, 1);
    
    % Select top m items based on weights
    sampled_indices = original_indices(1:k);
end