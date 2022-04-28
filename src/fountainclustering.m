function [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, runtime] = fountainclustering(mat_file, sample_interval)
    timerVal = tic;
    [d_target, d_OAR, ~, ~, num_beamlets, target_dose, ~] = integerdownsample(mat_file, sample_interval);
    state = false;
    count = 0;
    
    % go through each row, cluster each voxel based off its highest 
    while count < 2
        % we do this for both matrices
        if state==false
            P = d_target;
        else
            P = d_OAR;
        end
        
        numP = size(P,2); %beamlets
        dimP = size(P,1); %voxels

        % init cluster array
        cluster = zeros(1,numP);
        
        for idxP = 1:dimP
            [~,idx] = max(P(idxP,:));
            cluster(idxP) = idx;
        end
        
        % end up with a some#x1520 array where each column has the indices
        % of the voxels clustered to that specific beamlet
        % take the first voxel for each cluster, make a new (some#less than
        % 1520_x1520 matrix to pass to FMO
        
        if state==false
            d_target = cluster;
        else
            d_OAR = cluster;
        end
        
        state = true;
        count = count + 1;
    end
    num_target_voxels = size(d_target,1);
    num_OAR_voxels = size(d_OAR,1);
    runtime = toc(timerVal);
end
