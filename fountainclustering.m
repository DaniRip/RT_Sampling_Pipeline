function [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, runtime] = fountainclustering(mat_file, sample_interval)
    timerVal = tic;
    [d_target, d_OAR, ~, ~, num_beamlets, target_dose] = integerdownsample(mat_file, sample_interval);
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
        
        P = transpose(P);
        numP = size(P,2); %voxels
        dimP = size(P,1); %beamlets

        % init cluster array
        cluster = zeros(1,numP);
        [~, idx] = max(P, [], 2);
        for idxP = 1:numP
            cluster(idxP) = idx(idxP);
        end
    end
    num_target_voxels = size(d_target,1);
    num_OAR_voxels = size(d_OAR,1);
    runtime = toc(timerVal);
end
