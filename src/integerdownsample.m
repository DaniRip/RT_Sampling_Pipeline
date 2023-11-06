function [d_target, d_OAR, num_target_voxels, num_OAR_voxels, voxel_coord_target, voxel_coord_OAR, num_beamlets, target_dose, runtime] = integerdownsample(mat_file, sample_interval)
    timerVal = tic;
    input = load(mat_file);
    input = input.guiInput;
    num_beamlets = size(input.Dij,2);
    target_dose = input.targetDose;
    %beam_width = input.beamWidth;
    %d_beam_indicies = input.beamIndicies;
    voxel_coord_target = input.voxelIndices.target(1:sample_interval:end, :);
    voxel_coord_OAR = input.voxelIndices.OAR(1:sample_interval:end, :);
    
    d_target=input.Dij(input.structVoxels.target(1:sample_interval:end), :);
    d_OAR=input.Dij(input.structVoxels.OAR(1:sample_interval:end), :);
    
    num_target_voxels = size(d_target,1);
    num_OAR_voxels = size(d_OAR,1);
    runtime = toc(timerVal);
end
