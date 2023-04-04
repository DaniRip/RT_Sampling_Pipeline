function [avg_deviation, max_overage, voxels_underdosed] = metrics(mat_file, fmo_output)
    input = load(mat_file);
    input = input.guiInput;
    %num_beamlets = size(input.Dij,2);
    %target_dose = input.targetDose;
    %beam_width = input.beamWidth;
    %d_beam_indicies = input.beamIndicies;
    %voxel_indicies = input.voxelIndicies;

    %average deviation from limit
    fmo_output = transpose(fmo_output);
    D_target=input.Dij(input.structVoxels.target,:);
    D_OAR=input.Dij(input.structVoxels.OAR,:);
    dv_target = D_target*fmo_output;
    dv_OAR = D_OAR*fmo_output;

    % Average deviation from tumor dose:
    avg_deviation = mean(abs(dv_target-input.targetDose));

    % Maximum deviation from tumor dose:
    max_overage = max(abs(dv_target-input.targetDose));

    % Number of tumor voxels underdosed:
    voxels_underdosed = sum(dv_target<input.targetDose+.00001);
