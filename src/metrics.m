function [output] = metrics(mat_file)
    input = load(mat_file);
    input = input.guiInput;
    num_beamlets = size(input.Dij,2);
    target_dose = input.targetDose;
    beam_width = input.beamWidth;
    d_beam_indicies = input.beamIndicies;
    voxel_indicies = input.voxelIndicies;

    %average deviation from limit - how do I calculate this?

    %max overage - how do I calculate this?

    %# of constraints satisfied - handled in run_FMO.cpp

    %# of voxels underdosed - how do I calculate this?

end
