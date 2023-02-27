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
    D_vb_w_b = input.Dij*fmo_output;
    deviation = 0;
    for i = 1:numel(D_vb_w_b)
        deviation = deviation + abs((input.targetDose - D_vb_w_b(i)));
    end
    avg_deviation = deviation/numel(D_vb_w_b);

    %max overage (max dose of tumor and max dose of OAR) - D_vb*w_b 
    %gives total dose given to each voxel, then substract the prescribed (
    %the right hand side of all constraints which is the input of the FMO model (usually called targetDose and value of 42.4)) 
    %dosage for each voxel, then take the maximum of that vector
    D_vb_w_b = input.Dij*fmo_output;
    max_overage = 0;
    for i = 1:numel(D_vb_w_b)
        if (D_vb_w_b(i)-input.targetDose) > max_overage
            max_overage = D_vb_w_b(i)-input.targetDose; 
        end
    end

    %# of constraints satisfied - leave it out for now
    %all the ones that don't satisfy the ___

    %# of voxels underdosed - how do I calculate this?
    %Multiply the original matrix with the FMO output, then step through
    %that result to find any voxels that are below the target dosage
    %CHANGE THIS TO ONLY COUNT TARGET VOXELS
    D_vb_w_b = input.Dij*fmo_output;
    voxels_underdosed = 0;
    for i = 1:numel(D_vb_w_b)
        if D_vb_w_b(i) < input.targetDose
            voxels_underdosed = voxels_underdosed + 1;
        end
    end
end
