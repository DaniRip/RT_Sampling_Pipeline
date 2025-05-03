The fields of the guiInput structure 
  > Dij: voxels x beamlets matrix
  > structVoxels: voxel indices of each organ in Dij
  > targets: target organ id(s) in structVoxels ordering
  > OAR: sensitive organ id(s) in structVoxels ordering
  > targetDose: Prescribed dose
  > beamWidth: horizontal (and vertical) dimension of beamlet
  > beamIndicies: [row, col, angle] rowwise and columnwise location of beamlet (BEV) centre (cm), angle #
  > voxelIndices: I believe the [x,y,z] coord of each voxel
