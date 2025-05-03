%% COLLECTMATRADDATAFORPLANNING_MULTI 
% is a script that takes output from the RT planning interface matRad 
% (https://github.com/e0404/matRad) and transforms it to be used for the 
% sampling interface (or other % applications!)
% 
% Use:
% Step1. run matRad in this instance of matlab and load a patient file 
%           a) This can be a DICOM .dcm from a local clinic or found online
%           b) Or one of the CORT /phantoms already processed for matRad
%           c) Or any other dataset processed into a matRad .mat file
% Step2. Run matRad's Dij calculator (i.e., Calc. Influence Mx) and note:
%           a) set the parameters how you want them, re:bixel width, # and
%           location of angles, etc.
%           b) IF you want open beams, matRad_generateStf approx line 154 
%           will need to be editted, as by default it fits a target 
%           projection.
%           c) You need to run the dose calculator on ALL the structures of
%           interest. Populate them into the objectives/constraints before
%           running or they might not show up when you need them.
% Step3. hardcode the extraction components below (names should match cst)
%           Variable: target_name, oar_names, outputFileName, targetDose
%           See matRad's prostate, head & neck, liver examples below
% Step4. >> collectmatRadDataForPlanning_multi
%
% The fields of the output guiInput structure are: 
%  > Dij: voxels x beamlets matrix
%  > structVoxels: substruct with voxel indices of each organ in Dij
%  > targets: target organ id(s) in structVoxels ordering
%  > OAR: sensitive organ id(s) in structVoxels ordering
%  > targetDose: Prescribed dose
%  > beamWidth: horizontal dimension of beamlet
%  > beamIndicies: [row, col, angle] indices for each voxel, sorted by row
%  > voxelIndices: I believe the [x,y,z] coord of each voxel
%  > voxelDim:[x,y,z] dimensions of the voxels (seems like we dont need
%  this anymore)
%
% This script was prepared by Stoyan Hristov (stoyan.hristov@uwaterloo.ca) 
% and updated by Danielle Ripsman (danielle.ripsman@rotman.utoronto.ca)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract structVoxels:

% Requirements: A loaded matRad instance, with a cst, dij, stf, pln struct
%               Access to the matRad codebase for functions, 
%               e.g. matRad_resizeCstToGrid
%               Note, the dij must be run with the desired organs, or they
%               may not exist!
% Uncomment the desired phantom to run, or try your own!

% LIVER Phantom 
target_name = 'PTV';
oar_names = {'Liver', 'SpinalCord'};
outputFileName = 'craftLIV_guiInput_multioar';
targetDose = 45.0000;

% PROSTATE Phantom
% target_name = 'PTV_68';
% oar_names = {'Rectum', 'Penile_bulb', 'Rt femoral head','Bladder','Lt femoral head'};
% outputFileName = 'craftPROS_guiInput_multioar';
% targetDose = 68.0000;

% HEAD AND NECK Phantom
% target_name = 'PTV70';
% oar_names = {'BRAIN_STEM', 'CEREBELLUM', 'CHIASMA', 'LARYNX', 'LENS_LT',...
%              'LENS_RT', 'LIPS', 'OPTIC_NRV_LT', 'OPTIC_NRV_RT',...
%              'PAROTID_LT', 'PAROTID_RT', 'SPINAL_CORD', 'TEMP_LOBE_LT',...
%              'TEMP_LOBE_RT', 'TM_JOINT_LT', 'TM_JOINT_RT'};  
% outputFileName = 'craftHN_guiInput_multioar';
% targetDose = 70.0000;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remap structures:

% First, map the cst to the space of the Dij 
vXgridOld = dij.ctGrid.x;  
vYgridOld = dij.ctGrid.y;  
vZgridOld = dij.ctGrid.z;  
vXgridNew = dij.doseGrid.x; 
vYgridNew = dij.doseGrid.y;  
vZgridNew = dij.doseGrid.z;  
cst_mapped = matRad_resizeCstToGrid(cst, vXgridOld, vYgridOld, vZgridOld, vXgridNew, vYgridNew, vZgridNew);

% Now, find the appropriate indices for each structure 
targetRow = find(strcmp(cst_mapped(:,2),target_name));
if isempty(targetRow)
    error('No row in cst_mapped has target_name in the second column.');
end 
targetData = cst_mapped{targetRow,4};
if ~iscell(targetData) || ~isvector(targetData{1}) || size(targetData{1},2) ~= 1
    error("Fourth column in the target row does not contain the expected x-by-1 double format")
end 
targetIndices = targetData{1};

% Initialize
oarIndices = [];
structVoxels = struct();
% Store absolute voxel indices
for i = 1:length(oar_names)
    oar_name = oar_names{i};
    oarRows = find(strcmp(cst_mapped(:,2), oar_name));
    if isempty(oarRows)
        warning('No rows in cst_mapped have "%s" in the second column.', oar_name);
        continue;
    end
    oarData = cst_mapped{oarRows,4};
    if ~iscell(oarData) || ~isvector(oarData{1}) || size(oarData{1},2) ~= 1
        warning('Unexpected format for voxel data for "%s". Skipping.', oar_name);
        continue;
    end
    fprintf("oar name: %s\n", oar_name);
    fprintf("oar rows: %d\n", oarRows);
    fprintf("oar data length: %d\n", numel(oarData{1}));
    theseIndices = oarData{1};
    oarIndices = [oarIndices; theseIndices];
    structVoxels.(matlab.lang.makeValidName(oar_name)) = theseIndices(:)';  % store absolute voxel indices
end

fprintf("oar indices length: %d\n", numel(oarIndices));

%%% Extract Dij matrix and slice for appropriate indices %%%
allIndices = unique([targetIndices; oarIndices],'stable');
oldDij = dij.physicalDose{1};

% Identify non-empty voxel indices
nonEmptyMask = any(oldDij(allIndices, :), 2);
filteredIndices = allIndices(nonEmptyMask);

for i = 1:length(oar_names)
    name = oar_names{i};
    orgVoxels = structVoxels.(matlab.lang.makeValidName(name));  % these are the voxels for that organ
    nMatched = sum(ismember(filteredIndices, orgVoxels));
    fprintf("%s: %d voxels matched\n", name, nMatched);
end

% Build a mapping from old index to new continuous index
newIndexMap = zeros(max(allIndices), 1);
newIndexMap(filteredIndices) = 1:length(filteredIndices);

% Map target indices
isTarget = ismember(allIndices, targetIndices);
validTargetIndices = unique(newIndexMap(allIndices(isTarget & nonEmptyMask)));
structVoxels.target = 1:numel(validTargetIndices);  % Sequentially enumerate target

% Map OARs
oarCombined = [];
for i = 1:length(oar_names)
    fName = matlab.lang.makeValidName(oar_names{i});
    if isfield(structVoxels, fName)
        originalVoxelIDs = structVoxels.(fName);  % absolute IDs
        keptMask = ismember(originalVoxelIDs, filteredIndices);
        newIDs = newIndexMap(originalVoxelIDs(keptMask));
        structVoxels.(fName) = newIDs(:)';
        oarCombined = [oarCombined, newIDs(:)']; % Choosing to keep duplicates here, because if they're in 2 organs, should be penalized twice!
    end
end
structVoxels.OAR = oarCombined(:)';

% Overwrite Dij matrix
newDij = oldDij(filteredIndices, :);

Dij = struct();

if isfield(dij, 'physicalDose') && iscell(dij.physicalDose)
    Dij = full(newDij);
else
    error('dij.physicalDose does not exist or is not in cell format')
end

%%% Extract targets and OAR %%%
%  > targets: target organ id(s) in structVoxels ordering
targets = 1;

%  > OAR: sensitive organ id(s) in structVoxels ordering
% will probably need to edit this to include multiple OARs
OAR = 2;

%%% Extract beamWidth %%%
%  > beamWidth: horizontal dimension of beamlet
if isfield(pln, 'propStf') && isfield(pln.propStf, 'bixelWidth')
    beamWidth = pln.propStf.bixelWidth / 10; % matRad is in mm, sampling app wants cm
else
    error('The field pln.propStf.bixelWidth does not exist.');
end

angles=[stf.totalNumOfBixels];
numAngles = numel(angles); % The number of angles corresponds to rows in stf

beamIndicies = zeros(sum(angles),3);

start = 1;

for a = 1:numAngles
    nRays = numel(stf(a).ray);
    beamCoords = zeros(nRays, 3);
    for i = 1:nRays
        beamCoords(i, :) = stf(a).ray(i).rayPos_bev;
    end
    % % Preallocate arrays for row and column indices
    % rowIndices = zeros(nRays, 1);
    % colIndices = zeros(nRays, 1);
    % 
    % % uniqueX = sort(unique(beamCoords(:, 1)), 'ascend'); % Columns: lowest to highest x
    % % uniqueY = sort(unique(beamCoords(:, 3)), 'descend'); % Rows: highest to lowest y
    % 
    % % Map each coordinate to its grid position
    % for i = 1:nRays
    %     rowIndices(i) = find(uniqueY == beamCoords(i, 3)); % Match y to row
    %     colIndices(i) = find(uniqueX == beamCoords(i, 1)); % Match x to column
    % end
    % This ordering might not carry properly, since the beams are different
    % sizes at different angles. Furthermore, if we're not resequencing the
    % dij columns by the beamlet locations, this will have to be
    % post-processed anyhow. Report accurate locations in cm instead!
    
    % Combine row and column indices
    beamIndicies(start:start+angles(a)-1,:,:) = [beamCoords(:, 3)/10, beamCoords(:, 1)/10, repmat(a,nRays,1)];%[rowIndices, colIndices,repmat(a,nRays,1)];
    start = start+angles(a);
end


%%% Generate voxelIndices %%%
originalTargetIndices = allIndices(isTarget & nonEmptyMask);  % This gives the actual 1D dose grid indices
targetTable = matRad_cubeIndex2worldCoords(originalTargetIndices, dij.doseGrid);

isOAR = ~isTarget;
originalOARIndices = allIndices(isOAR & nonEmptyMask);
oarTable = matRad_cubeIndex2worldCoords(originalOARIndices, dij.doseGrid);

% write to structure 
voxelIndices = struct();
voxelIndices.target = targetTable; 
voxelIndices.OAR = oarTable;

%%% Create guiInput structure and add substructures, save to .mat %%%
guiInput = struct();
guiInput.Dij = Dij;
guiInput.structVoxels = structVoxels;
guiInput.targets = targets;
guiInput.OAR = OAR;
guiInput.targetDose = targetDose;
guiInput.beamWidth = beamWidth;
guiInput.beamIndicies = beamIndicies;
guiInput.voxelIndices = voxelIndices;

%%% Save guiInput to .mat %%%
save(outputFileName, 'guiInput');
disp(['Data saved to ' outputFileName]);
