function visualizeStructures(structs)
% VISUALIZESTRUCTURES plots the structures of interest on a grid
%
% Input: 
%   > structs: A structure containing the structures with a list of voxel
%   indices
%
% Written by: Danielle Ripsman daripsman@uwaterloo.ca
% Last edited: April 7, 2022
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

fn = fieldnames(structs);
numFields = length(fn);

cols = [1,0,0;0,0,1;0,1,0];

figure
for ii=1:numFields
    scatter3(structs.(fn{ii})(:,1), structs.(fn{ii})(:,2),structs.(fn{ii})(:,3), 'CData', cols(mod(ii,3)+1,:), 'Marker', '.')
    hold on
end
