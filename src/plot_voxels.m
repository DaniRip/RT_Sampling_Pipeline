function plot_voxels(vox_coordinates, type, col)
% INTEGERDOWNSAMPLE A quick function to return every n voxels
%
% Written by Gatik Gola Nov. 2023
% Updated by Dani June 16, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


if (~exist('col','var')|| isempty(col)); col = 'b'; end
if (~exist('type','var')|| isempty(type)); col = ''; end

figure

X = vox_coordinates(:, 1);
Y = vox_coordinates(:, 2);
Z = vox_coordinates(:, 3);

plot3(X, Y, Z, '.')%, 'Color', col);
xlabel('X');
ylabel('Y');
zlabel('Z');
title(strcat(['Plot of Sampled Voxels',type]));

end