function plot_voxels(vox_coordinates)
% INTEGERDOWNSAMPLE A quick function to return every n voxels
%
% Written by Gatik Gola Nov. 2023
% Updated by Dani June 16, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

X = vox_coordinates(:, 1);
Y = vox_coordinates(:, 2);
Z = vox_coordinates(:, 3);

plot3(X, Y, Z, '.');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Plot of Voxels');

end