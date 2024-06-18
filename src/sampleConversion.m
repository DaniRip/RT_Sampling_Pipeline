function [n,k,p] = sampleConversion(nToAll, kToAll, pToAll, value, totalVoxels)
% SAMPLECONVERSION A function to convert values from one sampling rate to
% all others, except the layered approach, since this involves
% post-processing.
%
% Inputs:
%   > _ToAll: Binary value to indicate which metric is being input (only 1)
%   > value: The value of that k/n/p input 
%   > totalVoxels: Total number of voxels under sampling consideration
%
% Updated by Danielle Ripsman June 16, 2024
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% Input check:
if sum([nToAll, kToAll, pToAll]) > 1 || sum([kToAll, nToAll, pToAll]) < 1
    error('invalid input');
end

if nToAll
    n = value;
    k = floor(totalVoxels/n);
    p = round(1/n*100,2);
end

if pToAll
    p = value;
    n = floor(100/value);
    k = floor(totalVoxels/value/100);
end

if pToAll
    k = value;
    n = floor(totalVoxels/k);
    p = round(k/totalVoxels*100);
end