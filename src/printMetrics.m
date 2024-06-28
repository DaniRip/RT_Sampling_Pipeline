function printMetrics(sol, d_target, d_OAR, voxelDim)
% METRICS extracts dose and time metrics from output plans
% 
% Input:
%   > d_target: The target dose-influence matrix
%   > d_OAR: The OAR dose-influence matrix
%   > sol: solution or set of solutions of the form sol.samplingtype.outputs
%   > voxelDim: [x,y,z] in mm if exists for the DXcc calculations
%
% Output:
%   > Two formulated-for-LaTeX sets of outputs
%   Table 1: Run times (sampling and FMO model)
%   Table 2: Dose info
%       - Avg_H: Average OAR dose
%       - Max_H: Max OAR dose
%       - Avg_T: Average tumor dose
%       - Min_T: Minimum tumor dose
%       - Max_T: Maximum tumor dose
%       - D0.5%: 99.5 percentile - dose delivered to the hottest .5% of the target
%       - D95%: 5th percentile - dose delivered to 95% of the target volume
%       - D99%: 1st percentile - dose delivered to 99% of the target volume
%
% Edited By: Danielle Ripsman    
% Last Edited: June 18, 2024
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-%
  
% Get all the field names (output names)
outputNames = fieldnames(sol);

fprintf('The output times for LaTeX are:\n')
fprintf('Name & Sample Time (sec) & Model Time (sec) & Actual Num Vox \\\\ \\hline\n');
% Loop through all the output names
for i = 1:length(outputNames)   
    % Extract the numeric value from the string
    FMOruntime = sscanf(sol.(outputNames{i}).FMOruntime, 'FMO Runtime = %f sec');
    sampRuntime = sol.(outputNames{i}).sampRuntime;
    if isfield(sol, 'target_voxels')
        numClusters = sol.target_voxels;
    else
        numClusters = size(d_target,1);
    end
    fprintf('%s & %.3f & %.3f & %d \\\\ \n', outputNames{i}, round(sampRuntime, 4), round(FMOruntime, 4), numClusters);
end

if ~isempty(voxelDim)
    rowNames = {'& Avg$_\cH$', '& D1cc', '& D10cc', '& Avg$_\cT$', ...
            '& Min$_\cT$', '& D$0.5\%$', '& D$95\%$', '& D$99\%$'};
else
    rowNames = {'& Avg$_\cH$', '& Max$_\cH$', '& Avg$_\cT$', ...
            '& Min$_\cT$', '& Max$_\cT$', '& D$0.5\%$', '& D$95\%$', '& D$99\%$'};
end

title = 'Metric';
matSummary = zeros(length(rowNames),length(outputNames));

fprintf('\n\nThe output metrics for LaTeX are:\n')
% For each type of method
for i = 1:length(outputNames)   
    % Table header:    
    title = strcat([title, ' & ' outputNames{i}]);

    count = 1;

    % Assign row-values
    v = d_OAR*sol.(outputNames{i}).w';
    matSummary(count,i) = mean(v); count = count+1;
    if isempty(voxelDim)
        matSummary(count,i) = max(v); count = count+1;
    else
        ccLoc = ceil(1/(voxelDim(1)*voxelDim(2)*voxelDim(3)));
        tenccLoc = ceil(10/(voxelDim(1)*voxelDim(2)*voxelDim(3)));
        v = sort(v,'descend');
        matSummary(count,i) = v(ccLoc); count = count+1;
        matSummary(count,i) = v(tenccLoc); count = count+1;
    end
        
    v = d_target*sol.(outputNames{i}).w';
    matSummary(count,i) = mean(v); count = count+1;
    matSummary(count,i) = min(v); count = count+1;

    if isempty(voxelDim)
        matSummary(count,i) = max(v); count = count+1;
    end

    numTvoxels = size(d_target,1);
    v = sort(v,'descend');
    %D0.5% 
    matSummary(count,i) = v(ceil(0.005*numTvoxels)-1); count = count+1;
    %D95% 
    matSummary(count,i) = v(ceil(0.95*numTvoxels)-1); count = count+1;
    %D99% 
    matSummary(count,i) = v(ceil(0.99*numTvoxels)-1);
end

fprintf(strcat(title,'\\\\\n'));

for i = 1:length(rowNames)
    % Start each row with the column name
    fprintf('%s & ', rowNames{i});
    
    % Add the numerical values, rounded to 4 decimal places
    for j = 1:length(outputNames)
        fprintf('%.3f', matSummary(i, j));
        if j < length(outputNames)
            fprintf(' & ');
        end
    end
    
    % End the row
    fprintf(' \\\\ \n');
end









