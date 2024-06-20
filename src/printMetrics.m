function printMetrics(sol, d_target, d_OAR)
% METRICS extracts dose and time metrics from output plans
% 
% Input:
%   > d_target: The target dose-influence matrix
%   > d_OAR: The OAR dose-influence matrix
%   > sol: solution or set of solutions of the form sol.samplingtype.outputs
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
fprintf('Name & Sample Time (sec) & Model Time (sec) \\\\ \\hline\n');
% Loop through all the output names
for i = 1:length(outputNames)   
    % Extract the numeric value from the string
    FMOruntime = sscanf(sol.(outputNames{i}).FMOruntime, 'FMO Runtime = %f sec');
    sampRuntime = sol.(outputNames{i}).sampRuntime;
    fprintf('%s & %.4f & %.4f \\\\ \\hline\n', outputNames{i}, round(sampRuntime, 4), round(FMOruntime, 4));
end

rowNames = {'& Avg$_\cH$', '& Max$_\cH$', '& Avg$_\cT$', ...
            '& Min$_\cT$', '& Max$_\cT$', '& D$0.5\%$', '& D$95\%$', '& D$99\%$'};

title = 'Metric';
matSummary =zeros(8,length(outputNames));

fprintf('\n\nThe output metrics for LaTeX are:\n')
% For each type of method
for i = 1:length(outputNames)   
    % Table header:    
    title = strcat([title, ' & ' outputNames{i}]);

    % Assign row-values
    v = d_OAR*sol.(outputNames{i}).w';
    matSummary(1,i) = mean(v);
    matSummary(2,i) = max(v);
        
    v = d_target*sol.(outputNames{i}).w';
    matSummary(3,i) = mean(v);
    matSummary(4,i) = min(v);
    matSummary(5,i) = max(v);

    numTvoxels = size(d_target,1);
    v = sort(v,'descend');
    %D0.5% 
    matSummary(6,i) = v(ceil(0.005*numTvoxels)-1);
    %D95% 
    matSummary(7,i) = v(ceil(0.95*numTvoxels)-1);
    %D99% 
    matSummary(8,i) = v(ceil(0.99*numTvoxels)-1);
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









