classdef samplingBenchmarkInterface_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        VisualizeSamplingCheckBox     matlab.ui.control.CheckBox
        VisualizeButton               matlab.ui.control.Button
        numofclustersEditField        matlab.ui.control.EditField
        ofclustersEditFieldLabel      matlab.ui.control.Label
        SampleintervalEditField       matlab.ui.control.EditField
        SampleintervalEditFieldLabel  matlab.ui.control.Label
        PlanFileTextArea              matlab.ui.control.TextArea
        PlanFileLabel                 matlab.ui.control.Label
        BrowseButton                  matlab.ui.control.Button
        OutputTextArea                matlab.ui.control.TextArea
        OutputLabel                   matlab.ui.control.Label
        voxelsEditField               matlab.ui.control.EditField
        voxelsEditFieldLabel          matlab.ui.control.Label
        layersEditField               matlab.ui.control.EditField
        layersEditFieldLabel          matlab.ui.control.Label
        RunButton                     matlab.ui.control.Button
        FMOConstraintsLabel           matlab.ui.control.Label
        BasicdoseCheckBox             matlab.ui.control.CheckBox
        CVaRCheckBox                  matlab.ui.control.CheckBox
        UpperboundsCheckBox           matlab.ui.control.CheckBox
        UniformityCheckBox_2          matlab.ui.control.CheckBox
        FMOObjectivesLabel            matlab.ui.control.Label
        MinhealthyCheckBox            matlab.ui.control.CheckBox
        MintargetCheckBox             matlab.ui.control.CheckBox
        UniformityCheckBox            matlab.ui.control.CheckBox
        SamplingTypeDropDown          matlab.ui.control.DropDown
        SamplingTypeDropDownLabel     matlab.ui.control.Label
        PlanSelectionPanel            matlab.ui.container.Panel
        PlaninfoLabel                 matlab.ui.control.Label
    end

    properties (Access = private)
        file % Description
        inputData  % Property to store input data
    end
    
    methods (Access = private)
        
        function [sol,nLayers] = runPlans(app, sol, type, optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetVCoords, p, k, n, l)
        % A separate function to actually run the plans
        % Note: the returned target_voxels is relative to the input tumor dij, not the absolute locations in the overall structure!

            visSamp = app.VisualizeSamplingCheckBox.Value;
            nLayers = [];

            switch type
               case "Downsample"
                    [d_target, target_voxels, sol.sampRuntime] = integerdownsample(d_target, n);
                    if (visSamp); plot_voxels(targetVCoords(target_voxels,:), type); end
                case "k-Means (coord)"
                    [d_target, target_voxels, targetVCoords, sol.clusters, sol.sampRuntime] = kmeans_coordinates(targetVCoords, d_target, k);
                    if (visSamp); plot_voxels(targetVCoords, type); end
                case "k-Means (Dij)"
                    if k < 20 % For a small enough number of clusters, use the regular algorithm
                        [d_target, target_voxels, targetVCoords, sol.clusters, sol.sampRuntime] = kmeans_Dij(targetVCoords, d_target, k);
                    else % However, if there are more, shortcut it
                        [d_target, target_voxels, targetVCoords, sol.clusters, sol.sampRuntime] = kmeans_DijHeuristicSpeedup(targetVCoords, d_target, k);
                    end
                    if (visSamp); plot_voxels(targetVCoords, type); end
                case "k-Neighbour"
                    [d_target, target_voxels, targetVCoords, sol.clusters, sol.sampRuntime] = kmeans_DijNeighbour(targetVCoords, d_target, k);
                    if (visSamp); plot_voxels(targetVCoords, type); end
                case "Beamlet-Based (max)"
                    [d_target, target_voxels, sol.clusters, sol.sampRuntime] = fountainclustering(d_target, k, 'max');
                    if (visSamp); plot_voxels(targetVCoords(target_voxels,:),type); end
                case "Beamlet-Based (med)"
                    [d_target, target_voxels, sol.clusters, sol.sampRuntime] = fountainclustering(d_target, k, 'med');
                    if (visSamp); plot_voxels(targetVCoords(target_voxels,:), type); end
                case "Layered" 
                    [d_target, target_voxels, sol.sampRuntime] = surfaceSampler_Gatik_v3(d_target, targetVCoords, l, k);
                    if (visSamp); plot_voxels(targetVCoords(target_voxels,:), type); end
                    disp(strcat(['Number of Voxels Selected: ',num2str(numel(target_voxels))]));
                    nLayers = numel(target_voxels);
            end
            
            nTargetV = numel(target_voxels);
            sol.target_voxels = target_voxels;
            [sol.status, sol.obj, sol.FMOruntime, sol.w] = run_FMO(optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose);
            sol.inputs = [n, k, p, l];
            app.OutputTextArea.Value = strcat([app.OutputTextArea.Value; type; num2str(sol.status); num2str(sol.obj); sol.FMOruntime; strcat(['subtime: ',num2str(sol.sampRuntime)])]);
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
    
            constraints = [app.BasicdoseCheckBox.Value, app.CVaRCheckBox.Value, app.UpperboundsCheckBox.Value, app.UniformityCheckBox_2.Value];
            optObjs = [app.MinhealthyCheckBox.Value, app.MintargetCheckBox.Value, app.UniformityCheckBox.Value];
            samplingType = app.SamplingTypeDropDown.Value;
            layers = app.layersEditField.Value;
            voxelsPercent = app.voxelsEditField.Value;
            nSample = app.SampleintervalEditField.Value;
            clusters = app.numofclustersEditField.Value;
            if isfield(app.inputData, 'voxelDim')
                voxelDim = app.inputData.voxelDim; 
            else
                voxelDim = [];
            end

            % Don't send things through if nothing is selected
            if sum(constraints)==0 || sum(optObjs)==0
                app.OutputTextArea.Value = 'Choose a set of objectives and constraints!';
                return
            end
            
            % Yell at users if no data has been uploaded!
            if isempty(app.inputData)
                app.OutputTextArea.Value = 'There''s nothing to run, upload some data!';
                return
            end

            if (isempty(layers) && isempty(nSample) && isempty(clusters) && isempty(voxelsPercent))
                app.OutputTextArea.Value = 'Select a sampling rate (put 1 for the complete plan)';
                return
            end

            app.OutputTextArea.Value = 'Running';

            % Set original inputs
            nBeamlets = size(app.inputData.Dij,2);
            pdose = app.inputData.targetDose; %Prescribed dose
            target_voxels = app.inputData.structVoxels.target;
            nTargetV = numel(app.inputData.structVoxels.target);
            nOARV = numel(app.inputData.structVoxels.OAR);
            d_target = app.inputData.Dij(app.inputData.structVoxels.target(:), :);
            d_OAR = app.inputData.Dij(app.inputData.structVoxels.OAR(:), :);
            targetCoords = app.inputData.voxelIndices.target(target_voxels,:);
            if ~isempty(clusters); k = str2double(clusters); else k =[]; end
            if ~isempty(nSample); n = str2double(nSample); else n =[];  end
            if ~isempty(voxelsPercent); p = str2double(voxelsPercent); else p =[];  end
            if ~isempty(layers); l = str2double(layers); else l =[];  end
            
            % If we're running the full FMO, there's no reason for any sampling 
            if ((str2double(nSample)==1 || str2double(clusters) == nTargetV || str2double(voxelsPercent) ==100))
                disp('Running the full model')
                [sol.Full.status, sol.Full.obj, sol.Full.FMOruntime, sol.Full.w] = run_FMO(optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose);
                sol.Full.sampRuntime = 0;
            else % Otherwise, we need to see what type of sampling we're doing  
                fprintf('Running in the mode %s\n', samplingType);
                switch samplingType
                    case "Downsample"
                        sol.Integer =[];
                        [sol.Integer, ~] = runPlans(app, sol.Integer, "Downsample", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                    case "k-Means (coord)"
                        sol.kmeanC =[];
                        [sol.kmeanC, ~] = runPlans(app, sol.kmeanC, "k-Means (coord)", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                    case "k-Means (Dij)"
                        sol.kmeanD = [];
                        [sol.kmeanD, ~] = runPlans(app, sol.kmeanD, "k-Means (Dij)", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                    case "k-Neighbour"
                        sol.kmeanNN = [];
                        [sol.kmeanNN, ~] = runPlans(app, sol.kmeanNN, "k-Neighbour", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                    case "Beamlet-Based (max)"
                        sol.BBmax = [];
                        [sol.BBmax, ~] = runPlans(app, sol.BBmax, "Beamlet-Based (max)", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                    case "Beamlet-Based (med)"
                        sol.BBmed = [];
                        [sol.BBmed, ~] = runPlans(app, sol.BBmed, "Beamlet-Based (med)", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                    case "Layered" 
                        sol.Layered = [];
                        [sol.Layered, ~] = runPlans(app, sol.Layered, "Layered", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                    case "All"
                        sol.Integer =[]; sol.kmeanC =[];sol.kmeanD = [];sol.kmeanNN = [];sol.BBmax = [];sol.BBmed = [];sol.Layered = [];
                        % If sample-first run conversions for the others
                        if ~isempty(layers) 
                             % Need to run layered approach first if its determining the sampling. 
                            [sol.Layered, ~] = runPlans(app, sol.Layered, "Layered", optObjs, d_target, d_OAR, target_voxels, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                            [n,k,p] = sampleConversion(0, 1, 0, numel(target_voxels), target_voxels);
                            app.SampleintervalEditField.Value = num2str(n);
                            app.voxelsEditField.Value = num2str(p);
                            app.numofclustersEditField.Value = num2str(k);
                            app.layersEditField.Value = num2str(k);
                        end
                        % Then run all others:  
                        [sol.Integer, ~] = runPlans(app, sol.Integer, "Downsample", optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                        [sol.kmeanC, ~] = runPlans(app, sol.kmeanC, "k-Means (coord)", optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                        [sol.kmeanD, ~] = runPlans(app, sol.kmeanD, "k-Means (Dij)", optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                        [sol.kmeanNN, ~] = runPlans(app, sol.kmeanNN, "k-Neighbour", optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                        [sol.BBmax, ~] = runPlans(app, sol.BBmax, "Beamlet-Based (max)", optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                        [sol.BBmed, ~] = runPlans(app, sol.BBmed, "Beamlet-Based (med)", optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);
                        if isempty(layers) % Keep the sequence for the sake of our tables!
                            [sol.Layered, ~] = runPlans(app, sol.Layered, "Layered", optObjs, d_target, d_OAR, nTargetV, nOARV, nBeamlets, pdose, targetCoords, p, k, n, l);        
                        end
                end
            end

            [file, path] = uiputfile('*.mat', 'Save As');
            % Check if the user did not cancel the dialog
            if ischar(file)
                fullpath = fullfile(path, file);
                % Save the structure variable 'sol' to the specified file
                save(fullpath, 'sol');
                disp(['Saved sol to ', fullpath]);
            end

            % Based on the solution, lets get all the outputs we need
            printMetrics(sol, d_target, d_OAR, voxelDim);
        end

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)
            %update file selection text field with file path chosen
            %call uigetfile here
            [app.file, path] = uigetfile('*.mat');
            app.file = fullfile(path,app.file);
            app.PlanFileTextArea.Value = app.file;
            % Load the input as guidata
            % I wanted to do this directly with input = load(... this might
            % be a 2022a-specific bug. Check versions
            tempInput = load(app.file);
            % Get the field names of the temporary structure
            fieldNames = fieldnames(tempInput);
            input = tempInput.(fieldNames{1});
            app.inputData = input;
            app.PlaninfoLabel.Text = strcat(['Plan info: ', num2str(size(input.Dij,1)), ' Voxels; ', num2str(numel(input.structVoxels.target)), ' target, ', num2str(numel(input.structVoxels.OAR)), ' OAR. ',num2str(size(input.Dij,2)), ' Beamlets.']);
        end

        % Button pushed function: VisualizeButton
        function VisualizeButtonPushed(app, event)
            visualizeStructures(app.inputData.voxelIndices);
        end

        % Value changed function: SamplingTypeDropDown
        function SamplingTypeChange(app, event)
            samplingType = app.SamplingTypeDropDown.Value;
            switch samplingType
                case "Layered"
                    app.SampleintervalEditField.Editable='on';
                    app.numofclustersEditField.Editable = 'on';
                    app.voxelsEditField.Editable = 'on';
                    app.layersEditField.Editable='on';    
                case "All"
                    app.SampleintervalEditField.Editable='on';
                    app.numofclustersEditField.Editable = 'on';
                    app.voxelsEditField.Editable = 'on';
                    app.layersEditField.Editable='on';
                    app.VisualizeSamplingCheckBox.Value = false;
                otherwise
                    app.SampleintervalEditField.Editable='on';
                    app.numofclustersEditField.Editable = 'on';
                    app.voxelsEditField.Editable = 'on';
                    app.layersEditField.Editable='off';
                    app.layersEditField.Value='';
             end
        end

        % Value changed function: SampleintervalEditField, 
        % ...and 3 other components
        function UpdateValuesSamp(app, event)
            if ~isempty(app.inputData)
                num_target_voxels = numel(app.inputData.structVoxels.target);
                nToAll = 0;
                kToAll=0;
                pToAll=0;
                layersChosen = false;
                if event.Source==app.SampleintervalEditField % If sample field
                    value = app.SampleintervalEditField.Value;
                    numericValue = str2double(value);
                    nToAll = 1;
                elseif event.Source==app.voxelsEditField % If sample field
                    value = app.voxelsEditField.Value;
                    numericValue = str2double(value);
                    pToAll = 1;
                elseif event.Source==app.numofclustersEditField % If sample field
                    value = app.numofclustersEditField.Value;
                    numericValue = str2double(value);
                    kToAll = 1;
                elseif event.Source == app.layersEditField % layers incompatible with all others
                    value = app.layersEditField.Value;
                    numericValue = str2double(value);
                    app.SampleintervalEditField.Value = '';
                    app.voxelsEditField.Value = '';
                    app.numofclustersEditField.Value = '';
                    layersChosen = true;
                end
                if isnan(numericValue)
                    app.TextField.Value = '';  % Reset the text field
                    app.OutputTextArea.Value = 'Please enter a legal value';
                else
                    if ((~pToAll && numericValue<1) || (numericValue<0)  || (pToAll && numericValue>100) || numericValue> num_target_voxels)
                        app.OutputTextArea.Value = 'Sampling rate either too high or too low';
                        app.OutputTextArea.Value = 'Please enter a legal value';
                    elseif (~layersChosen)    
                        [n,k,p] = sampleConversion(nToAll, kToAll, pToAll, numericValue, num_target_voxels);
                        app.SampleintervalEditField.Value = num2str(n);
                        app.voxelsEditField.Value = num2str(p);
                        app.numofclustersEditField.Value = num2str(k);
                        app.OutputTextArea.Value = '';
                        app.layersEditField.Value = '';
                    end
                end
            else
                app.SampleintervalEditField.Value='';
                app.OutputTextArea.Value = 'Please upload data first';
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 610 495];
            app.UIFigure.Name = 'MATLAB App';

            % Create PlanSelectionPanel
            app.PlanSelectionPanel = uipanel(app.UIFigure);
            app.PlanSelectionPanel.Title = 'Plan Selection';
            app.PlanSelectionPanel.Position = [15 354 583 129];

            % Create PlaninfoLabel
            app.PlaninfoLabel = uilabel(app.PlanSelectionPanel);
            app.PlaninfoLabel.Position = [9 9 565 22];
            app.PlaninfoLabel.Text = 'Plan info:';

            % Create SamplingTypeDropDownLabel
            app.SamplingTypeDropDownLabel = uilabel(app.UIFigure);
            app.SamplingTypeDropDownLabel.HorizontalAlignment = 'right';
            app.SamplingTypeDropDownLabel.Position = [15 314 85 22];
            app.SamplingTypeDropDownLabel.Text = 'Sampling Type';

            % Create SamplingTypeDropDown
            app.SamplingTypeDropDown = uidropdown(app.UIFigure);
            app.SamplingTypeDropDown.Items = {'Downsample', 'k-Means (coord)', 'k-Means (Dij)', 'k-Neighbour', 'Beamlet-Based (max)', 'Beamlet-Based (med)', 'Layered', 'All'};
            app.SamplingTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SamplingTypeChange, true);
            app.SamplingTypeDropDown.Position = [116 314 118 22];
            app.SamplingTypeDropDown.Value = 'Downsample';

            % Create UniformityCheckBox
            app.UniformityCheckBox = uicheckbox(app.UIFigure);
            app.UniformityCheckBox.Enable = 'off';
            app.UniformityCheckBox.Text = 'Uniformity';
            app.UniformityCheckBox.Position = [287 251 76 22];

            % Create MintargetCheckBox
            app.MintargetCheckBox = uicheckbox(app.UIFigure);
            app.MintargetCheckBox.Text = 'Min. target';
            app.MintargetCheckBox.Position = [287 272 79 22];
            app.MintargetCheckBox.Value = true;

            % Create MinhealthyCheckBox
            app.MinhealthyCheckBox = uicheckbox(app.UIFigure);
            app.MinhealthyCheckBox.Text = 'Min. healthy';
            app.MinhealthyCheckBox.Position = [287 293 87 22];
            app.MinhealthyCheckBox.Value = true;

            % Create FMOObjectivesLabel
            app.FMOObjectivesLabel = uilabel(app.UIFigure);
            app.FMOObjectivesLabel.Position = [287 314 91 22];
            app.FMOObjectivesLabel.Text = 'FMO Objectives';

            % Create UniformityCheckBox_2
            app.UniformityCheckBox_2 = uicheckbox(app.UIFigure);
            app.UniformityCheckBox_2.Enable = 'off';
            app.UniformityCheckBox_2.Text = 'Uniformity';
            app.UniformityCheckBox_2.Position = [423 230 76 22];

            % Create UpperboundsCheckBox
            app.UpperboundsCheckBox = uicheckbox(app.UIFigure);
            app.UpperboundsCheckBox.Enable = 'off';
            app.UpperboundsCheckBox.Text = 'Upper bounds';
            app.UpperboundsCheckBox.Position = [423 251 97 22];

            % Create CVaRCheckBox
            app.CVaRCheckBox = uicheckbox(app.UIFigure);
            app.CVaRCheckBox.Enable = 'off';
            app.CVaRCheckBox.Text = 'CVaR';
            app.CVaRCheckBox.Position = [423 272 53 22];

            % Create BasicdoseCheckBox
            app.BasicdoseCheckBox = uicheckbox(app.UIFigure);
            app.BasicdoseCheckBox.Text = 'Basic dose';
            app.BasicdoseCheckBox.Position = [423 293 81 22];
            app.BasicdoseCheckBox.Value = true;

            % Create FMOConstraintsLabel
            app.FMOConstraintsLabel = uilabel(app.UIFigure);
            app.FMOConstraintsLabel.Position = [423 314 96 22];
            app.FMOConstraintsLabel.Text = 'FMO Constraints';

            % Create RunButton
            app.RunButton = uibutton(app.UIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Position = [506 139 87 22];
            app.RunButton.Text = 'Run';

            % Create layersEditFieldLabel
            app.layersEditFieldLabel = uilabel(app.UIFigure);
            app.layersEditFieldLabel.HorizontalAlignment = 'right';
            app.layersEditFieldLabel.Position = [51 220 48 22];
            app.layersEditFieldLabel.Text = '# layers';

            % Create layersEditField
            app.layersEditField = uieditfield(app.UIFigure, 'text');
            app.layersEditField.ValueChangedFcn = createCallbackFcn(app, @UpdateValuesSamp, true);
            app.layersEditField.Editable = 'off';
            app.layersEditField.Position = [116 220 118 22];

            % Create voxelsEditFieldLabel
            app.voxelsEditFieldLabel = uilabel(app.UIFigure);
            app.voxelsEditFieldLabel.HorizontalAlignment = 'right';
            app.voxelsEditFieldLabel.Position = [48 190 53 22];
            app.voxelsEditFieldLabel.Text = '% voxels';

            % Create voxelsEditField
            app.voxelsEditField = uieditfield(app.UIFigure, 'text');
            app.voxelsEditField.ValueChangedFcn = createCallbackFcn(app, @UpdateValuesSamp, true);
            app.voxelsEditField.Position = [116 190 118 22];

            % Create OutputLabel
            app.OutputLabel = uilabel(app.UIFigure);
            app.OutputLabel.Position = [16 118 45 22];
            app.OutputLabel.Text = 'Output:';

            % Create OutputTextArea
            app.OutputTextArea = uitextarea(app.UIFigure);
            app.OutputTextArea.Position = [15 21 583 92];

            % Create BrowseButton
            app.BrowseButton = uibutton(app.UIFigure, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.Position = [488 404 100 22];
            app.BrowseButton.Text = 'Browse';

            % Create PlanFileLabel
            app.PlanFileLabel = uilabel(app.UIFigure);
            app.PlanFileLabel.HorizontalAlignment = 'right';
            app.PlanFileLabel.Position = [18 433 56 22];
            app.PlanFileLabel.Text = 'Plan File:';

            % Create PlanFileTextArea
            app.PlanFileTextArea = uitextarea(app.UIFigure);
            app.PlanFileTextArea.Position = [86 434 501 20];

            % Create SampleintervalEditFieldLabel
            app.SampleintervalEditFieldLabel = uilabel(app.UIFigure);
            app.SampleintervalEditFieldLabel.HorizontalAlignment = 'right';
            app.SampleintervalEditFieldLabel.Position = [11 281 88 22];
            app.SampleintervalEditFieldLabel.Text = 'Sample interval';

            % Create SampleintervalEditField
            app.SampleintervalEditField = uieditfield(app.UIFigure, 'text');
            app.SampleintervalEditField.ValueChangedFcn = createCallbackFcn(app, @UpdateValuesSamp, true);
            app.SampleintervalEditField.Position = [116 281 118 22];

            % Create ofclustersEditFieldLabel
            app.ofclustersEditFieldLabel = uilabel(app.UIFigure);
            app.ofclustersEditFieldLabel.HorizontalAlignment = 'right';
            app.ofclustersEditFieldLabel.Position = [30 250 70 22];
            app.ofclustersEditFieldLabel.Text = '# of clusters';

            % Create numofclustersEditField
            app.numofclustersEditField = uieditfield(app.UIFigure, 'text');
            app.numofclustersEditField.ValueChangedFcn = createCallbackFcn(app, @UpdateValuesSamp, true);
            app.numofclustersEditField.Position = [116 250 118 22];

            % Create VisualizeButton
            app.VisualizeButton = uibutton(app.UIFigure, 'push');
            app.VisualizeButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizeButtonPushed, true);
            app.VisualizeButton.Position = [376 404 100 22];
            app.VisualizeButton.Text = 'Visualize';

            % Create VisualizeSamplingCheckBox
            app.VisualizeSamplingCheckBox = uicheckbox(app.UIFigure);
            app.VisualizeSamplingCheckBox.Text = 'Visualize Sampling';
            app.VisualizeSamplingCheckBox.Position = [260 190 133 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = samplingBenchmarkInterface_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end