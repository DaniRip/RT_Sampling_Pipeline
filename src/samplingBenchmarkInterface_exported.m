classdef samplingBenchmarkInterface_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
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

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            constants = [app.BasicdoseCheckBox.Value, app.CVaRCheckBox.Value, app.UpperboundsCheckBox.Value, app.UniformityCheckBox_2.Value];
            objectives = [app.MinhealthyCheckBox.Value, app.MintargetCheckBox.Value, app.UniformityCheckBox.Value];
            layers = app.layersEditField.Value;
            voxelsField = app.voxelsEditField.Value;
            samplingType = app.SamplingTypeDropDown.Value;
            fileSelected = app.PlanFileTextArea.Value;
            sampleinterval = app.SampleintervalEditField.Value;
            clusters = app.numofclustersEditField.Value;

            if sum(constants)==0 || sum(objectives==0)
                app.OutputTextArea.Value = 'Choose a set of objectives and constraints!';
                break
            end
            
            switch samplingType
                case "Downsample"
                    [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, sampling_runtime] = integerdownsample(app.file, str2double(sampleinterval));
                case "K Clustering"
                    [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, sampling_runtime] = kmeansclustering(app.file, str2num(clusters), str2double(sampleinterval));
                case "Manham K Clustering"
                    [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, sampling_runtime] = kmeansclusteringmanham(app.file, str2num(clusters), str2double(sampleinterval));
                case "Fountain Clustering"
                    [d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose, sampling_runtime] = fountainclustering(app.file, str2double(sampleinterval));
            end
            
            %[a, b, c] = run_FMO(samplingType, voxelsField, layers, objectives, d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose);
            %app.OutputTextArea.Value = [a,newline,b,newline,"Sampling runtime = "+num2str(sampling_runtime)+" sec",newline,c];
            
            [a, b, c, d] = run_FMO(samplingType, voxelsField, layers, objectives, d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose);
            %import run_FMO matrix output txt file = d
            [avg_deviation, max_overage, voxels_underdosed] = metrics(app.file, d);
            app.OutputTextArea.Value = ["Average deviation = "+num2str(avg_deviation), newline, "Max overage = "+num2str(max_overage), newline, "Voxels underdosed = "+num2str(voxels_underdosed)];
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

        % Callback function
        function SamplingTypeChange(app, event)
            value = app.SamplingTypeDropDown.Value;
             switch samplingType
                case "Downsample"
                    
                case "K Clustering"

                case "Manham K Clustering"

                case "Fountain Clustering"

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
            app.SamplingTypeDropDown.Items = {'Downsample', 'K Clustering', 'Manham K Clustering', 'Fountain Clustering', 'All'};
            app.SamplingTypeDropDown.Position = [116 314 118 22];
            app.SamplingTypeDropDown.Value = 'Downsample';

            % Create UniformityCheckBox
            app.UniformityCheckBox = uicheckbox(app.UIFigure);
            app.UniformityCheckBox.Text = 'Uniformity';
            app.UniformityCheckBox.Position = [287 251 76 22];

            % Create MintargetCheckBox
            app.MintargetCheckBox = uicheckbox(app.UIFigure);
            app.MintargetCheckBox.Text = 'Min. target';
            app.MintargetCheckBox.Position = [287 272 79 22];

            % Create MinhealthyCheckBox
            app.MinhealthyCheckBox = uicheckbox(app.UIFigure);
            app.MinhealthyCheckBox.Text = 'Min. healthy';
            app.MinhealthyCheckBox.Position = [287 293 87 22];

            % Create FMOObjectivesLabel
            app.FMOObjectivesLabel = uilabel(app.UIFigure);
            app.FMOObjectivesLabel.Position = [287 314 91 22];
            app.FMOObjectivesLabel.Text = 'FMO Objectives';

            % Create UniformityCheckBox_2
            app.UniformityCheckBox_2 = uicheckbox(app.UIFigure);
            app.UniformityCheckBox_2.Text = 'Uniformity';
            app.UniformityCheckBox_2.Position = [423 230 76 22];

            % Create UpperboundsCheckBox
            app.UpperboundsCheckBox = uicheckbox(app.UIFigure);
            app.UpperboundsCheckBox.Text = 'Upper bounds';
            app.UpperboundsCheckBox.Position = [423 251 97 22];

            % Create CVaRCheckBox
            app.CVaRCheckBox = uicheckbox(app.UIFigure);
            app.CVaRCheckBox.Text = 'CVaR';
            app.CVaRCheckBox.Position = [423 272 53 22];

            % Create BasicdoseCheckBox
            app.BasicdoseCheckBox = uicheckbox(app.UIFigure);
            app.BasicdoseCheckBox.Text = 'Basic dose';
            app.BasicdoseCheckBox.Position = [423 293 81 22];

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
            app.layersEditFieldLabel.Position = [51 217 48 22];
            app.layersEditFieldLabel.Text = '# layers';

            % Create layersEditField
            app.layersEditField = uieditfield(app.UIFigure, 'text');
            app.layersEditField.Editable = 'off';
            app.layersEditField.Position = [116 217 118 22];

            % Create voxelsEditFieldLabel
            app.voxelsEditFieldLabel = uilabel(app.UIFigure);
            app.voxelsEditFieldLabel.HorizontalAlignment = 'right';
            app.voxelsEditFieldLabel.Position = [48 186 53 22];
            app.voxelsEditFieldLabel.Text = '% voxels';

            % Create voxelsEditField
            app.voxelsEditField = uieditfield(app.UIFigure, 'text');
            app.voxelsEditField.Position = [116 186 118 22];

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
            app.SampleintervalEditFieldLabel.Position = [11 280 88 22];
            app.SampleintervalEditFieldLabel.Text = 'Sample interval';

            % Create SampleintervalEditField
            app.SampleintervalEditField = uieditfield(app.UIFigure, 'text');
            app.SampleintervalEditField.Position = [116 280 118 22];

            % Create ofclustersEditFieldLabel
            app.ofclustersEditFieldLabel = uilabel(app.UIFigure);
            app.ofclustersEditFieldLabel.HorizontalAlignment = 'right';
            app.ofclustersEditFieldLabel.Position = [30 248 70 22];
            app.ofclustersEditFieldLabel.Text = '# of clusters';

            % Create numofclustersEditField
            app.numofclustersEditField = uieditfield(app.UIFigure, 'text');
            app.numofclustersEditField.Editable = 'off';
            app.numofclustersEditField.Position = [116 248 118 22];

            % Create VisualizeButton
            app.VisualizeButton = uibutton(app.UIFigure, 'push');
            app.VisualizeButton.ButtonPushedFcn = createCallbackFcn(app, @VisualizeButtonPushed, true);
            app.VisualizeButton.Position = [376 404 100 22];
            app.VisualizeButton.Text = 'Visualize';

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