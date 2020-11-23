classdef Transistor < Component
    
    properties (SetAccess = protected)
        %%%% Properties %%%%
        type = '';
        types
        material = '';
        materials
        MLver
        
        %%%% Inherited Publicly Viewable Properties %%%%
        % componentType
            % 'Transistor'
        % partNumber
            % String, the part number. Must be passed in constructor.
        % tableFields 
            % Cell array of strings, each is the name for a table property which may contain Min, Typ, and Max values.
            % See constructor for values.
        % tableUnits
            % Cell array of strings, each is the unit for the respective property from tableFields.
            % See constructor for values.
        % tableParameters
            % Structure containing information on table properties.
            % Each element corresponds to a tableField and each of those may contain a Min, Typ, and Max field.
        % graphParameters
            % Structure containing information on graph properties.
            % Each elelement corresponds to a graph parameter and has the following fields:
                % Data - 2 column numeric array, each row is an (x,y) coordinate.
                % Title - String, shown as title() in plot.
                % xLabel - String, shown as xlabel() in plot. Indicate units here!
                % yLabel - String, shown as yLabel() in plot. Indicate units here!
                % Log - 1x2 boolean array. If Log(1), x-axis will use log scale. Log(2) for y-axis.
        
        %%%% Inherited Protected Properites %%%%     
        % SIkeys
            % Cell array of strings - Keys for SI prefixes dict ordered properly
        % SIprefixes
            % Dictionary string->double - For graph parameter data scaling
    end
    
    
    % Public Methods
    methods     
        % Constructor - Must pass partNumber string
        function obj = Transistor(partNumber, noGUIs)   
            obj.componentType = 'Transistor';
            
            [obj.types,obj.materials,...
            obj.tableFields,obj.tableUnits,obj.tableDefaultMultipliers,...
            obj.graphFields,obj.graphUnits,obj.graphDefaultMultipliers,...
            obj.conditionFields,obj.conditionUnits,obj.conditionDefaultMultipliers] = Transistor.loadMetaData();
            
            % Given device number
            if nargin >= 1
                devices = obj.listDevices();
                obj.partNumber = Transistor.checkDuplicateDevices(strip(partNumber), devices);
                obj.loadData();
            end
            
            if nargin < 2 || ~noGUIs
                obj.addTableParametersGUI();
                obj.addAdditionalTableParametersGUI();
                obj.addGraphParametersGUI();
            end
            
            if ~isempty(obj.partNumber) && ~isempty(obj.type) && ~isempty(obj.material)
                obj.saveData();
            end
            
            MLver = version;
            [si, ei] = regexp(MLver, '\(.*\)');
            obj.MLver = {str2double(MLver(si+2:ei-2)), MLver(ei-1)};
        end
   
        % GUI for entering datasheet table parameters information
        function addTableParametersGUI(obj)
            % Create 3 row grid layout, rescales when window size is changed
            UIWidth = 550;
            UIHeight = 800;
            fig = uifigure('Position',[100 100 UIWidth UIHeight],'Visible', 'off'); % Do not display until all information is added.
            fig.Name = 'Transistor Datasheet Table Parameters';
            grid = uigridlayout(fig,'Scrollable','on');
            grid.RowHeight = {150,'1x',40};
            grid.ColumnWidth = {'1x'};
                            
            % Create Title
            titleGrid = uigridlayout(grid,[4,2]);
            titleGrid.ColumnWidth = {300,'1x'};
            titleLabel = uilabel(titleGrid);
            titleLabel.Text = 'Enter information for transistor';
            titleLabel.FontSize = 24;
            titleLabel.FontWeight = 'bold';
            titleLabel.HorizontalAlignment = 'center';
            titleLabel.Layout.Column = [1,2];

            % partNumber, type, material
            partNumberLabel = uilabel(titleGrid);
            partNumberLabel.Text = 'Manufacturer Part Number (Ex. C2M0040120D)';
            partNumberLabel.Layout.Column = 1;
            partNumberBox = uitextarea(titleGrid);
            partNumberBox.Layout.Column = 2;
            if ~isempty(obj.partNumber)
                partNumberBox.Value = obj.partNumber;
            end
            
            typeLabel = uilabel(titleGrid);
            typeLabel.Text = 'Transistor Type';
            typeLabel.Layout.Column = 1;
            typeBox = uidropdown(titleGrid,'Items',[{'-Select Type-'}; (obj.types)]);
            typeBox.Layout.Column = 2;
            if ~isempty(obj.type)
                typeBox.Value = obj.type;
            end
            
            materialLabel = uilabel(titleGrid);
            materialLabel.Text = 'Transistor Material';
            materialLabel.Layout.Column = 1;
            materialBox = uidropdown(titleGrid,'Items',[{'-Select Material-'}; (obj.materials)]);
            materialBox.Layout.Column = 2;
            if ~isempty(obj.material)
                materialBox.Value = obj.material;
            end
            
            
            % Table Parameters
            tableParamsGrid = uigridlayout(grid,'Scrollable','on'); 
            gridRows = cell(1,1+numel(obj.tableFields));
            if obj.MLver > 2019
                tableParamsGrid.ColumnWidth = {'fit','1x','1x','1x','fit'};
                gridRows(1) = {'fit'};
            else
                tableParamsGrid.ColumnWidth = {'1x','1x','1x','1x','1x'};
                gridRows(1) = {'1x'};
            end
            gridRows(2:end) = {20};
            tableParamsGrid.RowHeight = gridRows;
            
            minLabel = uilabel(tableParamsGrid);
            minLabel.Text = 'Min.';
            minLabel.HorizontalAlignment = 'center';
            minLabel.FontWeight = 'bold';
            minLabel.FontSize = 14;
            minLabel.Layout.Column = 2;
            
            typLabel = uilabel(tableParamsGrid);
            typLabel.Text = 'Typ.';
            typLabel.HorizontalAlignment = 'center';
            typLabel.FontWeight = 'bold';
            typLabel.FontSize = 14;
            typLabel.Layout.Column = 3;
            
            maxLabel = uilabel(tableParamsGrid);
            maxLabel.Text = 'Max.';
            maxLabel.HorizontalAlignment = 'center';
            maxLabel.FontWeight = 'bold';
            maxLabel.FontSize = 14;
            maxLabel.Layout.Column = 4;
            
            for i = 1:numel(obj.tableFields)
                field = obj.tableFields{i};
                
                rowLabel = uilabel(tableParamsGrid);
                rowLabel.Text = obj.tableFields{i};
                rowLabel.FontSize = 14;
                rowLabel.HorizontalAlignment = 'center';
                rowLabel.Layout.Column = 1;
                rowLabel.Layout.Row = 1+i;
                
                minBox = uitextarea(tableParamsGrid);
                minBox.Layout.Column = 2;
                minBox.Layout.Row = 1+i;                
                if isfield(obj.tableParameters,field) && isfield(obj.tableParameters.(field),'Min') && ~isempty(obj.tableParameters.(field).Min)
                    minBox.Value = num2str(obj.tableParameters.(field).Min);
                end
                minBoxes(i) = minBox;
                
                typBox = uitextarea(tableParamsGrid);
                typBox.Layout.Column = 3;
                typBox.Layout.Row = 1+i;
                if isfield(obj.tableParameters,field) && isfield(obj.tableParameters.(field),'Typ') && ~isempty(obj.tableParameters.(field).Typ)
                    typBox.Value = num2str(obj.tableParameters.(field).Typ);
                end                
                typBoxes(i) = typBox;
                
                maxBox = uitextarea(tableParamsGrid);
                maxBox.Layout.Column = 4;
                maxBox.Layout.Row = 1+i;
                if isfield(obj.tableParameters,field) && isfield(obj.tableParameters.(field),'Max') && ~isempty(obj.tableParameters.(field).Max)
                    maxBox.Value = num2str(obj.tableParameters.(field).Max);
                end                      
                maxBoxes(i) = maxBox;       
                
                unitLabel = uilabel(tableParamsGrid);
                unit = obj.tableUnits{i};
                multiplier = obj.tableDefaultMultipliers{i};
                if strcmp(multiplier, '1')
                    multiplier = '';
                end
                unitLabel.Text = [multiplier unit];
                unitLabel.FontSize = 14;
            end
            
            % Save button
            uibutton(grid,'Text','Save Data','FontSize',16,'FontWeight','bold','ButtonPushedFcn',...
                @(btn,event) obj.saveTableBtnPushed(partNumberBox,typeBox,materialBox,minBoxes,typBoxes,maxBoxes));
            
            % Make datasheet visible and pause further code execution until it is closed
            fig.Visible = 'on';
            uiwait(fig)
        end
      
        % GUI for entering additional table parameters information
        function addAdditionalTableParametersGUI(obj)
            % Create 3 row grid layout, rescales when window size is changed
            UIWidth = 1300;
            UIHeight = 600;
            fig = uifigure('Position',[500 200 UIWidth UIHeight],'Visible', 'off'); % Do not display until all information is added.
            fig.Name = ['Transistor Additional Table Parameters: ' obj.partNumber];
            grid = uigridlayout(fig,'Scrollable','on');
            grid.RowHeight = {150,'1x',50};
            grid.ColumnWidth = {'1x'};
            
            topGrid = uigridlayout(grid,[1,2]);
            topGrid.ColumnWidth = {'2x','3x'};
            
            topLeftGrid = uigridlayout(topGrid);
            topLeftGrid.RowHeight = {'2x','1x'};
            topLeftGrid.ColumnWidth = {'2x','1x',10,'1x',60};
            paramDropdown = uidropdown(topLeftGrid,'Items',['--';strcat(obj.graphFields,{' ('},...
                strrep(obj.graphDefaultMultipliers,'1',''),obj.graphUnits,{')'})]);
            paramDropdown.FontSize = 20;
            typeDropdown = uidropdown(topLeftGrid,'Items',{'--';'Min';'Typ';'Max'});
            typeDropdown.FontSize = 20;
            lbl = uilabel(topLeftGrid,'Text','=');
            lbl.FontSize = 20;
            valueBox = uitextarea(topLeftGrid);
            lbl = uilabel(topLeftGrid,'Text', 'when:');
            lbl.FontSize = 20;
            ExpBox = uicheckbox(topLeftGrid,'text','Experimental Data');
            ExpBox.Layout.Column = [2 4];
            
            topRightGrid = uigridlayout(topGrid);
            topRightGrid.ColumnWidth = {'2x',10,'1x',5,'2x',10,'1x',5,'2x',10,'1x',5};
            topRightGrid.RowHeight = {'1x','1x','1x'};
            for i = 1:6
                conditionDropdowns(i) = uidropdown(topRightGrid,'Items',['--';strcat(obj.conditionFields,{' ('},...
                    strrep(obj.conditionDefaultMultipliers,'1',''),obj.conditionUnits,{')'})]);
                lbl = uilabel(topRightGrid,'Text','=');
                lbl.FontSize = 20;
                conditionBoxes(i) = uitextarea(topRightGrid);
                lbl = uilabel(topRightGrid,'Text',',');
                lbl.FontWeight = 'bold';
                lbl.FontSize = 20;
            end
            extraBox = uitextarea(topRightGrid,'Value','Add additional (non-sortable) conditions if needed.');
            extraBox.Layout.Column = [1 12];
            
            mainText = uitextarea(grid,'Editable','off');
            mainText.FontSize = 14;
            
            % Save button
            saveBtn = uibutton(grid,'Text','Add/Save Row','ButtonPushedFcn',...
                @(btn,event) obj.saveAdditionalBtnPushed(mainText,paramDropdown,typeDropdown,valueBox,conditionDropdowns,conditionBoxes,extraBox,ExpBox));
            saveBtn.FontWeight = 'bold';
            saveBtn.FontSize = 20;
            
            % Make datasheet visible and pause further code execution until it is closed
            obj.loadAdditionalDataGUI(mainText);
            fig.Visible = 'on';
            uiwait(fig)            
        end
        
        % GUI for entering graph parameters information
        function addGraphParametersGUI(obj)
            
            obj.current_graphPage = 1;
            obj.temp_graphParameters = obj.graphParameters;
            % If there are already graph parameters, choose the first as the current parameter
            if numel(fieldnames(obj.temp_graphParameters)) > 0
                fields = fieldnames(obj.temp_graphParameters);
                fields = [fields; 'NEWPARAMETER'];
                obj.current_graphParameter = fields{1};
            else
                fields = {'NEWPARAMETER'};
                obj.current_graphParameter = 'NEWPARAMETER';
            end
            obj.temp_graphParameters.NEWPARAMETER.Data = [];
            obj.temp_graphParameters.NEWPARAMETER.Multiplier = ['1','1'];
            obj.temp_graphParameters.NEWPARAMETER.Log = [0,0];
            obj.temp_graphParameters.NEWPARAMETER.Conditions = struct();
            obj.temp_graphParameters.NEWPARAMETER.ExtraConditions = '';
            
            % Create 3 row grid layout, rescales when window size is changed
            UIWidth = 1600;
            UIHeight = 700;
            fig = uifigure('Position',[100 300 UIWidth UIHeight],'Visible', 'off'); % Do not display until all information is added.
            fig.Name = ['Transistor Graph Parameters: ' obj.partNumber];
            grid = uigridlayout(fig,'Scrollable','on');
            grid.RowHeight = {60,'1x',60};
            grid.ColumnWidth = {'4x','1x'};
            
            % Top Grid
            topGrid = uigridlayout(grid);
            topGrid.Layout.Column = [1 2];
            topGrid.RowHeight = {'1x'};
            topGrid.ColumnWidth = {'1x','1x','1x','1x'};
            
            % Data Grid
            dataGrid = uigridlayout(grid);
            dataGrid.RowHeight = {'1x','1x','1x',60};
            dataGrid.ColumnWidth = {'1x',120,'1x','1x','1x'};
            dataBox = uitextarea(dataGrid);
            dataBox.Layout.Column = 1;
            dataBox.Layout.Row = [1 4];
            
            dataPlot = uiaxes(dataGrid);
            dataPlot.Layout.Column = [3 5];
            dataPlot.Layout.Row = [1 3];
            disableDefaultInteractivity(dataPlot)
            
            y_fieldBox = uidropdown(dataGrid);
            y_fieldBox.FontSize = 16;
            y_fieldBox_items = cell(1,numel(obj.graphFields)+1);
            y_fieldBox_items(1) = {'--Select Y-axis Field--'};
            for i = 1:numel(obj.graphFields)
                y_fieldBox_items(i+1) = {[obj.graphFields{i} ' (' obj.graphUnits{i} ')']};
            end
            y_fieldBox.Items = y_fieldBox_items;
            y_fieldBox.Layout.Row = 1;
            
            y_multiplierBox = uidropdown(dataGrid, 'Items', Transistor.SIkeys);
            y_multiplierBox.FontSize = 16;
            y_multiplierBox.Value = {'1'};
            y_multiplierBox.Layout.Row = 2;
            
            y_logBox = uidropdown(dataGrid, 'Items', [{'Linear'},{'Log'}]);
            y_logBox.FontSize = 16;
            y_logBox.Layout.Row = 3;

            x_fieldBox = uidropdown(dataGrid);
            x_fieldBox.FontSize = 16;
            x_fieldBox_items = cell(1,numel(obj.graphFields)+1);
            x_fieldBox_items(1) = {'--Select X-axis Field--'};
            for i = 1:numel(obj.graphFields)
                x_fieldBox_items(i+1) = {[obj.graphFields{i} ' (' obj.graphUnits{i} ')']};
            end
            x_fieldBox.Items = x_fieldBox_items;
            x_fieldBox.Layout.Row = 4;
            x_fieldBox.Layout.Column = 3;
            
            x_multiplierBox = uidropdown(dataGrid, 'Items', Transistor.SIkeys);
            x_multiplierBox.FontSize = 16;
            x_multiplierBox.Value = {'1'};
            x_multiplierBox.Layout.Row = 4;
            x_multiplierBox.Layout.Column = 4;
            
            x_logBox = uidropdown(dataGrid, 'Items', [{'Linear'},{'Log'}]);
            x_logBox.FontSize = 16;
            x_logBox.Layout.Row = 4;
            x_logBox.Layout.Column = 5;
            
            updateBtn = uibutton(dataGrid,'Text','Update Plot','ButtonPushedFcn',...
                @(btn,event) obj.updateBtnPushed(dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, x_logBox, y_multiplierBox, x_multiplierBox));
            updateBtn.FontWeight = 'bold';
            updateBtn.FontSize = 18;
            updateBtn.BackgroundColor = [.9 .92 .92];
            updateBtn.Layout.Column = 2;
            updateBtn.Layout.Row = 4;  
            
            % Conditions Grid
            conditionsGrid = uigridlayout(grid);
            conditionsGrid.ColumnWidth = {'1x','1x',70};
            if obj.MLver > 2019
                conditionsGrid.RowHeight = {'fit','fit',50,50,50,50,50,'1x'};
            else
                conditionsGrid.RowHeight = {'1x','1x',50,50,50,50,50,'1x'};
            end
            
            conditionsLabel = uilabel(conditionsGrid);
            conditionsLabel.Text = 'Conditions';
            conditionsLabel.FontSize = 24;
            conditionsLabel.HorizontalAlignment = 'Center';
            conditionsLabel.Layout.Column = [1 3];
    
            expCheckbox = uicheckbox(conditionsGrid, 'Text', 'Experimental Data', 'Value', 0);
            expCheckbox.Layout.Column = [1 3];
            
            conditionBoxes = cell(5,2);
            for i = 1:5
                h = uidropdown(conditionsGrid);
                h.Items = [{'--Add Condition--'} ; strcat(obj.conditionFields, {' ('}, strrep(obj.conditionDefaultMultipliers,'1',''), obj.conditionUnits, {')'})];
                h.Layout.Column = [1 2];
                h.FontSize = 16;
                conditionBoxes(i,1) = {h};
                t = uitextarea(conditionsGrid);
                t.FontSize = 16;
                conditionBoxes(i,2) = {t};
            end
            conditionExtras = uitextarea(conditionsGrid, 'Value', 'Add additional (non-sortable) conditions if needed.');
            conditionExtras.Layout.Column = [1 3];
            
            % Dropdown for top grid
            selectDropdown = uidropdown(topGrid, 'ValueChangedFcn', ...
                @(selectDropdown,event) obj.selection(selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, ...
                y_logBox, x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras));
            selectDropdown.Items = strrep(fields,'VERSUS',' vs. ');
            selectDropdown.FontSize = 20;
            selectDropdown.FontWeight = 'bold';            
            selectDropdown.Layout.Column = [1 4];

            % Buttons Grid
            buttonsGrid = uigridlayout(grid,[1 4]);
            buttonsGrid.Layout.Column = [1 2];
            prevBtn = uibutton(buttonsGrid,'Text','Previous Page','ButtonPushedFcn',...
                @(btn,event) obj.prevBtnPushed(selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, ...
                y_logBox, x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras));
            prevBtn.FontSize = 20;
            prevBtn.Layout.Column = 1;
            prevBtn.HorizontalAlignment = 'center';
            
            saveBtn = uibutton(buttonsGrid,'Text','Save Data','ButtonPushedFcn',...
                @(btn,event) obj.saveGraphBtnPushed(selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, ...
                y_logBox, x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras));
            saveBtn.FontSize = 20;
            saveBtn.Layout.Column = [2 3];
            saveBtn.HorizontalAlignment = 'center';
            
            nextBtn = uibutton(buttonsGrid,'Text','Next/Add Page','ButtonPushedFcn',...
                @(btn,event) obj.nextBtnPushed(selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, ...
                y_logBox, x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras));
            nextBtn.FontSize = 20;
            nextBtn.Layout.Column = 4;
            nextBtn.HorizontalAlignment = 'center';
            
            
            % Load current data and plot
            if ~strcmp(obj.current_graphParameter,'NEWPARAMETER')
                obj.loadGraphDataGUI(dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, x_multiplierBox, ...
                    y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);
                obj.updateBtnPushed(dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, ...
                    x_logBox, y_multiplierBox, x_multiplierBox);
            end            
            % Make datasheet visible and pause further code execution until it is closed
            fig.Visible = 'on';
            uiwait(fig)
        end
        
        % Create datasheet UI Figure with all table and graph properties.
        function datasheet(obj)
         
            % Create grid layout, rescales when window size is changed
            UIWidth = 1500;
            UIHeight = 800;
            UIFigure = uifigure('Position',[100 100 UIWidth UIHeight],'Visible', 'off'); % Do not display until all information is added.
            UIFigure.Name = [obj.partNumber ' Datasheet'];
            grid = uigridlayout(UIFigure,[2,3]);
            grid.RowHeight = {100,'1x'};
            grid.ColumnWidth = {'1x','2x','1x'};
            titleLabel = uilabel(grid,'Text',[obj.partNumber newline obj.material ', ' obj.type],...
                'FontSize',24,'FontWeight','bold','HorizontalAlignment','center');
            titleLabel.Layout.Column = [1 3];
            
            % Create Table Parameter Table
            tableParams = uitable(grid);
            tableParams.ColumnName = {'Parameter'; 'Min'; 'Typ'; 'Max'; 'Unit'};
            tableParams.RowName = {};
            tableNames = fieldnames(obj.tableParameters);
            tableParams.Data = cell(length(tableNames),5);
            rowCounter = 1;
            for i = 1:numel(obj.tableFields)
                field = obj.tableFields(i);
                if isfield(obj.tableParameters, field{1})
                    current_param = obj.tableParameters.(field{1});
                    min = [];
                    typ = [];
                    max = [];
                    if isfield(current_param,'Min')
                        min = current_param.Min;
                    end
                    if isfield(current_param,'Typ')
                        typ = current_param.Typ;
                    end
                    if isfield(current_param,'Max')
                        max = current_param.Max;
                    end
                    multiplier = obj.tableDefaultMultipliers{i};
                    if strcmp(multiplier,'1')
                        multiplier = '';
                    end
                    unit = obj.tableUnits{i};
                    unit_str = [multiplier unit];
                    if ~(isempty(min) && isempty(typ) && isempty(max))
                        tableParams.Data(rowCounter,:) = {field{1}, min, typ, max, unit_str};
                        rowCounter = rowCounter + 1;
                    end
                end
            end
            
            % Create Graph Parameter Plots
            graphCount = 0;
            fields = fieldnames(obj.graphParameters);
            if numel(fields) > 0
                for i = 1:numel(fields)
                    field = fields{i};
                    graphCount = graphCount + numel(obj.graphParameters.(field));
                end
                graphGrid = uigridlayout(grid,[graphCount,2],'Scrollable','on');
                graphGrid.ColumnWidth = {'1x',150};
                graphGrid.RowHeight(1:graphCount) = {300};

                for i = 1:numel(fields)
                field = fields{i};
                for j = 1:numel(obj.graphParameters.(field))
                    current_param = obj.graphParameters.(field)(j);
                    Plot = uiaxes(graphGrid);
                    plot(Plot, current_param.Data(:,1), current_param.Data(:,2));
                    split_str = strsplit(field,'VERSUS');
                    xField = split_str{2};
                    yField = split_str{1};   
                    xUnit_index = strcmp(obj.graphFields,xField);
                    xUnit = obj.graphUnits{xUnit_index};
                    yUnit_index = strcmp(obj.graphFields,yField);
                    yUnit = obj.graphUnits{yUnit_index};                
                    disableDefaultInteractivity(Plot)
                    title(Plot, strrep([yField ' vs. ' xField ', Page ' num2str(j)], '_', '\_'));
                    xMultiplier = current_param.Multiplier(1);
                    if strcmp(xMultiplier,'1')
                        xMultiplier = '';
                    end
                    yMultiplier = current_param.Multiplier(2);
                    if strcmp(yMultiplier,'1')
                        yMultiplier = '';
                    end
                    xlabel(Plot, [xField ' (' xMultiplier xUnit ')']);
                    ylabel(Plot, [yField ' (' yMultiplier yUnit ')']);
                    Plot.Box = 'on';
                    if current_param.Log(1)
                        Plot.XScale = 'log';
                    end
                    if current_param.Log(2)
                        Plot.YScale = 'log';
                    end

                    textBox = uitextarea(graphGrid);
                    textStr = ['Conditions:' newline newline];
                    cFields = fieldnames(current_param.Conditions);
                    for j = 1:numel(cFields)
                        unit = obj.conditionUnits(strcmp(obj.conditionFields,cFields{j}));
                        textStr = [textStr newline cFields{j} ': ' num2str(current_param.Conditions.(cFields{j})) ' (' unit{1} ')']; 
                    end
                    if ~strcmp(current_param.ExtraConditions,'Add additional (non-sortable) conditions if needed.')
                        textStr = [textStr newline strjoin(current_param.ExtraConditions, newline)];
                    end
                    textBox.Value = textStr;
                    textBox.Editable = 'off';
                    end
                end
            end
            
            % Create Additional Parameters Text
            additionalParametersBox = uitextarea(grid, 'editable', 'off');
            additionalParametersBox.FontSize = 14;
            additionalParametersBox.Layout.Column = 3;
            obj.loadAdditionalDataGUI(additionalParametersBox)
            
            % Make datasheet visible
            UIFigure.Visible = 'on';
        end
        
        % Public method to load table and graph data from DBs.
        function loadData(obj)
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            try
                obj.load_rldb([path 'TransistorDB\Transistor_StandardTable.db']);
            catch
                warning(['Unable to load RLDB for ' obj.partNumber]);
                hasDBToolboxLicense = license('test','Database_Toolbox');
                v = ver;
                hasDBToolboxInstalled = ~isempty(intersect({v.Name}, 'Database Toolbox'));
                if ~all([hasDBToolboxLicense,hasDBToolboxInstalled])
                    error('Database Toolbox is required');
                end
            end
            try
                obj.additionalTableParameters = obj.load_nrdb([path 'TransistorDB\Transistor_AdditionalTable.json']);
            catch Error
                if ~(contains(Error.message, 'json does not exist') || contains(Error.message, 'not found in non-relational') || ...
                        contains(Error.message, 'JSON syntax error'))
                    rethrow(Error)
                else
                    % Part not found, do nothing
                end
            end
            try
                obj.graphParameters = obj.load_nrdb([path 'TransistorDB\Transistor_Graph.json']);
            catch Error
                if ~(contains(Error.message, 'json does not exist') || contains(Error.message, 'not found in non-relational') || ...
                        contains(Error.message, 'JSON syntax error'))
                    rethrow(Error)
                else
                    return
                end
            end
            
            % Workaround for MATLAB bug
            fields = fieldnames(obj.graphParameters);
            for i = 1:numel(fields)
                field = fields{i};
                oldData = obj.graphParameters.(field).Data;
                if isequal(oldData,[0;0])
                	obj.graphParameters.(field).Data = [0,0];
                end
                oldExtraConditions = obj.graphParameters.(field).ExtraConditions;
                if isempty(oldExtraConditions)
                    obj.graphParameters.(field).ExtraConditions = 'Add additional (non-sortable) conditions if needed.';
                end
            end
        end
        
        % Public method to save table and graph data to DBs.
        function saveData(obj)
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};           
            
            obj.save_rldb([path 'TransistorDB\Transistor_StandardTable.db']);
            obj.save_nrdb([path 'TransistorDB\Transistor_Graph.json'], ...
                          [path 'TransistorDB\Transistor_AdditionalTable.json']);
        end
   
    end
    
    
    % Private Methods for GUIs
    methods (Access = protected)
        %%% Datasheet Table Parameters GUI
        % Recieves GUI object handles and checks data before storing it into tableParameters
        function saveTableBtnPushed(obj, partNumberBox, typeBox, materialBox, minBoxes, typBoxes, maxBoxes)
            % Track issues
            issues = {};
            
            % Part Number
            new_partNumber = upper(strip(partNumberBox.Value{1}));
            if numel(partNumberBox.Value) > 1
                issues = [issues; 'Part number should be one line of text.'];
            elseif isempty(strip(partNumberBox.Value{1}))
                issues = [issues; 'Please add the manufacturer part number.'];
            elseif isempty(obj.partNumber) || ~strcmp(obj.partNumber, strip(partNumberBox.Value{1}))
                devices = obj.listDevices();
                obj.partNumber = Transistor.checkDuplicateDevices(strip(partNumberBox.Value{1}), devices); 
                partNumberBox.Value = obj.partNumber;
            end
            
            % Transistor Type and Material
            new_type = typeBox.Value;
            if strcmp(new_type(1),'-')
                issues = [issues; 'Please select a transistor type.'];
            end
            new_material = materialBox.Value;
            if strcmp(new_material(1),'-')
                issues = [issues; 'Please select a transistor material.'];
            end            
            
            new_minVals = cell(1,numel(obj.tableFields));
            new_typVals = cell(1,numel(obj.tableFields));
            new_maxVals = cell(1,numel(obj.tableFields));
            % Table Parameters
            for i = 1:numel(obj.tableFields)
                if numel(minBoxes(i).Value) > 1 || numel(typBoxes(i).Value) > 1 || numel(maxBoxes(i).Value) > 1
                    issues = [issues; ['Parameter values should be one line of text: ' obj.tableFields{i}]];
                else
                    minStr = minBoxes(i).Value{1};
                    if ~isempty(minStr)
                        if ~all(ismember(strip(minStr), '0123456789+-.'))
                            issues = [issues; ['Invalid min number for: ' obj.tableFields{i} '. Do not use exponentials']];
                        else
                            try
                                new_minVal = str2double(strip(minStr));
                                if isnan(new_minVal)
                                    error('NaN');
                                else
                                    new_minVals{i} = new_minVal;
                                end
                            catch
                                issues = [issues; ['Invalid min number for: ' obj.tableFields{i} '. Do not use exponentials']];
                            end
                        end
                    end
                    
                    typStr = typBoxes(i).Value{1};
                    if ~isempty(typStr)
                        if ~all(ismember(strip(typStr), '0123456789+-.'))
                            issues = [issues; ['Invalid typ number for: ' obj.tableFields{i} '. Do not use exponentials']];
                        else
                            try
                                new_typVal = str2double(strip(typStr));
                                if isnan(new_typVal)
                                    error('NaN');
                                else
                                    new_typVals{i} = new_typVal;
                                end
                            catch
                                issues = [issues; ['Invalid typ number for: ' obj.tableFields{i} '. Do not use exponentials']];
                            end
                        end
                    end  
                    
                    maxStr = maxBoxes(i).Value{1};
                    if ~isempty(maxStr)
                        if ~all(ismember(strip(maxStr), '0123456789+-.'))
                            issues = [issues; ['Invalid max number for: ' obj.tableFields{i} '. Do not use exponentials']];
                        else
                            try
                                new_maxVal = str2double(strip(maxStr));
                                if isnan(new_maxVal)
                                    error('NaN');
                                else
                                    new_maxVals{i} = new_maxVal;
                                end
                            catch
                                issues = [issues; ['Invalid max number for: ' obj.tableFields{i} '. Do not use exponentials']];
                            end
                        end
                    end                    
                end                
            end
            
            % Show issues or add data and close figure if there are no issues
            if numel(issues) > 0
                msgbox(issues);
            else
                obj.partNumber = new_partNumber;
                obj.type = new_type;
                obj.material = new_material;
                for i=1:numel(obj.tableFields)
                    if ~isempty(new_minVals{i})
                        obj.tableParameters.(obj.tableFields{i}).Min = new_minVals{i};
                    end
                    if ~isempty(new_typVals{i})
                        obj.tableParameters.(obj.tableFields{i}).Typ = new_typVals{i};
                    end
                    if ~isempty(new_maxVals{i})
                        obj.tableParameters.(obj.tableFields{i}).Max = new_maxVals{i};
                    end
                end
                msgbox('Table Parameters Saved');
            end
                
        end 
     
        
        %%% Additional Table Parameters GUI
        % Recieves GUI object handles and checks with user before adding new additional table parameter
        function saveAdditionalBtnPushed(obj,mainText,paramDropdown,typeDropdown,valueBox,conditionDropdowns,conditionBoxes,extraBox,ExpBox)
            % Return if misssing data
            if strcmp(paramDropdown.Value,'--') || strcmp(typeDropdown.Value,'--') || ...
                    (numel(valueBox.Value) == 1 && isempty(valueBox.Value{1}))
                return;
            end
            
            % Get data from boxes
            splitParam = strsplit(paramDropdown.Value,' ');
            new_param = splitParam{1};
            new_type = typeDropdown.Value;
            try 
                new_value = str2double(valueBox.Value{1});
                if isnan(new_value)
                    error('Invalid Parameter Value');
                end
            catch
                msgbox(['Invalid Parameter Value: ' valueBox.Value{1}],'Error');
                return;
            end
            new_struct.Value = new_value;
            new_struct.Conditions = struct();
            for i = 1:6
                try 
                    conditionValue = str2double(conditionBoxes(i).Value{1});
                    if isnan(conditionValue)
                        error('Invalid Parameter Value');
                    end
                catch
                    conditionValue = [];
                end
                splitCondition = strsplit(conditionDropdowns(i).Value,' ');
                condition = splitCondition{1};                
                if ~strcmp(condition,'--') && ~isempty(conditionValue)
                    new_struct.Conditions.(condition) = conditionValue;
                end
            end
            new_struct.ExtraConditions = extraBox.Value{1};
            new_struct.Exp = ExpBox.Value;
            
            % Check with user first
            newStr = [new_param ' ' typeDropdown.Value ' = ' num2str(new_value)];
            index = find(strcmp(obj.conditionFields, new_param));
            newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
            conditions = fieldnames(new_struct.Conditions);
            if numel(conditions) > 0
                newStr = [newStr ' when'];
            end
            for k = 1:numel(conditions)
                condition = conditions{k};
                if k ~= 1
                    newStr = [newStr ','];
                end
                newStr = [newStr ' ' condition ' = ' num2str(new_struct.Conditions.(condition))];
                index = find(strcmp(obj.conditionFields, condition));
                newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
            end
            if strcmp(new_struct.ExtraConditions,'Add additional (non-sortable) conditions if needed.')
                new_struct.ExtraConditions = '';
            end
            if ~isempty(new_struct.ExtraConditions)
                newStr = [newStr ', ' new_struct.ExtraConditions];
            end
            if new_struct.Exp
                newStr = [newStr ', (experimental)'];
            end
            answer = questdlg(['Add new additional table parameter? : ' newline newStr], 'Add New Parameter?', 'Yes', 'No', 'Yes');
            if ~strcmp(answer,'Yes')
                return
            end
                        
            % Add data to additional table parameters structure
            if isfield(obj.additionalTableParameters,new_param) && isfield(obj.additionalTableParameters.(new_param),new_type)
                flag = 0;
                for i = 1:numel(obj.additionalTableParameters.(new_param).(new_type))
                    if isequaln(obj.additionalTableParameters.(new_param).(new_type)(i),new_struct)
                        flag = 1;
                        break;
                    end
                end
                if flag
                    msgbox('Identical parameter already exists.','Duplicate Parameter');
                else
                    new_index = numel(obj.additionalTableParameters.(new_param).(new_type)) + 1;
                    obj.additionalTableParameters.(new_param).(new_type)(new_index) = new_struct;
                end
            else
                obj.additionalTableParameters.(new_param).(new_type) = new_struct;
            end
            
            % Update main text box on GUI
            obj.loadAdditionalDataGUI(mainText);
        end
        
        % Loads additional table parameters onto GUI text area for viewing in text form
        function loadAdditionalDataGUI(obj,mainText)
            mainText.Value = '';
            parameters = fieldnames(obj.additionalTableParameters);
            for i = 1:numel(parameters)
                parameter = parameters{i};
                if isfield(obj.additionalTableParameters.(parameter),'Min')
                    for j = 1:numel(obj.additionalTableParameters.(parameter).Min)
                        s = obj.additionalTableParameters.(parameter).Min(j);
                        newStr = [parameter ' minimum = ' num2str(s.Value)];
                        index = find(strcmp(obj.conditionFields, parameter));
                        newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
                        conditions = fieldnames(s.Conditions);
                        if numel(conditions) > 0
                            newStr = [newStr ' when'];
                        end
                        for k = 1:numel(conditions)
                            condition = conditions{k};
                            if k ~= 1
                                newStr = [newStr ','];
                            end
                            newStr = [newStr ' ' condition ' = ' num2str(s.Conditions.(condition))];
                            index = find(strcmp(obj.conditionFields, condition));
                            newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
                        end
                        if ~isempty(s.ExtraConditions)
                            newStr = [newStr ', ' s.ExtraConditions];
                        end
                        if s.Exp
                            newStr = [newStr ', (experimental)'];
                        end
                        mainText.Value = [mainText.Value ; newStr];
                    end
                end
                if isfield(obj.additionalTableParameters.(parameter),'Typ')
                    for j = 1:numel(obj.additionalTableParameters.(parameter).Typ)
                        s = obj.additionalTableParameters.(parameter).Typ(j);
                        newStr = [parameter ' typical = ' num2str(s.Value)];
                        index = find(strcmp(obj.conditionFields, parameter));
                        newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
                        conditions = fieldnames(s.Conditions);
                        if numel(conditions) > 0
                            newStr = [newStr ' when'];
                        end
                        for k = 1:numel(conditions)
                            condition = conditions{k};
                            if k ~= 1
                                newStr = [newStr ','];
                            end
                            newStr = [newStr ' ' condition ' = ' num2str(s.Conditions.(condition))];
                            index = find(strcmp(obj.conditionFields, condition));
                            newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
                        end
                        if ~isempty(s.ExtraConditions)
                            newStr = [newStr ', ' s.ExtraConditions];
                        end
                        if s.Exp
                            newStr = [newStr ', (experimental)'];
                        end                        
                        mainText.Value = [mainText.Value ; newStr];
                    end
                end
                if isfield(obj.additionalTableParameters.(parameter),'Max')
                    for j = 1:numel(obj.additionalTableParameters.(parameter).Max)
                        s = obj.additionalTableParameters.(parameter).Max(j);
                        newStr = [parameter ' maximum = ' num2str(s.Value)];
                        index = find(strcmp(obj.conditionFields, parameter));
                        newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
                        conditions = fieldnames(s.Conditions);
                        if numel(conditions) > 0
                            newStr = [newStr ' when'];
                        end
                        for k = 1:numel(conditions)
                            condition = conditions{k};
                            if k ~= 1
                                newStr = [newStr ','];
                            end
                            newStr = [newStr ' ' condition ' = ' num2str(s.Conditions.(condition))];
                            index = find(strcmp(obj.conditionFields, condition));
                            newStr = [newStr strrep(obj.conditionDefaultMultipliers{index},'1','') obj.conditionUnits{index}];
                        end
                        if ~isempty(s.ExtraConditions)
                            newStr = [newStr ', ' s.ExtraConditions];
                        end
                        if s.Exp
                            newStr = [newStr ', (experimental)'];
                        end
                        mainText.Value = [mainText.Value ; newStr];
                    end
                end
            end
            if isempty(mainText.Value{1}) && numel(mainText.Value) > 1
                mainText.Value = mainText.Value(2:end);
            end
        end
        
        
        %%% Graph Parameters GUI
        % Dropdown selection function in graph parameters GUI
        function selection(obj, selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras)
            % Check page
            obj.pageChecks(selectDropdown, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);
            
            % Change to new page selection
            obj.current_graphPage = 1;
            obj.current_graphParameter = strrep(selectDropdown.Value,' vs. ','VERSUS');
            % Add new parameter if new page selected
            if strcmp(obj.current_graphParameter, 'NEWPARAMETER')
                obj.temp_graphParameters.NEWPARAMETER.Data = [0,0];
                obj.temp_graphParameters.NEWPARAMETER.Multiplier = ['1','1'];
                obj.temp_graphParameters.NEWPARAMETER.Log = [0,0];
                obj.temp_graphParameters.NEWPARAMETER.Conditions = struct();
                obj.temp_graphParameters.NEWPARAMETER.ExtraConditions = '';
                obj.temp_graphParameters.NEWPARAMETER.Exp = 0;
            end
            loadGraphDataGUI(obj, dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, x_multiplierBox, ...
                y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras)
            obj.updateBtnPushed(dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, y_multiplierBox, x_multiplierBox);        
        end
        
        % Previous button pushed in graph parameters GUI
        function prevBtnPushed(obj, selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras)
            % Check page
            obj.pageChecks(selectDropdown, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);

            if obj.current_graphPage == 1
                return
            end
            
            obj.current_graphPage = obj.current_graphPage - 1;
            obj.loadGraphDataGUI(dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, x_multiplierBox, ...
                y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);
            obj.updateBtnPushed(dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, y_multiplierBox, x_multiplierBox);             
        end
        
        % Save data for graph parameters GUI from temp_graphParameters to graphParameters
        function saveGraphBtnPushed(obj, selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras)
            % Check page
            obj.pageChecks(selectDropdown, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);
                                    
            % Update Plot
            obj.loadGraphDataGUI(dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, x_multiplierBox, ...
                y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);
            obj.updateBtnPushed(dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, y_multiplierBox, x_multiplierBox);  
            
            % Save data
            fields = fieldnames(obj.temp_graphParameters);
            for i = 1:numel(fields)
                field = fields{i};
                if ~strcmp(field,'NEWPARAMETER')
                    for j = 1:numel(obj.temp_graphParameters.(field))
                        param = obj.temp_graphParameters.(field)(j);
                        if ~isfield(param, 'Valid')
                            obj.graphParameters.(field)(j) = param;
                        elseif param.Valid
                            obj.graphParameters.(field)(j) = rmfield(param,'Valid');
                        end
                    end
                end
            end
            msgbox('Data saved.')
        end
        
        % Next button pushed in graph parameters GUI
        function nextBtnPushed(obj, selectDropdown, dataPlot, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras)
            % Check page
            checkResult = obj.pageChecks(selectDropdown, dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, ...
                x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);
            % Don't add another page until the current one is valid
            if ~checkResult
                return
            end 
    
            % Create a new page
            if obj.current_graphPage+1 > numel(obj.temp_graphParameters.(obj.current_graphParameter))
                if strcmp(obj.current_graphParameter,'NEWPARAMETER')
                    msgbox('Enter information for current page before adding another page.');
                    return
                end
                
                obj.current_graphPage = obj.current_graphPage + 1;
                % Clear page for new graph parameter
                dataBox.Value = '0,0';
                expCheckbox.Value = 0;
                for i = 1:5
                   conditionBoxes{i,2}.Value = '';
                end
                conditionExtras.Value = 'Add additional (non-sortable) conditions if needed.';
                
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Data = [0,0];
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Multiplier = ['1','1'];
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Log = [0,0];
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Conditions = struct();
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).ExtraConditions = '';
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Exp = 0;
                
                obj.updateBtnPushed(dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, ...
                    x_logBox, y_multiplierBox, x_multiplierBox);
            
            % Load data for next page
            else
                obj.current_graphPage = obj.current_graphPage + 1;
                obj.loadGraphDataGUI(dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, x_multiplierBox, ...
                    y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras);
                obj.updateBtnPushed(dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, ...
                    x_logBox, y_multiplierBox, x_multiplierBox);
            end
        end
        
        % Load data from temp graph parameters onto GUI
        function loadGraphDataGUI(obj, dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, x_multiplierBox, ...
                y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras)
            data = obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Data;
            if length(data(:,1)) == 1 && data(1,1) == 0
                dataBox.Value = '';
            else
                dataBox.Value = Transistor.dataMatToStr(data);
            end
            
            if strcmp(obj.current_graphParameter,'NEWPARAMETER')
                x_fieldBox.Value = '--Select X-axis Field--';
                y_fieldBox.Value = '--Select Y-axis Field--';
            else
                split_str = strsplit(obj.current_graphParameter,'VERSUS');
                xField = split_str{2};
                yField = split_str{1}; 
                
                xUnit_index = strcmp(obj.graphFields,xField);
                xUnit = obj.graphUnits{xUnit_index};
                yUnit_index = strcmp(obj.graphFields,yField);
                yUnit = obj.graphUnits{yUnit_index};                    

                x_fieldBox.Value = {[xField ' (' xUnit ')']};
                y_fieldBox.Value = {[yField ' (' yUnit ')']};
            end

            x_multiplierBox.Value = {obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Multiplier(1)};
            y_multiplierBox.Value = {obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Multiplier(2)};

            log = obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Log;
            if log(1)  
                x_logBox.Value = 'Log';
            else
                x_logBox.Value = 'Linear';
            end
            if log(2)
                y_logBox.Value = 'Log';
            else
                y_logBox.Value = 'Linear';
            end
            
            if obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Exp
                expCheckbox.Value = 1;
            else
                expCheckbox.Value = 0;
            end
            
            Conditions = fieldnames(obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Conditions);
            for i = 1:numel(Conditions)
                Condition = Conditions{i};
                new_value = obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Conditions.(Condition);
                index = find(strcmp(obj.conditionFields,Condition));
                mult = strrep(obj.conditionDefaultMultipliers{index},'1','');
                unit = obj.conditionUnits{index};
                conditionField = [Condition ' (' mult unit ')'];
                conditionBoxes{i,1}.Value = conditionField;
                conditionBoxes{i,2}.Value = num2str(new_value);
            end
            if (numel(Conditions) < 5)
                for i = numel(Conditions)+1:5
                    conditionBoxes{i,1}.Value = '--Add Condition--';
                    conditionBoxes{i,2}.Value = '';
                end
            end
            
            conditionExtras.Value = obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).ExtraConditions;
        end
        
        % Store data for current page and check if the NEWPARAMETER page can be replaced with a real page
        function checkResult = pageChecks(obj, selectDropdown, dataBox, y_fieldBox, x_fieldBox, y_logBox, ...
                x_logBox, x_multiplierBox, y_multiplierBox, expCheckbox, conditionBoxes, conditionExtras)
                    
            issues = {};
            % Store temp data for current page
            obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Data = Transistor.dataStrToMat(dataBox);
            obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Multiplier = [x_multiplierBox.Value,y_multiplierBox.Value];
            xLog = strcmp(x_logBox.Value,'Log');
            yLog = strcmp(y_logBox.Value,'Log');
            obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Log = [xLog,yLog];
            obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Exp = expCheckbox.Value;
            if isfield(obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage),'Conditions')
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Conditions = struct;
            end
            for i = 1:5
                conditionLabel = conditionBoxes{i,1};
                if ~strcmp(conditionLabel.Value,'--Add Condition--')
                    split_str = strsplit(conditionLabel.Value,' (');
                    field = split_str{1};
                    data = conditionBoxes{i,2}.Value;
                    if numel(data) > 1
                        issues = [issues; ['Condition values should be one line of text: ' field]];
                    else
                        dataStr = data{1};
                        if ~isempty(dataStr)
                            if ~all(ismember(strip(dataStr), '0123456789+-.e'))
                                issues = [issues; ['Invalid condition number for: ' field '. Exponentials allowed (ex. 1.4e-9)']];
                            else
                                try
                                    obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Conditions.(field) = str2double(strip(dataStr));
                                catch
                                    issues = [issues; ['Invalid condition number for: ' field '. Exponentials allowed (ex. 1.4e-9)']];
                                end
                            end
                        end
                    end
                end
            end
            obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).ExtraConditions = conditionExtras.Value;
            
            % Check for missing fields
            if strcmp(y_fieldBox.Value,'--Select Y-axis Field--') || strcmp(x_fieldBox.Value,'--Select X-axis Field--')
                msgbox('Missing Y or X fields', 'Missing fields')
                checkResult = 0;
                return
            end
            
            % Check if NEWPARAMETER page updated
            if strcmp(obj.current_graphParameter,'NEWPARAMETER') && ...
                ~strcmp(y_fieldBox.Value,'--Select Y-axis Field--') && ...
                ~strcmp(x_fieldBox.Value,'--Select X-axis Field--')
                split_yField = strsplit(y_fieldBox.Value,' ');
                split_xField = strsplit(x_fieldBox.Value,' ');
                newFieldName = [split_yField{1} 'VERSUS' split_xField{1}];
                if any(strcmp(fieldnames(obj.temp_graphParameters),newFieldName))
                    checkResult = 0;
                    obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Valid = 0;
                    msgbox(['Parameter ' split_yField{1} ' vs. ' split_xField{1} ' already exists. Select its dropdown to add pages.']);
                    return
                end
                obj.temp_graphParameters.(newFieldName) = obj.temp_graphParameters.NEWPARAMETER;
                obj.temp_graphParameters = rmfield(obj.temp_graphParameters,'NEWPARAMETER');
                obj.current_graphParameter = newFieldName;
                % If NEWPARAMETER updated to a parameter, add a new NEWPARAMETER
                selectDropdown.Items = strrep(fieldnames(obj.temp_graphParameters),'VERSUS',' vs. ');
                selectDropdown.Items = [selectDropdown.Items, 'NEWPARAMETER'];
                selectDropdown.Value = strrep(obj.current_graphParameter,'VERSUS',' vs. ');
                obj.loadGraphDataGUI(dataBox, y_fieldBox, x_fieldBox, y_logBox, x_logBox, x_multiplierBox, y_multiplierBox);
            end
            
            % Check for changed fields
            parts = strsplit(obj.current_graphParameter,'VERSUS');
            xField = strip(strsplit(x_fieldBox.Value(1:end-1), '('));
            yField = strip(strsplit(y_fieldBox.Value(1:end-1), '('));
            if ~strcmp(obj.current_graphParameter,'NEWPARAMETER') && (~strcmp(parts{1},yField{1}) || ~strcmp(parts{2},xField{1}))
                split_str = strsplit(obj.current_graphParameter,'VERSUS');
                xField = split_str{2};
                yField = split_str{1}; 
                
                xUnit_index = strcmp(obj.graphFields,xField);
                xUnit = obj.graphUnits{xUnit_index};
                yUnit_index = strcmp(obj.graphFields,yField);
                yUnit = obj.graphUnits{yUnit_index};                    

                x_fieldBox.Value = {[xField ' (' xUnit ')']};
                y_fieldBox.Value = {[yField ' (' yUnit ')']};
                
                checkResult = 0;
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Valid = 0;
                msgbox('Parameter fields cannot be changed once set. If this is a new parameter, leave this blank and add a new corrected parameter.');
                return;
            end
            
            % Check for data points
            dataMat = Transistor.dataStrToMat(dataBox);
            x = dataMat(:,1);
            y = dataMat(:,1);
            if numel(x) <= 9 || numel(y) <= 9
                issues = [issues; 'Need at least 10 data points'];
            end
            
            % Check ascending X values
            if ~issorted(x)
                issues = [issues; 'X values should be in ascending order'];
            end
            
            % Check for points every 10% of domain
            xMin = min(x);
            xMax = max(x);
            linRange = xMax-xMin;
            linThresholds = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1]*linRange+xMin;
            linSpacing = true;
            linMissing = [];
            for i = 1:numel(linThresholds)-1
               if ~any(x>=linThresholds(i) & x<=linThresholds(i+1))
                   linSpacing = false;
                   linMissing = [linMissing; linThresholds(i), linThresholds(i+1)];
               end
            end
            logSpacing = false;
            if all(x>0)
                logMin = log10(xMin);
                logMax = log10(xMax);
                logRange = logMax-logMin;
                logThresholds = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1]*logRange+logMin;
                logSpacing = true;
                logMissing = [];
                for i = 1:numel(logThresholds)-1
                    if ~any(log10(x)>=logThresholds(i) & log10(x)<=logThresholds(i+1))
                        logSpacing = false;
                        logMissing = [logMissing; 10^logThresholds(i), 10^logThresholds(i+1)];
                    end
                end
            end
            if ~linSpacing && ~logSpacing
                if xLog
                    missingStr = strrep(mat2str(logMissing),';',']  [');
                    issues = [issues; ['Need at least one point in every 10% of domain. For log domain, need points in ranges: ' newline missingStr]];
                else
                    missingStr = strrep(mat2str(linMissing),';',']  [');
                    issues = [issues; ['Need at least one point in every 10% of domain. For linear domain, need points in ranges: ' newline missingStr]];
                end
            end
            
            % Show issues if needed and return
            if ~isempty(issues)
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Valid = 0;
                msgbox([['Invalid parameter page cannot be saved- ' strrep(obj.current_graphParameter,'VERSUS',' vs. ') ...
                    ', page ' num2str(obj.current_graphPage) newline] ; issues]);
                checkResult = 0;
                return
            else
                % No errors, check for unexpected multipliers (eg. Coss (F) instead of Coss (pF))
                yMult = y_multiplierBox.Value;
                xMult = x_multiplierBox.Value;
                expectedYMult = obj.graphDefaultMultipliers{strcmp(obj.graphFields,yField{1})};
                expectedXMult = obj.graphDefaultMultipliers{strcmp(obj.graphFields,xField{1})};
                if ~strcmp(yMult,expectedYMult)
                    answer = questdlg(['Unexpected multiplier: ' yField{1} ' (' strrep(yMult,'1','') yField{2} ').' ...
                        ' Did you intend: ' yField{1} ' (' strrep(expectedYMult,'1','') yField{2} ')?'], ...
                        'Unexpected Y Field Multiplier', ...
                        [yField{1} ' (' strrep(yMult,'1','') yField{2} ')'], ...
                        [yField{1} ' (' strrep(expectedYMult,'1','') yField{2} ')'], ...
                        [yField{1} ' (' strrep(expectedYMult,'1','') yField{2} ')']);
                    if strcmp(answer, [yField{1} ' (' strrep(expectedYMult,'1','') yField{2} ')'])
                        y_multiplierBox.Value = expectedYMult;
                    end
                end
                if ~strcmp(xMult,expectedXMult)
                    answer = questdlg(['Unexpected multiplier: ' xField{1} ' (' strrep(xMult,'1','') xField{2} ').' ...
                        ' Did you intend: ' xField{1} ' (' strrep(expectedXMult,'1','') xField{2} ')?'], ...
                        'Unexpected X Field Multiplier', ...
                        [xField{1} ' (' strrep(xMult,'1','') xField{2} ')'], ...
                        [xField{1} ' (' strrep(expectedXMult,'1','') xField{2} ')'], ...
                        [xField{1} ' (' strrep(expectedXMult,'1','') xField{2} ')']);
                    if strcmp(answer, [xField{1} ' (' strrep(expectedXMult,'1','') xField{2} ')'])
                        x_multiplierBox.Value = expectedXMult;
                    end
                end                    
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Multiplier = [x_multiplierBox.Value,y_multiplierBox.Value];
                % Return
                obj.temp_graphParameters.(obj.current_graphParameter)(obj.current_graphPage).Valid = 1;
                checkResult = 1;
                return
            end
        end
    end
   
    
    % Other private methods
    methods (Access = protected)
        % Protected method to save relational DB information to .db file
        function save_rldb(obj, rldb_file)
            
            % Create DB file if it does not exist
            if ~isfile(rldb_file)
                conn = sqlite(rldb_file, 'create');
                rldbInit = ['create table Transistor ' ...
                            '(partNumber varchar(255) NOT NULL UNIQUE, ' ...
                            'type varchar(255) NOT NULL, ' ...
                            'material varchar(255) NOT NULL'];
                for parameter = obj.tableFields'
                    rldbInit = [rldbInit ', ' ...
                                parameter{1} '_Min NUMERIC, ' ...
                                parameter{1} '_Typ NUMERIC, ' ...
                                parameter{1} '_Max NUMERIC'];
                end
                rldbInit = [rldbInit ')'];
                exec(conn, rldbInit);
                close(conn);
                
            else
                conn = sqlite(rldb_file);
                
                % Check if table exists
                rldb_tblCheck = 'SELECT name FROM sqlite_master WHERE type="table" AND name="Transistor"';
                data = fetch(conn,rldb_tblCheck);
                if isempty(data)
                    rldbInit = ['create table Transistor ' ...
                            '(partNumber varchar(255) NOT NULL UNIQUE, ' ...
                            'type varchar(255) NOT NULL, ' ...
                            'material varchar(255) NOT NULL'];
                    for parameter = obj.tableFields'
                        rldbInit = [rldbInit ', ' ...
                                    parameter{1} '_Min NUMERIC, ' ...
                                    parameter{1} '_Typ NUMERIC, ' ...
                                    parameter{1} '_Max NUMERIC'];
                    end
                    rldbInit = [rldbInit ')'];
                    exec(conn, rldbInit);
                end
                
                data = fetch(conn,'SELECT sql FROM sqlite_master WHERE tbl_name="Transistor" AND type="table"');
                if isempty(data)
                    columnNames = {};
                else
                    columnCells = strip(split(data,','));
                    columnNames = {};
                    for i = 4:numel(columnCells)
                       columnCell = columnCells(i);
                       splitCell = split(columnCell{1}, ' ');
                       columnNames(i-1) = splitCell(1);
                    end
                end
                for parameter = obj.tableFields'
                    if ~any(strcmp(columnNames,[parameter{1},'_Typ']))
                        alter_rldb = ['ALTER TABLE Transistor ADD ' parameter{1} '_Min NUMERIC'];
                        exec(conn,alter_rldb);
                        alter_rldb = ['ALTER TABLE Transistor ADD ' parameter{1} '_Typ NUMERIC'];
                        exec(conn,alter_rldb);
                        alter_rldb = ['ALTER TABLE Transistor ADD ' parameter{1} '_Max NUMERIC'];
                        exec(conn,alter_rldb);
                    end
                end
                
                close(conn)
            end
            
            % Connect and check if row exists for this device
            conn = sqlite(rldb_file);
            rldb_check = ['SELECT 1 FROM Transistor WHERE partNumber="' obj.partNumber '"'];
            result = fetch(conn,rldb_check);
            % If it doesn't exist, add new device row
            if size(result,1) == 0
                rldb_insertNewEntry = ['INSERT INTO Transistor (partNumber, type, material) VALUES ("' obj.partNumber '", "' obj.type '", "' obj.material '")'];
                exec(conn, rldb_insertNewEntry);
            % If it does exist, update type/material if needed
            else
                result = fetch(conn,['SELECT type,material FROM Transistor WHERE partNumber="' obj.partNumber '"']);
                if ~strcmp(result{1},obj.type)
                    exec(conn,['UPDATE Transistor SET type="' obj.type '" WHERE partNumber="' obj.partNumber '"']);
                end
                if ~strcmp(result{2},obj.material)
                    exec(conn,['UPDATE Transistor SET material="' obj.material '" WHERE partNumber="' obj.partNumber '"']);
                end                
            end
            
            % Execute SQL UPDATE command for each cell for which we have data
            for parameter = obj.tableFields'
                if isfield(obj.tableParameters, parameter{1})
                    field = obj.tableParameters.(parameter{1});
                    if isfield(field, 'Min')
                        rldb_insert = ['UPDATE Transistor SET ' parameter{1} '_Min=' num2str(field.Min) ' WHERE partNumber="' obj.partNumber '"'];
                        exec(conn, rldb_insert);
                    end
                    if isfield(field, 'Typ')
                        rldb_insert = ['UPDATE Transistor SET ' parameter{1} '_Typ=' num2str(field.Typ) ' WHERE partNumber="' obj.partNumber '"'];
                        exec(conn, rldb_insert);
                    end
                    if isfield(field, 'Max')
                        rldb_insert = ['UPDATE Transistor SET ' parameter{1} '_Max=' num2str(field.Max) ' WHERE partNumber="' obj.partNumber '"'];
                        exec(conn, rldb_insert);
                    end
                end
            end
            close(conn)
        end
        
        % Protected method to save non-relational DB information to .json file
        function save_nrdb(obj, graph_nrdb_file, additionalTable_nrdb_file)
            % Create graph parameters JSON DB if it does not exist and add/replace information
            if ~isfile(graph_nrdb_file)
                fid = fopen(graph_nrdb_file, 'w');
                Graph_Structure = struct();
                Graph_Structure.(obj.partNumber) = obj.graphParameters;
                json_text = jsonencode(Graph_Structure);
                fprintf(fid, json_text);
                fclose(fid);
            % Add/replace information on graph parameters
            else
                old_json_text = fileread(graph_nrdb_file);
                Graph_Structure = jsondecode(old_json_text);
                Graph_Structure.(obj.partNumber) = obj.graphParameters;
                json_text = jsonencode(Graph_Structure);
                fid = fopen(graph_nrdb_file, 'w');
                fprintf(fid, json_text);
                fclose(fid);
            end
            
            % Create additional table parameters JSON DB if it does not exist and add/replace information
            if ~isfile(additionalTable_nrdb_file)
                fid = fopen(additionalTable_nrdb_file, 'w');
                additionalTable_Structure = struct();
                additionalTable_Structure.(obj.partNumber) = obj.additionalTableParameters;
                json_text = jsonencode(additionalTable_Structure);
                fprintf(fid, json_text);
                fclose(fid);
            % Add/replace information on additional table parameters
            else
                old_json_text = fileread(additionalTable_nrdb_file);
                additionalTable_Structure = jsondecode(old_json_text);
                additionalTable_Structure.(obj.partNumber) = obj.additionalTableParameters;
                json_text = jsonencode(additionalTable_Structure);
                fid = fopen(additionalTable_nrdb_file, 'w');
                fprintf(fid, json_text);
                fclose(fid);
            end            
        end
        
        % Protected method to load relational DB information from .db file
        function load_rldb(obj, rldb_file)
            
            % Connect to SQL DB
            if ~isfile(rldb_file)
                error(['Relational database file ' rldb_file ' does not exist.']);
            end
            conn = sqlite(rldb_file);
            
            % Fetch transistor type and material. This will error and be caught if part does not exist in DB
            response = fetch(conn, ['SELECT type,material from Transistor WHERE partNumber="' obj.partNumber '"']);
            obj.type = response{1};
            obj.material = response{2};
            
            % Fetch Min, Typ, and Max for each table parameter.
            % sqlite.fetch() throws error on NULL cells, so try each individually and catch errors.
            for parameter = obj.tableFields'
                try
                    rldb_fetchMin = ['SELECT ' parameter{1} '_Min from Transistor WHERE partNumber="' obj.partNumber '"'];
                    min = fetch(conn, rldb_fetchMin);
                    obj.tableParameters.(parameter{1}).Min = min{1};
                catch Error % Ignore null cells
                    if ~(contains(Error.message,'Unexpected NULL') || contains(Error.message,'no such column'))
                        rethrow(Error);
                    end
                end
                try
                    rldb_fetchTyp = ['SELECT ' parameter{1} '_Typ from Transistor WHERE partNumber="' obj.partNumber '"'];
                    typ = fetch(conn, rldb_fetchTyp);
                    obj.tableParameters.(parameter{1}).Typ = typ{1};
                catch Error % Ignore null cells
                    if ~(strcmp(Error.message(1:15),'Unexpected NULL') || contains(Error.message,'no such column'))
                        rethrow(Error);
                    end
                end
                try
                    rldb_fetchMax = ['SELECT ' parameter{1} '_Max from Transistor WHERE partNumber="' obj.partNumber '"'];
                    max = fetch(conn, rldb_fetchMax);
                    obj.tableParameters.(parameter{1}).Max = max{1};
                catch Error % Ignore null cells
                    if ~(strcmp(Error.message(1:15),'Unexpected NULL') || contains(Error.message,'no such column'))
                        rethrow(Error);
                    end
                end
            end
            close(conn);
        end
        
        % Protected method to load non-relational DB information from .json file
        function data_struct = load_nrdb(obj, nrdb_file)
            % Connect to nrdb parameters JSON DB
            if ~isfile(nrdb_file)
                error(['Non-relational database file ' nrdb_file ' does not exist.']);
            end
            % Open JSON object and load in object with key partNumber
            json_text = fileread(nrdb_file);
            top_struct = jsondecode(json_text);
            if ~isfield(top_struct, obj.partNumber)
                error(['Transistor part ' obj.partNumber ' not found in non-relational database ' nrdb_file]);
            end
            data_struct = top_struct.(obj.partNumber);
        end
    end
    
    
    % Private, Static Methods
    methods (Static, Access = protected) 
        % Read in meta-data from associated text files
        function [types,materials,...
                tableFields,tableUnits,tableDefaultMultipliers,...
                graphFields,graphUnits,graphDefaultMultipliers,...
                conditionFields,conditionUnits,conditionDefaultMultipliers] = loadMetaData()
            types = strip(split(fileread('Transistor_types.txt'),','));
            materials = strip(split(fileread('Transistor_materials.txt'),','));
            
            lines = splitlines(fileread('Transistor_tableFields.txt'));
            tableFields = strip(split(lines(1), ','));
            tableUnits = strip(split(lines(2), ','));
            tableDefaultMultipliers = strip(split(lines(3), ','));
            
            lines = splitlines(fileread('Transistor_graphFields.txt'));
            graphFields = [tableFields; strip(split(lines(1), ','))];
            graphUnits = [tableUnits; strip(split(lines(2), ','))];
            graphDefaultMultipliers = [tableDefaultMultipliers; strip(split(lines(3), ','))];
            
            lines = splitlines(fileread('Transistor_conditionFields.txt'));
            conditionFields = [graphFields; strip(split(lines(1), ','))];
            conditionUnits = [graphUnits; strip(split(lines(2), ','))];
            conditionDefaultMultipliers = [graphDefaultMultipliers; strip(split(lines(3), ','))]; 
        end
        
        % Merge another SQL table into the default SQL table for this class for table parameters 
        function changedStr = merge_tableParameters(rldb_file, other_rldb_file)
            changedStr = '';
            [~,~,tableFields,tableUnits,tableDefaultMultipliers,~,~,~,~,~,~] = Transistor.loadMetaData();
            
            
            % Connect to original SQL DB
            if ~isfile(rldb_file)
                error(['Relational database file ' rldb_file ' does not exist.']);
            end
            conn1 = sqlite(rldb_file);

            % Connect to other SQL DB
            if ~isfile(other_rldb_file)
                error(['Relational database file ' other_rldb_file ' does not exist.']);
            end
            conn2 = sqlite(other_rldb_file);

            % Get parameters list from original DB
            result = fetch(conn1,'SELECT sql FROM sqlite_master WHERE tbl_name="Transistor" AND type="table"');
            if isempty(result)
                columnNames = {};
            else
                columnCells = strip(split(result,','));
                columnNames = {};
                for i = 2:numel(columnCells)
                    columnCell = columnCells(i);
                    splitCell = split(columnCell{1}, ' ');
                    columnNames(i-1) = splitCell(1);
                end
            end
            parameterNames = {};
            for i = 1:numel(columnNames)
                name = columnNames{i};
                % Get the parameter names from the columns that end in '_Typ'
                if ~isempty(name) && length(name) > 4 && strcmp(name(end-3:end),'_Typ')
                    parameterNames = [parameterNames name(1:end-4)];
                end
            end
            
            % Get data for parameters from original DB
            oldData = struct();
            devices = fetch(conn1,'SELECT partNumber FROM Transistor');
            for i = 1:numel(devices)
                device = devices{i};
                type = fetch(conn1,['SELECT type FROM Transistor WHERE partNumber="' device '"']);
                oldData.(device).type = type{1};
                material = fetch(conn1,['SELECT material FROM Transistor WHERE partNumber="' device '"']);
                oldData.(device).material =  material{1};
                for j = 1:numel(parameterNames)
                    parameter = parameterNames{j}; 
                    
                    min = fetch(conn1, ...
                        ['SELECT ' parameter '_Min FROM Transistor WHERE partNumber="' device '" AND ' parameter '_Min IS NOT NULL']);
                    if ~isempty(min)
                        oldData.(device).(parameter).Min = min{1};
                    end
                    
                    typ = fetch(conn1, ...
                        ['SELECT ' parameter '_Typ FROM Transistor WHERE partNumber="' device '" AND ' parameter '_Typ IS NOT NULL']);
                    if ~isempty(typ)
                        oldData.(device).(parameter).Typ = typ{1};
                    end
                    
                    max = fetch(conn1, ...
                        ['SELECT ' parameter '_Max FROM Transistor WHERE partNumber="' device '" AND ' parameter '_Max IS NOT NULL']);
                    if ~isempty(max)
                        oldData.(device).(parameter).Max = max{1};
                    end                    
                end
            end
                
            % Get data for parameters from other DB
            otherData = struct();
            devices = fetch(conn2,'SELECT partNumber FROM Transistor');
            for i = 1:numel(devices)
                device = devices{i};
                type = fetch(conn2,['SELECT type FROM Transistor WHERE partNumber="' device '"']);
                otherData.(device).type = type{1}; 
                material = fetch(conn2,['SELECT material FROM Transistor WHERE partNumber="' device '"']);
                otherData.(device).material =  material{1};
                for j = 1:numel(parameterNames)
                    parameter = parameterNames{j}; 
                    
                    min = fetch(conn2, ...
                        ['SELECT ' parameter '_Min FROM Transistor WHERE partNumber="' device '" AND ' parameter '_Min IS NOT NULL']);
                    if ~isempty(min)
                        otherData.(device).(parameter).Min = min{1};
                    end
                    
                    typ = fetch(conn2, ...
                        ['SELECT ' parameter '_Typ FROM Transistor WHERE partNumber="' device '" AND ' parameter '_Typ IS NOT NULL']);
                    if ~isempty(typ)
                        otherData.(device).(parameter).Typ = typ{1};
                    end
                    
                    max = fetch(conn2, ...
                        ['SELECT ' parameter '_Max FROM Transistor WHERE partNumber="' device '" AND ' parameter '_Max IS NOT NULL']);
                    if ~isempty(max)
                        otherData.(device).(parameter).Max = max{1};
                    end                    
                end
            end                

            % Merge data
            oldDevices = fieldnames(oldData);
            otherDevices = fieldnames(otherData);
            for i = 1:numel(otherDevices)
                device = otherDevices{i};
                
                % Add any new devices
                if ~any(strcmp(oldDevices,device))
                    insert_cmd = ['INSERT INTO Transistor (partNumber, type, material) VALUES ("' ...
                        device '", "' otherData.(device).type '", "' otherData.(device).material '")'];
                    exec(conn1,insert_cmd);
                    fields = fieldnames(otherData.(device));
                    for j = 1:numel(fields)
                        field = fields{j};
                        if ~any(strcmp({'type','material'},field))
                            if isfield(otherData.(device).(field),'Min')
                                exec(conn1,['UPDATE Transistor SET ' field '_Min=' num2str(otherData.(device).(field).Min) ...
                                    ' WHERE partNumber="' device '"']);
                            end
                            if isfield(otherData.(device).(field),'Typ')
                                exec(conn1,['UPDATE Transistor SET ' field '_Typ=' num2str(otherData.(device).(field).Typ) ...
                                    ' WHERE partNumber="' device '"']);
                            end
                            if isfield(otherData.(device).(field),'Max')
                                exec(conn1,['UPDATE Transistor SET ' field '_Max=' num2str(otherData.(device).(field).Max) ...
                                    ' WHERE partNumber="' device '"']);
                            end  
                        end
                    end
                    
                % Compare existing devices
                else
                    % Different types?
                    if ~strcmp(oldData.(device).type, otherData.(device).type)
                        answer = questdlg(['Conflicting transistor types: ' device '. Which is correct?'], 'Conflicting Type', ...
                            oldData.(device).type, otherData.(device).type, oldData.(device).type);
                        exec(conn1, ['UPDATE Transistor SET type="' answer '" WHERE partNumber="' device '"']);
                        if strcmp(answer, otherData.(device).type)
                            changedStr = [changedStr ',' device ' ' oldData.(device).type '->' otherData.(device).type];
                        end
                    end
                    % Different materials?
                    if ~strcmp(oldData.(device).material, otherData.(device).material)
                        answer = questdlg(['Conflicting transistor materials: ' device '. Which is correct?'], 'Conflicting Material', ...
                            oldData.(device).material, otherData.(device).material, oldData.(device).material);
                        exec(conn1, ['UPDATE Transistor SET material="' answer '" WHERE partNumber="' device '"']);
                        if strcmp(answer, otherData.(device).material)
                            changedStr = [changedStr ',' device ' ' oldData.(device).material '->' otherData.(device).material];
                        end
                    end
                    % Compare parameters
                    for j = 1:numel(parameterNames)
                        parameter = parameterNames{j};
                        % Does other DB have info for this parameter on this device?
                        if isfield(otherData.(device), parameter)
                            % If this DB does not have info for this parameter on this device, add other DB's info
                            if ~isfield(oldData.(device), parameter)
                                if isfield(otherData.(device).(parameter), 'Min')
                                    exec(conn1, ['UPDATE Transistor SET "' parameter '_Min"=' ...
                                        num2str(otherData.(device).(parameter).Min) ' WHERE partNumber="' device '"']);
                                end
                                if isfield(otherData.(device).(parameter), 'Typ')
                                    exec(conn1, ['UPDATE Transistor SET "' parameter '_Typ"=' ...
                                        num2str(otherData.(device).(parameter).Typ) ' WHERE partNumber="' device '"']);
                                end
                                if isfield(otherData.(device).(parameter), 'Max')
                                    exec(conn1, ['UPDATE Transistor SET "' parameter '_Max"=' ...
                                        num2str(otherData.(device).(parameter).Max) ' WHERE partNumber="' device '"']);
                                end
                                
                            % If there is data for parameter in both DBs:
                            else
                                % Check min
                                if isfield(otherData.(device).(parameter), 'Min')
                                    if ~isfield(oldData.(device).(parameter), 'Min')
                                        exec(conn1, ['UPDATE Transistor SET "' parameter '_Min"=' ...
                                            num2str(otherData.(device).(parameter).Min) ' WHERE partNumber="' device '"']); 
                                    else
                                        a = double(oldData.(device).(parameter).Min);
                                        b = double(otherData.(device).(parameter).Min);
                                        if a ~= b
                                            index = strcmp(tableFields,parameter);
                                            unitStr = ['(' strrep(tableDefaultMultipliers{index},'1','') tableUnits{index} ')'];
                                            answer = questdlg(['Conflicting Parameter: ' parameter ' min.  in ' device '. Which is correct?'], ...
                                                'Conflicting Parameter', [num2str(a) unitStr], [num2str(b) unitStr], [num2str(a) unitStr]);
                                            answer = strsplit(answer, '(');
                                            answer = answer{1};
                                            exec(conn1, ['UPDATE Transistor SET "' parameter '_Min"=' ...
                                                num2str(answer) ' WHERE partNumber="' device '"']);
                                            if strcmp(answer, num2str(b))
                                                changedStr = [changedStr ',' device ' ' parameter ' ' num2str(a) '->' num2str(b)];
                                            end                                            
                                        end
                                    end
                                end
                                
                                % Check typ
                                if isfield(otherData.(device).(parameter), 'Typ')
                                    if ~isfield(oldData.(device).(parameter), 'Typ')
                                        exec(conn1, ['UPDATE Transistor SET "' parameter '_Typ"=' ...
                                            num2str(otherData.(device).(parameter).Typ) ' WHERE partNumber="' device '"']);  
                                    else
                                        a = double(oldData.(device).(parameter).Typ);
                                        b = double(otherData.(device).(parameter).Typ);
                                        if a ~= b
                                            index = strcmp(tableFields,parameter);
                                            unitStr = ['(' strrep(tableDefaultMultipliers{index},'1','') tableUnits{index} ')'];
                                            answer = questdlg(['Conflicting Parameter: ' parameter ' typ.  in ' device '. Which is correct?'], ...
                                                'Conflicting Parameter', [num2str(a) unitStr], [num2str(b) unitStr], [num2str(a) unitStr]);
                                            answer = strsplit(answer, '(');
                                            answer = answer{1};
                                            exec(conn1, ['UPDATE Transistor SET "' parameter '_Typ"=' ...
                                                num2str(answer) ' WHERE partNumber="' device '"']);                                          
                                            if strcmp(answer, num2str(b))
                                                changedStr = [changedStr ',' device ' ' parameter ' ' num2str(a) '->' num2str(b)];
                                            end
                                        end
                                    end
                                end
                                
                                % Check max
                                if isfield(otherData.(device).(parameter), 'Max')
                                    if ~isfield(oldData.(device).(parameter), 'Max')
                                        exec(conn1, ['UPDATE Transistor SET "' parameter '_Max"=' ...
                                            num2str(otherData.(device).(parameter).Max) ' WHERE partNumber="' device '"']);  
                                    else
                                        a = double(oldData.(device).(parameter).Max);
                                        b = double(otherData.(device).(parameter).Max);
                                        if a ~= b
                                            index = strcmp(tableFields,parameter);
                                            unitStr = ['(' strrep(tableDefaultMultipliers{index},'1','') tableUnits{index} ')'];
                                            answer = questdlg(['Conflicting Parameter: ' parameter ' max.  in ' device '. Which is correct?'], ...
                                                'Conflicting Parameter', [num2str(a) unitStr], [num2str(b) unitStr], [num2str(a) unitStr]);
                                            answer = strsplit(answer, '(');
                                            answer = answer{1};
                                            exec(conn1, ['UPDATE Transistor SET "' parameter '_Max"=' ...
                                                num2str(answer) ' WHERE partNumber="' device '"']);                                          
                                            if strcmp(answer, num2str(b))
                                                changedStr = [changedStr ',' device ' ' parameter ' ' num2str(a) '->' num2str(b)];
                                            end
                                        end
                                    end
                                end
                            end                                
                        end
                    end
                end
            end

            close(conn1);
            close(conn2);
        end
        
        % Merge another json structure into the default json structure for this class for additional table parameters
        function merge_additionalTable(nrdb_file, other_nrdb_file)
            % Connect to original json file
            if ~isfile(nrdb_file)
                error(['Non-relational database file ' nrdb_file ' does not exist.']);
            end
            oldData = jsondecode(fileread(nrdb_file));
            
            % Connect to other json file
            if ~isfile(other_nrdb_file)
                error(['Non-relational database file ' other_nrdb_file ' does not exist.']);
            end
            otherData = jsondecode(fileread(other_nrdb_file));
            
            % Import data
            newData = oldData;
            devices = fieldnames(otherData);
            for i = 1:numel(devices)
                device = devices{i};
                
                % Add new devices
                if ~any(strcmp(fieldnames(oldData),device))
                    newData.(device) = otherData.(device);
                    
                % Compare parameters
                else
                    parameters = fieldnames(otherData.(device));
                    for j = 1:numel(parameters)
                        parameter = parameters{j};
                        
                        % Add new parameters
                        if ~any(strcmp(fieldnames(oldData.(device)),parameter))
                            newData.(device).(parameter) = otherData.(device).(parameter);

                        % Compare existing parameters with new ones
                        elseif ~isequaln(oldData.(device).(parameter), otherData.(device).(parameter))
                            oldStruct = oldData.(device).(parameter);
                            otherStruct = otherData.(device).(parameter);
                            
                            if isfield(otherStruct,'Min') && ~isfield(oldStruct,'Min')
                                newData.(device).(parameter).Min = otherStruct.Min;
                            elseif isfield(otherStruct,'Min')
                                for k = 1:numel(otherStruct)
                                    duplicateFlag = 0;
                                    for l = 1:numel(oldStruct)
                                        if isequaln(otherStruct.Min(k),oldStruct.Min(l))
                                            duplicateFlag = 1;
                                        end
                                    end
                                    if ~duplicateFlag
                                        newIndex = numel(newData.(device).(parameter).Min) + 1;
                                        newData.(device).(parameter).Min(newIndex) = otherStruct.Min(k);
                                    end
                                end
                            end
                            
                            if isfield(otherStruct,'Typ') && ~isfield(oldStruct,'Typ')
                                newData.(device).(parameter).Typ = otherStruct.Typ;
                            elseif isfield(otherStruct,'Typ')
                                for k = 1:numel(otherStruct)
                                    duplicateFlag = 0;
                                    for l = 1:numel(oldStruct)
                                        if isequaln(otherStruct.Typ(k),oldStruct.Typ(l))
                                            duplicateFlag = 1;
                                        end
                                    end
                                    if ~duplicateFlag
                                        newIndex = numel(newData.(device).(parameter).Typ) + 1;
                                        newData.(device).(parameter).Typ(newIndex) = otherStruct.Typ(k);
                                    end
                                end
                            end
                            
                            if isfield(otherStruct,'Max') && ~isfield(oldStruct,'Max')
                                newData.(device).(parameter).Max = otherStruct.Max;
                            elseif isfield(otherStruct,'Max')
                                for k = 1:numel(otherStruct)
                                    duplicateFlag = 0;
                                    for l = 1:numel(oldStruct)
                                        if isequaln(otherStruct.Max(k),oldStruct.Max(l))
                                            duplicateFlag = 1;
                                        end
                                    end
                                    if ~duplicateFlag
                                        newIndex = numel(newData.(device).(parameter).Max) + 1;
                                        newData.(device).(parameter).Max(newIndex) = otherStruct.Max(k);
                                    end
                                end
                            end
                        end
                    end                    
                end
            end
            
            json_text = jsonencode(newData);
            fid = fopen(nrdb_file, 'w');
            fprintf(fid, json_text);
            fclose(fid);  
        end
        
        % Merge another json structure into the default json structure for this class for graph parameters
        function merge_graphParameters(nrdb_file, other_nrdb_file, direction)
            if ~strcmp(direction,'PUSH') && ~strcmp(direction,'PULL')
                error('Direction argument should be "PUSH" or "PULL"')
            end
            
            % Connect to original json file
            if ~isfile(nrdb_file)
                error(['Non-relational database file ' nrdb_file ' does not exist.']);
            end
            oldData = jsondecode(fileread(nrdb_file));
            
            % Connect to other json file
            if ~isfile(other_nrdb_file)
                error(['Non-relational database file ' other_nrdb_file ' does not exist.']);
            end
            otherData = jsondecode(fileread(other_nrdb_file));  
            
            % Import data
            newData = oldData;
            devices = fieldnames(otherData);
            for i = 1:numel(devices)
                device = devices{i};
                
                % Add new devices
                if ~any(strcmp(fieldnames(oldData),device))
                    newData.(device) = otherData.(device);
                    
                % Compare parameters
                else
                    parameters = fieldnames(otherData.(device));
                    for j = 1:numel(parameters)
                        parameter = parameters{j};
                        
                        % Add new parameters
                        if ~any(strcmp(fieldnames(oldData.(device)),parameter))
                            newData.(device).(parameter) = otherData.(device).(parameter);

                        else
                            % Add new pages of information on this parameter
                            other_param = otherData.(device).(parameter);
                            old_param = oldData.(device).(parameter);
                            for k = 1:numel(other_param)
                                isDuplicate = 0;
                                for l = 1:numel(old_param)
                                    if isequaln(other_param(k),old_param(l))
                                        isDuplicate = 1;
                                    end
                                end
                                if ~isDuplicate
                                    next_index = numel(newData.(device).(parameter)) + 1;
                                    newData.(device).(parameter)(next_index) = other_param(k);
                                end
                            end
                        end
                    end                    
                end
            end            
            
            json_text = jsonencode(newData);
            fid = fopen(nrdb_file, 'w');
            fprintf(fid, json_text);
            fclose(fid);   
        end
        
        % Scale/multiply data by multiplier (unit prefix)
        function scaledData = scaleData(data,multiplier)
            scaledData = data;
            if ~isempty(multiplier)
                if ~contains(Transistor.SIkeys,multiplier)
                    error(['Multiplier: ' string(multiplier) ' is not a valid SI unit prefix']);
                end
                scaledData = data*Transistor.SIprefixes(multiplier);
            end
        end
    end
    
    
    % Public, Static Methods
    methods (Static)
        %%% Methods for pulling and comparing device data
        % Get a table parameter value for a device. 
            % If there is no data, returns NaN. 
            % All other issues such as invalid devices will throw errors.
        function paramValue = getTableParam(device,param,min_typ_max)
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            % Check inputs
            if ~contains(Transistor.listDevices,strip(device))
                error(['Device ' device 'does not exist in DB']);
            end
            if ~contains({'min','typ','max'},lower(strip(min_typ_max)))
                error('Third argument must be Min Typ or Max');
            end
            
            min_typ_max = strip(min_typ_max);
            min_typ_max = [upper(min_typ_max(1)) lower(min_typ_max(2:end))];
            paramValue = NaN;
            
            % Connect to SQL DB
            if ~isfile([path 'TransistorDB\Transistor_StandardTable.db'])
                error(['Relational database file ' path 'TransistorDB\Transistor_StandardTable.db does not exist.']);
            end
            conn = sqlite([path 'TransistorDB\Transistor_StandardTable.db']);
            
            % Get parameter value
            try
                [~,~,tableFields,~,tableDefaultMultipliers,~,~,~,~,~,~] = Transistor.loadMetaData(); 
                index = find(strcmpi(tableFields,strip(param)));
                if isempty(index)
                    error('no such column');
                end
                param = tableFields{index};
                
                result = fetch(conn, ['SELECT ' param '_' min_typ_max ...
                    ' FROM Transistor WHERE partNumber="' strip(device) '"']);
                paramValue = double(result{1});
                
                mult = Transistor.SIprefixes(tableDefaultMultipliers{strcmp(tableFields,param)});
                paramValue = paramValue*mult;
            catch Exception
                if contains(Exception.message,'no such column')
                    error(['Param ' strip(param) ' does not exist in the DB']);
                elseif contains(Exception.message,'Unexpected NULL')
                    % No data, return NaN
                else
                    rethrow(Exception)
                end
            end
            
        end
        
        % Get graph curve for a parameter comparison for a device
            % Returns a cell array of Nx2 arrays.
            % If there is no data, returns a 0x0 cell array.
            % All other issues such as invalid devices will throw errors.
        function [graphData, ConditionsStruct] = getGraphData(device,yParam,xParam)
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            % Check if device exists
            if ~contains(Transistor.listDevices,strip(device))
                error(['Device ' device 'does not exist in DB']);
            end
            if ~isfile([path 'TransistorDB\Transistor_Graph.json'])
                error(['Non-relational database file ' path 'TransistorDB\Transistor_Graph.json does not exist.']);
            end
            json_text = fileread([path 'TransistorDB\Transistor_Graph.json']);
            top_struct = jsondecode(json_text);
            if ~isfield(top_struct, device)
                error(['Device ' device ' does not exist in DB']);
            end
            
            graphData = cell(0);
            ConditionsStruct = struct;
            data_struct = top_struct.(device);
            field = [strip(yParam) 'VERSUS' strip(xParam)];
            if ~isfield(data_struct,field)
                return
            end
            for i = 1:numel(data_struct.(field))
                data = data_struct.(field)(i).Data;
                data(:,1) = Transistor.scaleData(data(:,1), data_struct.(field)(i).Multiplier(1));
                data(:,2) = Transistor.scaleData(data(:,2), data_struct.(field)(i).Multiplier(2));
                graphData(i) = {data};
                ConditionsStruct.Conditions(i) = data_struct.(field)(i).Conditions;
            end
        end
        
        % Returns the y-values and associated conditions for given xParam, yParam, and xValue.
            % Will extrapolate if xValue is out of the domain if it is within a 10% margin of the domain
                % Will return NaN otherwise
            % Returns a ConditionsStruct structure which has one element- a structure array of conditions
            % corresponding to the curves from which the y-values were obtained.
        function [yValue, ConditionsStruct] = getGraphPoint(device,yParam,xParam,xValue)
             [graphData, ConditionsStruct] = Transistor.getGraphData(device,yParam,xParam);
             yValue = [];
             for i = 1:numel(graphData)
                data = graphData{i};
                x = data(:,1);
                y = data(:,2);
                if xValue > min(x) && xValue < max(x)
                    % Spline (3rd order) interpolation works well
                    yValue(i) = interp1(x,y,xValue,'spline');
                elseif xValue < min(x) && xValue > (min(x)-(max(x)-min(x))*0.1) || ...
                       xValue > max(x) && xValue < (max(x)+(max(x)-min(x))*0.1)
                    % Use linear extrapolation, it's the most consistent.
                    yValue(i) = interp1(x,y,xValue,'linear','extrap');
                else
                    % Data out of range, return []
                    yValue(i) = NaN;
                end
             end
        end
        
        % Compare curves for a parameter comparison for one device
            % Ex. Show Id vs Vds curves for different Vgs for one device
        function showGraphData(device,yParam,xParam,axis)
            [graphData, ConditionsStruct] = Transistor.getGraphData(device,yParam,xParam);
            strings = Transistor.ConditionsStructToStr(ConditionsStruct);
            [~,~,~,~,~,graphFields,graphUnits,graphDefaultMultipliers,~,~,~] = Transistor.loadMetaData();
        
            for i = 1:numel(graphData)
                data = graphData{i};
                hold on
                plot(axis,data(:,1),data(:,2),'LineWidth',3);
            end
            yIndex = find(strcmp(graphFields,yParam));
            ylabel(axis,[yParam ' (' strrep(graphDefaultMultipliers{yIndex},'1','') graphUnits{yIndex} ')'])
            xIndex = find(strcmp(graphFields,xParam));
            xlabel(axis,[xParam ' (' strrep(graphDefaultMultipliers{xIndex},'1','') graphUnits{xIndex} ')'])
            if numel(graphData) > 1
                lgnd = legend(axis,strings);
                set(lgnd,'FontSize',12);
            end
            title(axis,[device ': ' yParam ' vs. ' xParam])
        end
        
        % Returns a cell array of strings made from a ConditionsStruct structure from Transistor.getGraphPoint
            % For data visualization purposes
        function strings = ConditionsStructToStr(ConditionsStruct)
            strings = cell(0);
            [~,~,~,~,~,~,~,~,conditionFields,conditionUnits,conditionDefaultMultipliers] = Transistor.loadMetaData();
            for i = 1:numel(ConditionsStruct.Conditions)
                s = ConditionsStruct.Conditions(i);
                str = '';
                conditions = fieldnames(s);
                for j = 1:numel(conditions)
                    condition = conditions{j};
                    if j > 1 
                        str = [str ','];
                    end
                    index = find(strcmp(conditionFields,condition));
                    str = [' ' str condition ' = ' num2str(s.(condition)) '(' ,...
                        strrep(conditionDefaultMultipliers{index},'1',''),...
                        conditionUnits{index} ')'];
                end
                strings(i) = {strip(str)};
            end
        end
  
        % Numerically integrate a parameter curve over a given domain
            % Linearly extrapolates domain to xMin and xMax if out of the domain
        function result = integrateCurve(device, yParam, xParam, xMin, xMax)
            % Get Data
            [graphData, ~] = Transistor.getGraphData(device,yParam,xParam);
            
            result = {};
            for i = 1:numel(graphData)
                X = graphData{i};
                x = X(:,1);
                y = X(:,2);

                % Extrapolate to xMin if needed
                if min(x) > xMin
                    y = [interp1(x,y,0,'linear','extrap'); y];
                    x = [xMin; x];
                end
                % Extrapolate to xMax if needed
                if max(x) < xMax
                    y = [y; interp1(x,y,0,'linear','extrap')];
                    x = [x; xMax];
                end  

                % Find cumulative numerical intergral
                Q = cumtrapz(x,y);
                sumMin = interp1(x,Q,xMin,'spline');
                sumMax = interp1(x,Q,xMax,'spline');
                result(i) = {sumMax-sumMin};
            end
        end
 
        % Find the average value of a curve over a given domain
        function avgVal = getAvgVal(device, yParam, xParam, xMin, xMax)
            integral = Transistor.integrateCurve(device, yParam, xParam, xMin, xMax);
            avgVal = {};
            for i = 1:numel(integral)
                avgVal(i) = {integral{i}/(xMax-xMin)};
            end
        end
            
        % Find the charge in the output capacitance at a given Vds
            % If Vds is greater than Vds_Max, throws an error
        function Qoss = getChargeFromCoss(device, Vds)
            % Check for positive Vds
            if Vds <= 0
                error('Vds must be a positive numeric value.')
            end
            
            % Check for Vds < Vds Max
            vds_max = Transistor.getTableParam(device,'Vds','max');
            if isnan(vds_max)
                disp(['WARNING: No Vds max for ' device '. Data may be inaccurate']);
            elseif Vds > vds_max
                error(['Given Vds, ' num2str(Vds) '(V), is greater than the max voltage '...
                    'for this device: ' num2str(vds_max) '(V), ' device]);
            end

            % Get Qoss data
            Qoss = Transistor.integrateCurve(device, 'Coss','Vds',0,Vds);
            
        end
           
        % Returns a cell array of devices whom have a maximum parameter value greater than the given threshold
            % Uses SQL filtering
        function validDevices = filter(param, min_typ_max, operand, threshold)
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            % Check inputs
            if ~contains({'min','typ','max'},lower(strip(min_typ_max)))
                error('Second argument must be Min Typ or Max');
            end
            if ~contains({'=','>','<','>=','<='},strip(operand))
                error(['Third argument must be one of the following' newline '= > < >= <=']);
            end
            min_typ_max = strip(min_typ_max);
            min_typ_max = [upper(min_typ_max(1)) lower(min_typ_max(2:end))];
            
            % Connect to SQL DB
            if ~isfile([path 'TransistorDB\Transistor_StandardTable.db'])
                error(['Relational database file ' path 'TransistorDB\Transistor_StandardTable.db does not exist.']);
            end
            conn = sqlite([path '\TransistorDB\Transistor_StandardTable.db']);
            
            % Find valid devices
            validDevices = cell(0);
            try
                [~,~,tableFields,~,tableDefaultMultipliers,~,~,~,~,~,~] = Transistor.loadMetaData(); 
                index = find(strcmpi(tableFields,strip(param)));
                if isempty(index)
                    error('no such column');
                end
                param = tableFields{index};
                mult = Transistor.SIprefixes(tableDefaultMultipliers{strcmp(tableFields,param)});
                threshold = threshold/mult;
                
                query = ['SELECT "partNumber" FROM Transistor WHERE "' ... 
                    param '_' min_typ_max '" ' operand ' ' num2str(threshold)];
                validDevices = fetch(conn, query);
            catch Exception
                if contains(Exception.message,'no such column')
                    error(['Param ' strip(param) ' does not exist in the DB']);
                elseif contains(Exception.message,'Unexpected NULL')
                    % No data, return NaN
                else
                    rethrow(Exception)
                end
            end
        end

        
        
        %%% Other misc. static public methods
        % Delete a part
        function removeFromDatabase(partNumber)
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            % Connect to SQL DB
            rldb_file = [path 'TransistorDB\Transistor_StandardTable.db'];
            if ~isfile(rldb_file)
                error(['Relational database file ' rldb_file ' does not exist.']);
            end
            try
                conn = sqlite(rldb_file);
                exec(conn, ['DELETE FROM Transistor WHERE partNumber="' partNumber '"']);
                close(conn)
            catch 
            end
            
            % Connect to additional table parameters db
            nrdb_file = [path 'TransistorDB\Transistor_AdditionalTable.json'];
            if ~isfile(nrdb_file)
                error(['Non-relational database file ' nrdb_file ' does not exist.']);
            end
            try 
                paramsStruct = jsondecode(fileread(nrdb_file));
                paramsStruct = rmfield(paramsStruct,partNumber);
                fid = fopen(nrdb_file, 'w');
                fprintf(fid,jsonencode(paramsStruct));
                fclose(fid);
            catch
            end

            % Connect to graph parameters db
            nrdb_file = [path 'TransistorDB\Transistor_Graph.json'];
            if ~isfile(nrdb_file)
                error(['Non-relational database file ' nrdb_file ' does not exist.']);
            end
            try 
                paramsStruct = jsondecode(fileread(nrdb_file));
                paramsStruct = rmfield(paramsStruct,partNumber);
                fid = fopen(nrdb_file, 'w');
                fprintf(fid,jsonencode(paramsStruct));
                fclose(fid);
            catch
            end
        end 

        % Push (export) DBs to another folder
        function push(other_location, userID)     
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            msg = {};
            time = char(datetime);
            time = strrep(time,':','-');
            time = strrep(time,' ','_');
            fid = fopen([path 'TransistorDB/Log.txt'], 'a+');
            fprintf(fid, '%s,%s', time, userID);
            
            % Push Standard Table Parameters
            if ~isfile([other_location '\Transistor_StandardTable.db'])
                copyfile([path 'TransistorDB\Transistor_StandardTable.db'], ...
                    [other_location '\Transistor_StandardTable.db']);
                msg = [msg ; 'Standard table parameter database created.'];
            else
                try
                    changedStr = Transistor.merge_tableParameters([other_location '\Transistor_StandardTable.db'], ...
                        [path 'TransistorDB\Transistor_StandardTable.db']);
                    if ~isempty(changedStr)
                        fprintf(fid, '%s', changedStr);
                    end
                    msg = [msg ; 'Standard table parameters pushed.' newline];
                catch Error
                    msg = [msg ; 'WARNING: Standard table parameters NOT pushed:' ; Error.message newline];
                end
            end
            
            % Push Additional Table Parameters
            if ~isfile([other_location '\Transistor_AdditionalTable.json'])
                copyfile([path 'TransistorDB\Transistor_AdditionalTable.json'], ...
                    [other_location '\Transistor_AdditionalTable.json']);
                msg = [msg ; 'Additional table parameter database created.'];
            else
                try
                    Transistor.merge_additionalTable([other_location '\Transistor_AdditionalTable.json'], ...
                        [path 'TransistorDB\Transistor_AdditionalTable.json']);
                    msg = [msg ; 'Additional table parameters pushed.' newline];
                catch Error
                    msg = [msg ; 'WARNING: Additional table parameters NOT pushed:' ; Error.message newline];
                end
            end
            
            % Push Graph Parameters
            if ~isfile([other_location '\Transistor_Graph.json'])
                copyfile([path 'TransistorDB\Transistor_Graph.json'], ...
                    [other_location '\Transistor_Graph.json']);
                msg = [msg ; 'Graph parameters database created.'];
            else
                try
                    Transistor.merge_graphParameters([other_location '\Transistor_Graph.json'], ...
                        [path 'TransistorDB\Transistor_Graph.json'], 'PUSH');
                    msg = [msg ; 'Graph parameters pushed.' newline];
                catch Error
                    msg = [msg ; 'WARNING: Graph parameters NOT pushed:' ; Error.message newline];
                end
            end            
            
            fprintf(fid, newline);
            fclose(fid);
            copyfile([path 'TransistorDB/Log.txt'], [other_location '/Log.txt']);
            
            msgbox(msg, 'Push Finished');            
        end
        
        % Pull (import) DBs from another folder
        function pull(other_location) 
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            msg = {};
            copyfile([other_location '/Log.txt'], [path 'TransistorDB/Log.txt']);
            
            try
                Transistor.merge_tableParameters([path 'TransistorDB\Transistor_StandardTable.db'], ...
                    [other_location '\Transistor_StandardTable.db']);
                msg = [msg ; 'Standard table parameters pulled.' newline];
            catch Error
                msg = [msg ; 'WARNING: Standard table parameters NOT pulled:' ; Error.message newline];
            end
            
            try
                Transistor.merge_additionalTable([path 'TransistorDB\Transistor_AdditionalTable.json'], ...
                    [other_location '\Transistor_AdditionalTable.json']);
                msg = [msg ; 'Additional table parameters pulled.' newline];
            catch Error
                msg = [msg ; 'WARNING: Additional table parameters NOT pulled:' ; Error.message newline];
            end
            
            try
                Transistor.merge_graphParameters([path 'TransistorDB\Transistor_Graph.json'], ...
                        [other_location '\Transistor_Graph.json'], 'PULL');
                msg = [msg ; 'Graph parameters pulled.' newline];                    
            catch Error
                msg = [msg ; 'WARNING: Graph parameters NOT pulled:' ; Error.message newline];
            end
            
            msgbox(msg, 'Pull Finished');
        end
        
        % List all devices
        function devices = listDevices()
            path_parts = strsplit(which('Transistor'),'Transistor.m');
            path = path_parts{1};
            
            devices = {};
            % Connect to SQL DB
            rldb_file = [path 'TransistorDB\Transistor_StandardTable.db'];
            if ~isfile(rldb_file)
                return
            end
            try
                conn = sqlite(rldb_file);
                devices = fetch(conn, 'SELECT partNumber FROM Transistor');
                close(conn)
            catch 
            end
        end
    end
    
    
end


