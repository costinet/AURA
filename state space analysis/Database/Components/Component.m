classdef Component < handle
    
    properties (SetAccess = protected)
        componentType;
        partNumber = '';
        
        tableParameters = struct();
        tableFields;
        tableUnits;
        tableDefaultMultipliers;
        
        additionalTableParameters = struct();
        
        graphParameters = struct();
        graphFields;
        graphUnits;
        graphDefaultMultipliers;
        
        conditionFields;
        conditionUnits;
        conditionDefaultMultipliers;
    end
    
    properties (Access = protected)
        temp_graphParameters = struct();
        current_graphParameter = '';
        current_graphPage = 1;
    end
    
    properties (Access = protected, Constant)
        SIkeys = {'f','p','n','u','m','1','k','M','G','T','P'};
        SIprefixes = containers.Map({'f','p','n','u','m','1','k','M','G','T','P'}, ...
                                    [1E-15 1E-12 1E-9 1E-6 1E-3 1E0 1E3 1E6 1E9 1E12 1E15]);        
    end
    
    
    methods (Access = protected)
        % Plot data, given handles to GUI object with plot information
        function updateBtnPushed(obj, dataBox, dataPlot, y_fieldBox, x_fieldBox, y_logBox, x_logBox, y_multiplierBox, x_multiplierBox)
            try
                X = Component.dataStrToMat(dataBox);
                x = X(:,1);
                y = X(:,2);
                if length(x) == 1 && x(1) == 0
                    dataBox.Value = '';
                else
                    try
                        dataPlot.XLim = [min(x), max(x)]; 
                    catch
                        dataPlot.XLim = [-inf inf];
                    end
                end
                plot(dataPlot,x,y,'-o','LineWidth',2);
                if strcmp(x_logBox.Value,'Linear') 
                    dataPlot.XScale = 'linear';
                else
                    dataPlot.XScale = 'log';
                end
                if strcmp(y_logBox.Value,'Linear') 
                    dataPlot.YScale = 'linear';
                else
                    dataPlot.YScale = 'log';
                end
                
                if ~strcmp(x_fieldBox.Value,'--Select X-axis Field--') && ~strcmp(y_fieldBox.Value,'--Select Y-axis Field--')
                    xField = strsplit(x_fieldBox.Value(1:end-1), '(');
                    yField = strsplit(y_fieldBox.Value(1:end-1), '(');            
                    xMultiplier = x_multiplierBox.Value;
                    yMultiplier = y_multiplierBox.Value;
                    if strcmp(xMultiplier, '1')
                        xMultiplier = '';
                    end
                    if strcmp(yMultiplier, '1')
                        yMultiplier = '';
                    end

                    xlabel = [strip(xField{1}) ' (' xMultiplier strip(xField{2}) ')'];
                    ylabel = [strip(yField{1}) ' (' yMultiplier strip(yField{2}) ')'];
                    title_str = [ylabel ' vs. ' xlabel ', Page ' num2str(obj.current_graphPage)];
                    title(dataPlot, strrep(title_str, '_', '\_'), 'FontSize', 16);
                else
                    title(dataPlot, 'Select Fields for X and Y axes', 'FontSize', 16)
                end
                
            catch Error
%                 rethrow(Error);
                Component.dialogBox({'Invalid data points, use default data output from WebPlotDigitizer'}, 12);
            end
        end          
    end
    
    methods (Access = protected, Static)
        % For graph parameters GUI- Convert formatted string for text area to data matrix        
        function X = dataStrToMat(dataBox)
            line_cells = dataBox.Value;
            non_empty_cell_indices = ~cellfun(@isempty, dataBox.Value);
            non_empty_cells = line_cells(non_empty_cell_indices);
            text = strjoin(non_empty_cells, ';');
            rows = strsplit(text,';');
            
            % If there is no data, return a (0,0) point
            if numel(rows) == 1 && isempty(rows{1})
                X = [0, 0];
                return
            end
            
            % Get data
            X = zeros(numel(rows),2);
            for i = 1:numel(rows)
                split_str = strsplit(rows{i},',');
                X(i,1) = str2double(strip(split_str{1}));
                X(i,2) = str2double(strip(split_str{2}));
            end    
        end
        
        % For graph parameters GUI- Convert data matrix to formatted string for text area
        function S = dataMatToStr(X)
            data_str = strrep(mat2str(X),';',newline);
            S = strrep(data_str(2:end-1),' ',', ');
        end
        
        function deviceNumber = checkDuplicateDevices(givenDevice, devices)
            deviceNumber = givenDevice;
            
            % Exact match, use previous data
            if any(strcmp(devices,givenDevice))
                deviceNumber = givenDevice;
                return
            end

            while length(givenDevice) > 3 && strcmpi(givenDevice(end-2:end),'-nd')
                answer = inputdlg('Please enter the Manufacturer Part Number, not the Digi-Key Part Number:');
                givenDevice = strip(answer{1});
            end
                
            % Check for near matches
            % First, check if given device number is substring of an existing device number
            for i = 1:numel(devices)
                device = devices{i};
                if contains(device, givenDevice)
                    answer = questdlg(['Your given device number, ' givenDevice ', is similar to an existing device: ' ...
                        device '. Double-check that the device you are entering is not the same device.'], 'Similar Device', ...
                        'It is the same device.', 'It is a new device.', 'It is the same device.');
                    if strcmp(answer, 'It is the same device.')
                        deviceNumber = device;
                    else
                        deviceNumber = givenDevice;
                    end
                    return
                end
            end
            % Next, check if given device number is a superstring of an existing device number
            for i = 1:numel(devices)
                device = devices{i};
                if contains(givenDevice, device)
                    answer = questdlg(['Your given device number, ' givenDevice ', is very similar to an existing device: ' ...
                        device '. Double-check that the device you are entering is not the same device.'], 'Similar Device', ...
                        'It is the same device.', 'It is a new device.', 'It is the same device.');
                    if strcmp(answer, 'It is the same device.')
                        deviceNumber = device;
                    else
                        deviceNumber = givenDevice;
                    end
                    return
                end
            end
            % Lastly, check for close similarity betweeen given device number and an existing device using an algorithmic function 
            for i = 1:numel(devices)
                device = devices{i};
                if EditDist(device,givenDevice) <= 2
                    answer = questdlg(['Your given device number, ' givenDevice ', is very similar to an existing device: ' ...
                        device '. Double-check that the device you are entering is not the same device.'], 'Similar Device', ...
                        'It is the same device.', 'It is a new device.', 'It is the same device.');
                    if strcmp(answer, 'It is the same device.')
                        deviceNumber = device;
                    else
                        deviceNumber = givenDevice;
                    end
                    return
                end
            end
        end
       
    end
    
    methods (Abstract)
        datasheet(obj);
        loadData(obj);
        saveData(obj);   
    end
    
    methods (Abstract, Access = protected)
        save_rldb(obj, rldb_file);
        save_nrdb(obj, nrdb_file);
        load_rldb(obj, rldb_file);
        load_nrdb(obj, nrdb_file);
    end
    
    methods (Abstract, Access = protected, Static)
        changedStr = merge_tableParameters(rldb_file, other_rldb_file);
        merge_additionalTable(nrdb_file, other_nrdb_file);
        merge_graphParameters(nrdb_file, other_nrdb_file);
    end
    
    methods (Abstract, Static)
        removeFromDatabase(partNumber); 
        push(other_location);
        pull(other_location);
        devices = listDevices(obj);
    end

end

