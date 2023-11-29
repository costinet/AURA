classdef transistorDB < componentDB
    %transistorDB is a database of transistor objects
    
%     Inherited properties:
%         tableProps
%         addlTabProps
%         graphProps
%         components

    properties (Dependent)
        transistors
    end
    
    properties (Hidden, Transient, Constant)
        componentType = transistor();
    end
    
    % componentDB required methods
    methods 
        [table, addlTable, graph] = json(obj)
        subTable(obj, varargin) % Return a nicely-formatted table with searchable parameters
    end
    
    methods (Hidden)

    end
    
    methods
        function obj = transistorDB()
            loadDB(obj)
        end
        
        function transistors = get.transistors(obj)
            transistors = obj.components;
        end    
        
        function saveDB(obj)
            fn = mfilename('fullpath');
            fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
            save(fn,'obj')
        end
        
        function loadDB(obj)
            fn = mfilename('fullpath');
            fn = [strrep(fn, '\transistorDB', '\') 'transistors.mat'];
            if isfile(fn)
                savedData = load(fn,'obj');
                obj.components = savedData.obj.transistors;
            end
        end

        function clearUpdated(obj)
            for i = 1:length(obj.components)
                obj.components(i).clearUpdated();
            end
        end

        function digitizeDatasheet(obj,datasheet)
            %% 
            plotDigitizer(datasheet);
        end
        
        function sync(obj)
            % Upload all transistors
            host = netConnect();
            response = host.postData(obj.transistors)
            
            [graphs, params] = host.pullData;
            
            for i = 1:length(graphs)
%                 graphs(i).plotData = jsondecode(graphs(i).plotData);
%                 graphs(i).testConditions = jsondecode(graphs(i).testConditions);
                
                cPD = componentPlotData(transistor);
                
                
                % jsondecode sometimes recognizes cells, but will return it as a 3D array
                % if all curves have the same number of points, so this will fix it to always
                % be cell arrays
                plotData = jsondecode(graphs(i).plotData);
                if isa(plotData, 'cell')
                    cPD.plotData = plotData';
                elseif isa(plotData, 'double')
                    if size(plotData,3) == 2
                        for j=1:size(plotData,1)
                            cPD.plotData{j} = [plotData(j,:,1)', plotData(j,:,2)'];
                        end
                    else
                        error('test2');
                    end
                else
                    error('test');
                end
                
                % also gets wonky -- need to reshape
                testConditions = jsondecode(graphs(i).testConditions)';
                cPD.testConditions = reshape(testConditions, [length(testConditions)/5, 5]);
                
                cPD.dataLabels = jsondecode(graphs(i).dataLabels);
                cPD.title = graphs(i).title;
                cPD.axisLabels{1} = graphs(i).xLabel;
                cPD.axisLabels{2} = graphs(i).yLabel;  
                
                loc = find(strcmp({obj.transistors.partNumber}, graphs(i).partNumber),1);
%                 if isempty(loc)
%                     loc = length(obj.transistors)+1;
%                 end

                if isempty(loc) %New part Number
                    loc = length(obj.transistors)+1;
                    obj.components(loc) = transistor(graphs(i).partNumber);
                    obj.components(loc).addGraph(cPD);
                else
                    obj.components(loc).addGraph(cPD);
                end
 
            end
            
            for i = 1:length(params)
                warning('ParamSync: you havent written this yet')
                %% REMEMBER component.addParameter(obj,param) -- componentTableData
                cTD = componentTableData(transistor)
%                 componentTableData(type, paramName,typVal,maxVal,minVal,testConditions, units)
                
            end

            obj.saveDB;
            obj.clearUpdated
        end

    end
    
end

