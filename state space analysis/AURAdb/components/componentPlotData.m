classdef componentPlotData 
    %componentPlotData Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        plotData    % x-y data extracted from plot
        testConditions
        dataLabels  % labels for individual traces
        title
    end
    
    properties (Dependent)
       xLabel
       yLabel
       zLabel
    end
    
    properties (Access = protected)
        componentType
    end
    
    properties (Hidden)
        axisLabels  % string axislabel
        axisType = {'lin', 'lin'}    % two-element {x,y}, values 'log', 'norm', or 'lin'
        SIUnits
        conversionEps = 1e-8; %relative change in data acceptable from sync
        dim = 2
    end
    
    properties (Dependent, Hidden)
        nTraces
        type
       
        xSignal
        ySignal 
        zSignal
        
        xUnit
        yUnit
        zUnit
        
        xSI
        ySI
        zSI
        
        xMult
        yMult
        zMult
    end
        
    
    methods
        function obj = componentPlotData(type, varargin)
            %componentPlotData Construct an instance of this class
            %   Detailed explanation goes here
            
            assert( isa(type, 'component'), 'input ''type'' must be a component class or subclass of component');
            obj.componentType = type;
            
            if nargin == 2 && isa(varargin{1}, 'digitizedPlot')
%                 obj = obj.copyDigitizedPlot(varargin{1});
                dP = varargin{1};
                obj.plotData = dP.plotData;
                obj.SIUnits  = dP.SIUnits;
                
                %%Convert to base units
                for i = 1:length(obj.plotData)
                    obj.plotData{i}(:,1) = obj.plotData{i}(:,1)*obj.componentType.SIprefixes(obj.SIUnits{1});
                    obj.plotData{i}(:,2) = obj.plotData{i}(:,2)*obj.componentType.SIprefixes(obj.SIUnits{2});
                end
                obj.SIUnits = {'1', '1'};
                
                for i = 1:2
                    if dP.logAxes(i)
                        obj.axisType{i} = 'log';
                    elseif dP.normAxes(i)
                        obj.axisType{i} = 'norm';
                    else
                        obj.axisType{i} = 'lin';
                    end
                end
                obj.axisLabels  = dP.axisLabels;
                obj.dataLabels  = dP.dataLabels;
                obj.testConditions = dP.testConditions;
                obj.title = dP.getPlotTitle('short');
            elseif nargin == 2 && isa(varargin{1}, 'componentPlotData')
                dP = varargin{1};
                obj.plotData = dP.plotData;
                obj.SIUnits  = dP.SIUnits;
                obj.axisType = dP.axisType;
                obj.axisLabels  = dP.axisLabels;
                obj.dataLabels  = dP.dataLabels;
                obj.testConditions = dP.testConditions;
                obj.title = dP.title;
            end

        end
        
        function plot(obj, axes)

           if nargin < 2
               axes = gca;
           end
           hold(axes, 'off');

          
                     
           xM = obj.componentType.defaultMultipliers(strcmpi(isParamOf(obj.componentType,obj.xSignal), obj.componentType.knownParams));
           xMn = obj.componentType.SIprefixes(xM{1});
           yM = obj.componentType.defaultMultipliers(strcmpi(isParamOf(obj.componentType,obj.ySignal),obj.componentType.knownParams));
           if isempty(yM)
               yM = obj.componentType.defaultMultipliers(strcmpi(obj.dataLabels{1}, obj.componentType.knownParams));
           end
           yMn = obj.componentType.SIprefixes(yM{1});
           
           for i=1:obj.nTraces

               if strcmp(obj.axisType{1}, 'log')
                   xData = abs(obj.plotData{i}(:,1)/xMn);
               else
                   xData = obj.plotData{i}(:,1)/xMn;
               end

               if strcmp(obj.axisType{2}, 'log')
                   yData = abs(obj.plotData{i}(:,2)/yMn);
               else
                   yData = obj.plotData{i}(:,2)/yMn;
               end

               if size(obj.plotData{i},2) == 2
                    plot(axes, xData, yData, '-o');
               else
                   zM = obj.componentType.defaultMultipliers(strcmpi(isParamOf(obj.componentType,obj.zSignal), obj.componentType.knownParams));
                   zMn = obj.componentType.SIprefixes(zM{1});
                   if numel(obj.axisType) >= 3 && strcmp(obj.axisType{3}, 'log')
                       zData = abs(obj.plotData{i}(:,3)/zMn);
                   else
                       zData = obj.plotData{i}(:,3)/zMn;
                   end
                   xAx = unique(xData);
                   yAx = unique(yData);

                   [Xg,Yg] = meshgrid(xAx,yAx);
                   F = scatteredInterpolant(xData,yData,zData);
                   Zg = F(Xg,Yg);

                   surf(axes, Xg, Yg, Zg);
               end
               hold(axes, 'on');
           end
           legend(axes, obj.dataLabels);
           
           if ~isempty(obj.testConditions)
               h = line(axes, nan, nan, 'Color', 'none');
               entry = '';
               for i = 1:size(obj.testConditions, 1)
                    unit = obj.testConditions{i,4};
                    if strcmp(unit, 'none')
                       unit = '';
                    end
                    entry = [entry obj.testConditions{i,2}, '=', num2str(obj.testConditions{i,3}), unit, obj.testConditions{i,5} newline];
               end
               legend(axes, {obj.dataLabels{:}, entry},'location','best');
           end
           
          if strcmp(obj.axisType{1}, 'log')
               axes.XScale = 'log';
               grid(axes, 'on');
           end
           
           if strcmp(obj.axisType{2}, 'log')
               axes.YScale = 'log';
               grid(axes, 'on');
           end
           
           xlabel(axes, [obj.xLabel, ' [', xM{1}, obj.xUnit, ']']);
           ylabel(axes, [obj.yLabel, ' [', yM{1}, obj.yUnit, ']']);

           if size(obj.plotData{1},2) >= 3
               zlabel(axes, [obj.zLabel, ' [', zM{1}, obj.zUnit, ']']);
           end
           
        end

        function obj = deleteTrace(obj,ind)
            assert(obj.nTraces > length(ind), 'Cannot delete all traces from plot')
            remainingCurves = setdiff(1:length(obj.plotData), ind);
            obj.plotData = {obj.plotData{remainingCurves}};
            obj.dataLabels = {obj.dataLabels{remainingCurves}};
            obj.title = obj.title(1:strfind(obj.title, "(")-2);
            if length(obj.dataLabels) > 1
                obj.title = [obj.title ' (', strjoin(obj.dataLabels, ', '), ')'];
            end

        end
        
        function cT = get.type(obj)
            cT = class(obj.componentType);
        end
        
        function xs = get.xLabel(obj)
            xs = obj.axisLabels{1};
        end
        
        function ys = get.yLabel(obj)
            ys = obj.axisLabels{2};
        end

        function zs = get.zLabel(obj)
            if numel(obj.axisLabels) >= 3
                zs = obj.axisLabels{3};
            else
                zs = [];
            end
        end
        
        function xs = get.xSignal(obj)
            indX = strfind(obj.xLabel, ',');
            if ~isempty(indX)
                xs = obj.xLabel(1:indX-1);
            else
                xs = obj.xLabel;
            end
        end
        
        function ys = get.ySignal(obj)
            indY = strfind(obj.yLabel, ',');
            if ~isempty(indY)
                ys = obj.yLabel(1:indY-1);
            else
               ys = obj.yLabel;
            end
        end

         function zs = get.zSignal(obj)
            indZ = strfind(obj.zLabel, ',');
            if ~isempty(indZ)
                zs = obj.zLabel(1:indZ-1);
            else
               zs = obj.zLabel;
            end
        end
        
        function nT = get.nTraces(obj)
            nT = size(obj.plotData,2);
        end
        
        function xU = get.xUnit(obj)
            xU = obj.componentType.defaultUnits(strcmpi(obj.componentType.knownParams, obj.xSignal));
            xU = xU{1};
        end

        function yU = get.yUnit(obj)
            yU = obj.componentType.defaultUnits(strcmpi(obj.componentType.knownParams, obj.ySignal));
            if isempty(yU)
               yU = obj.componentType.defaultUnits(strcmpi(obj.componentType.knownParams, obj.dataLabels{1}));
            end
            yU = yU{1};
        end    

        function zU = get.zUnit(obj)
            zU = obj.componentType.defaultUnits(strcmpi(obj.componentType.knownParams, obj.zSignal));
            if isempty(zU)
               zU = obj.componentType.defaultUnits(strcmpi(obj.componentType.knownParams, obj.dataLabels{1}));
            end
            zU = zU{1};
        end  
        
        function xU = get.xSI(obj)
            xU = obj.SIUnits{1}; 
        end
        
        function yU = get.ySI(obj)
            yU = obj.SIUnits{2}; 
        end

        function zU = get.zSI(obj)
            if numel(obj.SIUnits) >= 3
                zU = obj.SIUnits{3}; 
            else
                zU = '1';
            end
        end
        
        function xM = get.xMult(obj)
            if strcmp(obj.xSI, 'none') || strcmpi(obj.xSI, '1')
               xM = 1;
            else
               xM = obj.componentType.SIprefixes(obj.xSI);
            end
        end

        function yM = get.yMult(obj)
            if strcmp(obj.ySI, 'none') || strcmpi(obj.ySI, '1')
                yM = 1;
            else
                yM = obj.componentType.SIprefixes(obj.ySI);
            end
        end

        function zM = get.zMult(obj)
            if strcmp(obj.zSI, 'none') || strcmpi(obj.zSI, '1')
                zM = 1;
            else
                zM = obj.componentType.SIprefixes(obj.zSI);
            end
        end

    
        function [exact, plotTF] = eq(obj, graph)
            % [tf, plotTF] = eq(obj, graph)
            % tf is true if the two graphs are identical
            % plotTF is true if the two are the same plot, but with
            % different data
            if length(obj) > 1 && length(graph) == 1
                for i = 1:length(obj)
                    [xtf, xplotTF] = eq(obj(i), graph);
                    exact(i) = xtf;
                    plotTF(i) = xplotTF;
                end
                return
            elseif length(obj) == 1 && length(graph) > 1
                for i = 1:length(graph)
                    [xtf, xplotTF] = eq(obj, graph(i));
                    exact(i) = xtf;
                    plotTF(i) = xplotTF;
                end
                return
            elseif length(obj) > 1 && length(graph) > 1
                error('eq() is only defined when one object is singleton');
            else
                if ~isa(graph, class(obj))
                    plotTF = false; 
%                 elseif length(obj.plotData) ~= length(graph.plotData)
%                     tf = false; 
                elseif ~strcmp(obj.xSignal, graph.xSignal)
                    plotTF = false; 
                elseif ~strcmp(obj.ySignal, graph.ySignal)
                    plotTF = false; 
                elseif ~strcmp(obj.title, graph.title)
                    plotTF = false; 
                %elseif ~all(strcmp(obj.axisType, graph.axisType))
                    %tf = false; 
%                 elseif ~all(strcmp(obj.dataLabels, graph.dataLabels))
%                     tf = false; 
                else
                    plotTF = true;
                end

                if ~plotTF
                    exact = false; return
                else
                    exact = true;
                    for i = 1:length(obj.plotData)
                        % Due to storing online and returning, some
                        % resolution may be lost, so check for less then
                        % 1e-6 % difference in values
                        exact =  all(size(obj.plotData) == size(graph.plotData)) && ...
                            all(cellfun(@isequal,cellfun(@size,obj.plotData,'UniformOutput',false),cellfun(@size,graph.plotData,'UniformOutput',false)) )&& ...
                            all([exact ...
                            all(abs(obj.plotData{i} - graph.plotData{i})./max(abs(graph.plotData{i}),eps) < obj.conversionEps | ...
                                (isnan(obj.plotData{i}) & isnan(graph.plotData{i}))) ...
                            ]);
                    end
                end
            end
        end
        
        function merge(obj, newGraph, mode)
            if nargin ==2
                mode = 0;
            end
            
            %mode == 0 --> graphical
            %mode == 1 --> replace
            %mode == 2 --> keep old data
            %mode == 3 --> merge
            
            assert(isa(newGraph, class(obj)), ...
                ['merge can only be used on objects of the same type'] );
            [exact, plotTF] = obj.eq(newGraph);
            assert(any(plotTF), 'merge can only be used on two plots of the same type');
            if any(exact)
                warning('Merge should not be called on two identical copies');
                return
            end
            
            % so it is definitely same plot, different data
            for i = 1:length(obj.plotData)
                for j = 1:length(newGraph.plotData)
                    minL = min(length(obj.plotData{i}), length(newGraph.plotData{j}));
                    exact(i,j) =  all(abs(obj.plotData{i}(1:minL,:) - newGraph.plotData{j}(1:minL,:))./max(abs(newGraph.plotData{j}(1:minL,:)),eps) < obj.conversionEps, 1:2);
                end
            end
            newCurves = sum(exact,1) ~=1;
            
            %create a copy of the object for editing 
            %(NOT A HANDLE CLASS, so this was unecessary      
            mergedGraph = componentPlotData(obj.componentType, obj);
            
            for i = find(newCurves)
                sameLabel = find(strcmp(newGraph.dataLabels, obj.dataLabels(i)));
                if sameLabel
                   %We already have the curve, but the new graph has different data
                   mergedGraph.plotData{sameLabel} = unique(sortrows([obj.plotData{sameLabel}; newGraph.plotData{i}], 1), 'rows');
                else
                    % The label is different, so it may be a new curve on
                    % the same plot
                     mergedGraph.plotData{end+1} = newGraph.plotData{i};
                     mergedGraph.dataLabels{end+1} =  newGraph.dataLabels{end+1};
                end
                mergedTestConditions = [obj.testConditions;newGraph.testConditions];
                if ~isempty(mergedTestConditions)
                    mergedGraph.testConditions = table2cell(unique(cell2table(mergedTestConditions)));
                end
            end
            
            if mode == 1
                obj = newGraph;
            elseif mode == 2
                return
            elseif mode == 3
                obj = mergedGraph;
            elseif mode == 0
                figure;
                subplot(3,1,1)
                plot(obj, gca);
                subplot(3,1,2)
                plot(mergedGraph, gca);
                subplot(3,1,3)
                plot(newGraph, gca);
                while(1)
                    uI = input(['Conflicting data found for plot title: ' newGraph.title '\n' ...
                        'Top plot is old data, middle plot is merged data, bottom plot is new data' ... 
                        '\n Select Top [T], Middle [M], or Bottom [B] plot to save \n'], 's');
                    if strcmpi(uI,'T')
                        break
                    elseif strcmpi(uI,'B')
                        obj = mergedGraph;
                        break
                    elseif strcmpi(uI,'M')
                        obj = newGraph;
                        break
                    else
                        display('Invalid Input');
                    end     
                end
                close(gcf);
            end
        end
    end
end

