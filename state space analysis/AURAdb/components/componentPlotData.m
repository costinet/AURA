classdef componentPlotData 
    %componentPlotData Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        plotData    % x-y data extracted from plot
        testConditions   
        title
    end
    
    properties (Dependent)
       xLabel
       yLabel
    end
    
    properties (Access = protected)
        componentType
    end
    
    properties (Hidden)
        axisLabels  % string axislabel
        dataLabels  % labels for individual traces
        axisType    % two-element {x,y}, values 'log', 'norm', or 'lin'
        SIUnits
    end
    
    properties (Dependent, Hidden)
        nTraces
        type
       
        xSignal
        ySignal 
        
        xUnit
        yUnit
        
        xSI
        ySI
        
        xMult
        yMult
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
            end

        end
        
        function plot(obj, axes)
           hold(axes, 'off');
           
           
           xM = obj.componentType.defaultMultipliers(strcmp(obj.xSignal, obj.componentType.knownParams));
           xMn = obj.componentType.SIprefixes(xM{1});
           yM = obj.componentType.defaultMultipliers(strcmp(obj.ySignal,obj.componentType.knownParams));
           if isempty(yM)
               yM = obj.componentType.defaultMultipliers(strcmp(obj.dataLabels{1}, obj.componentType.knownParams));
           end
           yMn = obj.componentType.SIprefixes(yM{1});
           
           for i=1:obj.nTraces

               plot(axes, obj.plotData{i}(:,1)/xMn, obj.plotData{i}(:,2)/yMn);
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
               legend(axes, {obj.dataLabels{:}, entry});
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
        
        function nT = get.nTraces(obj)
            nT = size(obj.plotData,2);
        end
        
        function xU = get.xUnit(obj)
            xU = obj.componentType.defaultUnits(strcmp(obj.componentType.knownParams, obj.xSignal));
            xU = xU{1};
        end

        function yU = get.yUnit(obj)
            yU = obj.componentType.defaultUnits(strcmp(obj.componentType.knownParams, obj.ySignal));
            if isempty(yU)
               yU = obj.componentType.defaultUnits(strcmp(obj.componentType.knownParams, obj.dataLabels{1}));
            end
            yU = yU{1};
        end        
        
        function xU = get.xSI(obj)
            xU = obj.SIUnits{1}; 
        end
        
        function yU = get.ySI(obj)
            yU = obj.SIUnits{2}; 
        end
        
        function xM = get.xMult(obj)
            if strcmp(obj.xSI, 'none') || strcmp(obj.xSI, '1')
               xM = 1;
            else
               xM = obj.componentType.SIprefixes(obj.xSI);
            end
        end

        function yM = get.yMult(obj)
            if strcmp(obj.ySI, 'none') || strcmp(obj.ySI, '1')
                yM = 1;
            else
                yM = obj.componentType.SIprefixes(obj.ySI);
            end
        end

    
        function [tf, plotTF] = eq(obj, graph)
            % [tf, plotTF] = eq(obj, graph)
            % tf is true if the two grpahs are identical
            % plotTF is true if the two are the same plot, but with
            % different data
            if length(obj) > 1 && length(graph) == 1
                for i = 1:length(obj)
                    [xtf, xplotTF] = eq(obj(i), graph);
                    tf(i) = xtf;
                    plotTF(i) = xplotTF;
                end
                return
            elseif length(obj) == 1 && length(graph) > 1
                for i = 1:length(graph)
                    [xtf, xplotTF] = eq(obj, graph(i));
                    tf(i) = xtf;
                    plotTF(i) = xplotTF;
                end
                return
            elseif length(obj) > 1 && length(graph) > 1
                error('eq() is only defined when one object is singleton');
            else
                if ~isa(graph, class(obj))
                    tf = false; 
                elseif length(obj.plotData) ~= length(graph.plotData)
                    tf = false; 
                elseif ~strcmp(obj.xSignal, graph.xSignal)
                    tf = false; 
                elseif ~strcmp(obj.ySignal, graph.ySignal)
                    tf = false; 
                elseif ~all(strcmp(obj.axisType, graph.axisType))
                    tf = false; 
                elseif ~all(strcmp(obj.dataLabels, graph.dataLabels))
                    tf = false; 
                else
                    tf = true;
                end

                if ~tf
                    plotTF = false; return
                else
                    plotTF = true;
                    for i = 1:length(obj.plotData)
                        tf = all([tf, all(obj.plotData{i} == graph.plotData{i})]);
                    end
                end
            end
        end
    end
end

