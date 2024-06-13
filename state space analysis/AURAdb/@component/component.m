classdef component < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=protected)
        parameters = componentTableData.empty;
        graphs = componentPlotData.empty;
    end
    
    properties (Abstract)
        partNumber
        manufacturer
        material
        type
    end
    
    properties (Hidden, Transient, Abstract, Constant)
%         tableParams
        knownParams
        paramDict
        defaultUnits
        defaultMultipliers
        paramNames
        
        knownTypes
        knownMaterials
    end
    
    properties (Hidden, Transient, Constant)
        signalTypes = {'Capacitance', 'Voltage', 'Current', 'Resistance', 'Energy', 'Charge', 'Temperature', 'Thermal Resistance', 'Time', 'Transconductance', 'Inductance', 'Flux Density'};
        signalUnits = {'F','V','A','Ohm','J','C',[char(176) 'C'],[char(176) 'C/W'], 's', 'A/V', 'H', 'T'};
        SIprefixes = containers.Map({'f','p','n',[char(181)],'m','','k','M','G','T'}, ...
                                    [1e-15 1e-12 1e-9 1e-6 1e-3 1e0 1e3 1e6 1e9 1e12]);    
    end

    properties (Hidden, SetAccess=protected)
        upDated = 1;
    end
    
    methods (Abstract)
%         [guess, full] = subsref(obj, param)
    end
    
    methods
        UIFigure = datasheet(obj)
    end
    
    methods
%         function obj = component(inputArg1,inputArg2)
%             %UNTITLED Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
        
        function result = addParameter(obj,varargin)
            %obj.addParameter(param) add parameter to component
            %   param -> componentTableData
            %obj.addParameter(name, value)
            %  name -> parameter name as string
            %  value -> parameter typical value in units with SI base 1
            %obj.addParameter(name, value, type)
            %  type -> "min", "max", or "typ"
            %obj.addParameter(name, [typ max min])
            %obj.addParameter(..., conditions)
            % 
            % result is zero if nothing was added, one if a parameter was
            % added and -1 if the parameter is already present
            
%             addparameter(name, value)
%             addparameter(name, value, conditions)
%             
%             addparameter(name, value, type)
%             addparameter(name, value, type, conditions)
%             addparameter(name, value, conditions)
%             
%             addparameter(name, [value])
%             addparameter(name, [value], conditions)
%             
            typVal = [];
            maxVal = [];
            minVal = [];
            testConditions = {};
            typeStr = [];
%             
            result = false;

            if nargin == 1+1
                param = varargin{1};
            else
                paramName = varargin{1};
                [~, paramValid] = isParamOf(obj, paramName);
                assert(paramValid, [paramName ' is not a valid parameter of component type ' class(obj)]);
                
                %commented out and moved to componentTableData constructor
                %varargin{2} = convertUnitsToDefault(obj, varargin{1}, varargin{2}, '');
                
                if nargin == 4+1
                    testConditions = varargin{4};
                    typeStr = varargin{3};
                elseif nargin == 3+1
                    if any(strcmpi(varargin{3}, {'min', 'max', 'typ'}))
                        typeStr = varargin{3};
                    else
                        testConditions = varargin{3};
                    end
                end
                    
                if length(varargin{2}) == 1
                    if isempty(typeStr)
                        typVal = varargin{2};
                    else
                        if strcmpi(typeStr, 'min')
                            minVal = varargin{2};
                        elseif strcmpi(typeStr, 'max')
                            maxVal = varargin{2};
                        elseif strcmpi(typeStr, 'typ')
                            typVal = varargin{2};
                        else
                            error('Invalid type specified');
                        end
                    end
                elseif length(varargin{2}) == 2
                    typVal = varargin{2}(1);
                    maxVal = varargin{2}(2);
                    minVal = [];
                    assert(maxVal >= typVal, 'values must be specified such that min < typ < max');
                elseif length(varargin{2}) == 3
                    typVal = varargin{2}(1);
                    maxVal = varargin{2}(2);
                    minVal = varargin{2}(3);
                    assert(maxVal >= typVal && minVal <= typVal, 'values must be specified such that min < typ < max');
                else
                    error('parameter value must be singleton or a vector of three values');
                end
                
                param = componentTableData(obj, paramName,typVal,maxVal,minVal,testConditions);
                
            end
                
% %                 else
% %                     error('addParameter not defined for this input configuration');
% %                 end
% %                     
% %                 
% %                 
% %                 
% %                 
% %             elseif nargin == 3
% %                 if length(varargin{2}) == 1
% %                     param = componentTableData(obj, varargin{1}, varargin{2}, [], [], []);
% %                 elseif length(varargin{2}) == 3
% %                     assert(varargin{2}(1) < varargin{2}(2) && varargin{2}(1) > varargin{2}(3), 'values must be specified such that min < typ < max');
% %                     param = componentTableData(obj, varargin{1}, varargin{2}(1), varargin{2}(2), varargin{2}(3), []);
% %                 else
% %                     error('addParameter not defined for this input configuration');
% %                 end
% %             elseif nargin >=4
% %                 if length(varargin{2}) == 1
% %                     varargin{2} = convertUnitsToDefault(obj, varargin{1}, varargin{2}, '');
% %                     if strcmpi(varargin{3} == "typ")
% %                         param = componentTableData(obj, varargin{1}, varargin{2}, [], [], []);
% %                     elseif strcmpi(varargin{3} == "min")
% %                         param = componentTableData(obj, varargin{1}, [], [], varargin{2}, [], []);
% %                     elseif strcmpi(varargin{3} == "max")
% %                         param = componentTableData(obj, varargin{1}, [], varargin{2}, [], []);
% %                     end
% %                 elseif  length(varargin{2}) == 3
% %                     assert(varargin{2}(1) < varargin{2}(2) && varargin{2}(1) > varargin{2}(3), 'values must be specified such that min < typ < max');
% %                     varargin{2} = convertUnitsToDefault(obj, varargin{1}, varargin{2}, varargin{3});
% %                     param = componentTableData(obj, varargin{1}, varargin{2}(1), varargin{2}(2), varargin{2}(3), []);
% %             elseif nargin == 5
% %                 param = componentTableData(obj, varargin{1}, varargin{2}, [], [], []);
% %             elseif nargin == 6
% %                 param = componentTableData(obj, varargin{1}, varargin{2}, [], [], []);
% %             elseif nargin == 7
% %                 param = componentTableData(obj, varargin{1}, varargin{2}, [], [], []);
% % %                 componentTableData(type, paramName,typVal,maxVal,minVal,testConditions)
% %             end
            
            %assert(length(param)==1, 'addParameter() can only be used with a single parameter at a time');
            if length(param) > 1
                for i = 1:length(param)
                    result(i) = obj.addParameter(param(i));
                end
            end
            if isa(param, 'componentTableData')
                if isempty(obj.parameters)
                    obj.parameters = param;
                    obj.upDated = 1;
                    result = true;
                else
                    [sameParam, sameData] = eq(obj.parameters,param);
                    sameParam = find(sameParam,1);
                    if ~any(sameParam)
                        % New param
                        obj.parameters(length(obj.parameters)+1) = param;
                        obj.upDated = 1;
                        result = true;
                    elseif ~any(sameData)
                        % Existing parameter, but new data
                        if isempty(obj.parameters(sameParam).min)
                            obj.parameters(sameParam).min = param.min;
                        elseif ~isempty(param.min) && obj.parameters(sameParam).min ~= param.min
                            if nargout == 0
                                warning(['Conflict with exsiting parameter ', param.name, '. Use replaceParam instead to overwrite.  Capture function output to suppress this warning']);
                            end
                            result = false;
                            % obj.parameters(length(obj.parameters)+1) = param;
                        end
                        if isempty(obj.parameters(sameParam).typ)
                            obj.parameters(sameParam).typ = param.typ;
                        elseif ~isempty(param.typ) && obj.parameters(sameParam).typ ~= param.typ
                            if nargout == 0
                                warning(['Conflict with exsiting parameter ', param.name, '. Use replaceParam instead to overwrite.  Capture function output to suppress this warning']);
                            end
                            result = false;
                            % obj.parameters(length(obj.parameters)+1) = param;
                        end
                        if isempty(obj.parameters(sameParam).max)
                            obj.parameters(sameParam).max = param.max;
                        elseif ~isempty(param.max) && obj.parameters(sameParam).max ~= param.max
                            if nargout == 0
                                warning(['Conflict with exsiting parameter ', param.name, '. Use replaceParam instead to overwrite.  Capture function output to suppress this warning']);
                            end
                            result = false;
                            % obj.parameters(length(obj.parameters)+1) = param;
                        end
                        obj.upDated = 1;
                    else
                        % parameter already present, do nothing
                        result = -1;
                    end
%                     obj.parameters(length(obj.parameters)+1) = param; 
                end
            else
                error([class(obj) '.addParameter() not defined for inputs of type ' class(param) ]);
            end
        end
     
        function addGraph(obj, graph)
%             assert(length(graph)==1, 'addGraph() can only be used with a single graph at a time');
             if length(graph) > 1
                for i = 1:length(graph)
                    obj.addGraph(graph(i));
                end
            end
            if isa(graph, 'componentPlotData')
                if isempty(obj.graphs)
                    obj.graphs = graph;
                    obj.upDated = 1;
                else
                    [sameData, samePlots] = eq(obj.graphs,graph);
                    if ~any(samePlots)
                        % New plot
                        obj.graphs(length(obj.graphs)+1) = graph; 
                        obj.upDated = 1;
                    elseif ~any(sameData)
                        % Existing plot, but new data
                        % obj.graphs(length(obj.graphs)+1) = graph; 
                        % warning(['Adding a duplicate plot for ', graph.title]);
                        
                        obj.graphs(find(samePlots,1)).merge(graph);
                        obj.upDated = 1;
                        
                    else
                        % Plot already present, do nothing
                    end
                end
            else
                error([class(obj) '.addParameter() not defined for inputs of type ' class(graph) ]);
            end
        end

        function replaceGraph(obj, ind, graph)
            assert(length(obj.graphs) >= ind, 'Cannot replace nonexistant graph')
            obj.graphs(ind) = graph;
        end

         function replaceParam(obj, newParam)
            params = obj.parameters;
            names = {params.name};
            paramLoc = strcmp(newParam.name,names);
            assert(~isempty(paramLoc), 'Cannot replace nonexistant param')
            obj.parameters(paramLoc) = newParam;
         end

        function deleteGraphTrace(obj,iGraph, iTrace)
            newGraph = obj.graphs(iGraph).deleteTrace(iTrace);
            obj.replaceGraph(iGraph, newGraph);
        end
        
        function merge(obj, newComponent)
           assert( strcmp(obj.partNumber, newComponent.partNumber), ...
               'merge can only be used on two components with the same part number');
           assert(isa(newComponent, class(obj)), ...
                ['class ' class(obj) ' can only be merged with other objects of type ' class(obj)] );
           
           for i = 1:length(newComponent.parameters)
               obj.addParameter(newComponent.parameters(i))
           end
           
           for i = 1:length(newComponent.graphs)
               obj.addGraph(newComponent.graphs(i))
           end
           
%            obj.parameters = [obj.parameters, newComponent.parameters];
%            obj.graphs = [obj.graphs, newComponent.graphs];
           
           if isempty(obj.manufacturer)
               obj.manufacturer = newComponent.manufacturer;
           end
           if isempty(obj.type)
               obj.type = newComponent.type;
           end
           if isempty(obj.material)
               obj.material = newComponent.material;
           end
           
        end
        
        function convValues = convertUnitsToDefault(obj, param, values, units)
            if isempty(units) 
                mult = 1;
            elseif any(strcmp(units(1), obj.SIprefixes.keys))
                mult = obj.SIprefixes(units(1));
            elseif any(strcmp(units, obj.signalUnits)) || isempty(units)
                mult = 1;
            else
                error('Invalid unit specified');
            end
            defaultMult = obj.SIprefixes(obj.defaultMultipliers{strcmp(isParamOf(obj,param), obj.knownParams)});
            convValues = values*mult/defaultMult;
        end
    end
    
    methods (Hidden)
        function varargout = subsref(obj, s)
        %subsref overloads dot-indexing to give back parameters when
        %they are available.  Returns two paramters:
        %   res is a single value best-guess for the response
        %   full is a struct or array of the full relevant data
        % The value in res is determined based on the request and the
        % available data.
% %             [varargout{1:nargout}] = builtin('subsref',obj,s);
%             if (length(param) == 1 && length(obj) == 1) || ...
%                     (length(param) == 1 && length(obj) > 1 && strcmp(param(1).type, '.'))
%                 %% THIS IS STILL WIP -- can I do "tdb(1:3).ron"?
%                 %% added a bunch of obj(1) to try... still need approx() to handle multiple -- i.e. need component edits
%                 [pName, isParam] = isParamOf(obj, param.subs);
%                 if ~isprop(obj(1), param(1).subs) && ~ismethod(obj(1), param(1).subs) && isParam
%                     [~,IA,~] = intersect({obj(1).parameters.name}, pName);
%                     if ~isempty(IA)
%                         objParams = [obj.parameters];
%                         varargout{1} = approx(objParams(IA));
%                         varargout{2} = objParams(IA);
%                     else
%                         [varargout{1:nargout}] = builtin('subsref',obj,param);
%     %                    error(['no such parameter "' param(1).subs '" available for ' class(obj) ' object']);
%                     end
%                 else 
%                     [varargout{1:length(obj)}] = builtin('subsref',obj,param);
%                 end
%             
%             else
%                  [varargout{1:nargout}] = builtin('subsref',obj,param);
%             end 

            if strcmp(s(1).type, '.')
                try 
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                catch e
                    if strcmp(e.identifier,'MATLAB:noSuchMethodOrField')
                       [name,valid] = isParamOf(obj, s(1).subs);
                       if valid
                           loc = strcmpi({obj.parameters.name}, name);
                           if any(loc)
                               if length(s) == 1
                                    varargout{1} = obj.parameters(loc);
                               else
                                   [varargout{1}] = builtin('subsref',obj.parameters(loc),s(2:end));
                               end
                           else

                               blankParam = componentTableData(obj, name, [], [], [], [], '');
                               if length(s) == 1
                                    varargout{1} = blankParam;
                               else
                                   [varargout{1}] = builtin('subsref',blankParam,s(2:end));
                               end
                           end
                       else
                           rethrow(e)
                       end
                           
                    
                    else
                        rethrow(e)
                    end

                end
            else
                 [varargout{1:nargout}] = builtin('subsref',obj,s);
            end
                
        end
        
%          function n = numArgumentsFromSubscript(obj,s,indexingContext)
%             %% redefined to aid in subsref
% %             if strcmp(s(1).type, '()') 
% %                 if strcmp(s(1).subs{:}, ':')
% %                     n = length(obj.components(s(1).subs{:}));
% %                 else
% %                     if length(s) == 1
% %                         n = builtin('numArgumentsFromSubscript',obj.components,s,indexingContext);
% %                     else
% %                         n = 0;
% %                         comps = obj.components(s(1).subs{:});
% %                         for i = 1:length(comps)
% %                             n = n + builtin('numArgumentsFromSubscript',obj.components(i),s(2:end),indexingContext);
% %                         end
% %                     end
% %                 end
% %             else
%                 x=1;
%                 if x==1
%                     n = builtin('numArgumentsFromSubscript',obj,s,indexingContext);
%                 else
%                     n=2;
%                 end
% %             end
%         end
        
        function [name, valid] = isParamOf(obj, name)
            %isParamOf(obj, name) checks if the parameter name is
            %valid for the component type, or if it is an alias of a known
            %parameter. 
            
            if isempty(obj)
                refobj = eval(class(obj));
            else
                refobj = obj(1);
            end

        	%drop any underscores or parentheses
            name = erase(name, {'_', '(', '[', ')', ']'});
            
            %Uniform case
            name = [upper(name(1)), lower(name(2:end))];
            
            %Check if it is a known parameter
            if any(ismember(refobj.knownParams, name))
                valid = true;
                return;
            elseif any(strcmpi(refobj.knownParams, name))
                valid = true;
                return;
            end
            
            % Else, check if it is an aliased name
            [~,IA,~] = intersect(refobj.paramDict.keys, name);
            if ~isempty(IA)
                name = refobj.paramDict(name);
                name = refobj.knownParams{strcmpi(name, refobj.paramNames)};
                valid = true;
            else
                valid = false;
            end
        end
        
    end
    
    methods (Static, Hidden)
        getknownParams % never used at runtime, but makes it easier to type out known parameters for subclasses
    end
end

