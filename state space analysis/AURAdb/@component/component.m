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
    
    methods (Abstract)
%         [guess, full] = subsref(obj, param)
    end
    
    methods
        datasheet(obj)
    end
    
    methods
%         function obj = component(inputArg1,inputArg2)
%             %UNTITLED Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
        
        function addParameter(obj,param)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            assert(length(param)==1, 'addParameter() can only be used with a single parameter at a time');
            if isa(param, 'componentTableData')
                if isempty(obj.parameters)
                    obj.parameters = param;
                else
                    [sameData, sameParam] = eq(obj.parameters,graph);
                    if ~any(sameParam)
                        % New plot
                        obj.parameters(length(obj.parameters)+1) = param; 
                    elseif ~any(sameData)
                        % Existing plot, but new data
                        if isempty(obj.parameters(sameParam).min)
                            obj.parameters(sameParam).min = param.min;
                        elseif ~isempty(param.min) && obj.parameters(sameParam).min ~= param.min
                            warning(['Adding a duplicate param for ', param.name]);
                            obj.parameters(length(obj.parameters)+1) = param;
                        end
                        if isempty(obj.parameters(sameParam).typ)
                            obj.parameters(sameParam).typ = param.typ;
                        elseif ~isempty(param.typ) && obj.parameters(sameParam).typ ~= param.typ
                            warning(['Adding a duplicate param for ', param.name]);
                            obj.parameters(length(obj.parameters)+1) = param;
                        end
                        if isempty(obj.parameters(sameParam).max)
                            obj.parameters(sameParam).max = param.max;
                        elseif ~isempty(param.max) && obj.parameters(sameParam).max ~= param.max
                            warning(['Adding a duplicate param for ', param.name]);
                            obj.parameters(length(obj.parameters)+1) = param;
                        end
                    else
                        % Plot already present, do nothing
                    end
%                     obj.parameters(length(obj.parameters)+1) = param; 
                end
            else
                error([class(obj) '.addParameter() not defined for inputs of type ' class(param) ]);
            end
        end
        
        function addGraph(obj, graph)
            assert(length(graph)==1, 'addGraph() can only be used with a single graph at a time');
            if isa(graph, 'componentPlotData')
                if isempty(obj.graphs)
                    obj.graphs = graph;
                else
                    [sameData, samePlots] = eq(obj.graphs,graph);
                    if ~any(samePlots)
                        % New plot
                        obj.graphs(length(obj.graphs)+1) = graph; 
                    elseif ~any(sameData)
                        % Existing plot, but new data
                        obj.graphs(length(obj.graphs)+1) = graph; 
                        warning(['Adding a duplicate plot for ', graph.title]);
                    else
                        % Plot already present, do nothing
                    end
                end
            else
                error([class(obj) '.addParameter() not defined for inputs of type ' class(graph) ]);
            end
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
    end
    
    methods (Hidden)
        function varargout = subsref(obj, param)
        %subsref overloads dot-indexing to give back parameters when
        %they are available.  Returns two paramters:
        %   res is a single value best-guess for the response
        %   full is a struct or array of the full relevant data
        % The value in res is determined based on the request and the
        % available data.
            [varargout{1:nargout}] = builtin('subsref',obj,param);
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
        end
        
        function [name, valid] = isParamOf(obj, name)
            %isParamOf(obj, name) checks if the parameter name is
            %valid for the component type, or if it is an alias of a known
            %parameter. 
            
        	%drop any underscores or parentheses
            name = erase(name, {'_', '(', '[', ')', ']'});
            
            %Uniform case
            name = [upper(name(1)), lower(name(2:end))];
            
            %Check if it is a known parameter
            if any(ismember(obj(1).knownParams, name))
                valid = true;
                return;
            end
            
            % Else, check if it is an aliased name
            [~,IA,~] = intersect(obj(1).paramDict.keys, name);
            if ~isempty(IA)
                name = obj(1).paramDict(name);
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

