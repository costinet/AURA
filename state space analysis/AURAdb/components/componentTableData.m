classdef componentTableData
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = ''
        max = []
        typ = []
        min = []
        conditions = []
    end
    
    properties (Dependent)
        unit
        mult
    end
    
    properties (Hidden)
        baseUnit
        multiplier
    end    
    
    properties (Access = protected)
        componentType
    end
    
    properties (Dependent, Hidden)
       type 
    end
    
    methods
        function obj = componentTableData(type, paramName,typVal,maxVal,minVal,testConditions, units)
            %componentTableData Construct an instance of this class
            %   Detailed explanation goes here
            
            %Make sure the type input is a valid component class
            assert( isa(type, 'component'), 'input ''type'' must be a component class or subclass of component');
            obj.componentType = type;
            
            %Check that paramName is a valid paramter of component class
            [name, valid] = isParamOf(obj.componentType, paramName);
            if valid
%                 knownParamIndex = strcmp(type.paramNames, fullName);
%                 obj.name = type.knownParams{knownParamIndex};
            	obj.name = name;
            else
                warning(['Unknown parameter ' paramName ' for components of type ' class(obj.componentType) '.' newline...
                    'Known parameters are ' strjoin(obj.componentType.knownParams, ', ')]);
            end
            
            if nargin < 7
                units = '';
            end
            
            obj.typ = convertUnitsToDefault(type, obj.name, typVal, units);
            obj.max = convertUnitsToDefault(type, obj.name, maxVal, units);
            obj.min = convertUnitsToDefault(type, obj.name, minVal, units);
            obj.conditions = testConditions;
            
            obj.baseUnit = type.defaultUnits(strcmp(type.knownParams, obj.name));
            obj.multiplier = type.defaultMultipliers(strcmp(type.knownParams, obj.name));
        end
        
        function out = approx(obj)
            %approx Summary of this method goes here
            %   Detailed explanation goes here
            unitScale = obj.componentType.SIprefixes(obj.unit{1});
            if ~isempty(obj.typ)
                out = obj.typ*unitScale;
            elseif ~isempty(obj.max)
                out = obj.max*unitScale;
            elseif ~isempty(obj.min)
                out = obj.min*unitScale;
            else
                out = [];
            end 
        end
        
        function unit = get.unit(obj)
           unit = [obj.multiplier, obj.baseUnit];
        end
        
        function cT = get.type(obj)
            cT = class(obj.componentType);
        end

        function m = get.mult(obj)
            m = obj.componentType.SIprefixes(obj.unit{1});
        end
        
        function [tf, paramTF] = eq(obj, param)
            % [tf, paramTF] = eq(obj, param)
            % tf is true if the two table entries are identical
            % paramTF is true if the two are the same plot, but with
            % different data
            if length(obj) > 1 && length(param) == 1
                for i = 1:length(obj)
                    [xtf, xparamTF] = eq(obj(i), param);
                    tf(i) = xtf;
                    paramTF(i) = xparamTF;
                end
                return
            elseif length(obj) == 1 && length(param) > 1
                for i = 1:length(param)
                    [xtf, xparamTF] = eq(obj, param(i));
                    tf(i) = xtf;
                    paramTF(i) = xparamTF;
                end
                return
            elseif length(obj) > 1 && length(param) > 1
                error('eq() is only defined when one object is singleton');
            else
                if ~isa(param, class(obj))
                    tf = false; 
                elseif ~strcmp(obj.name, param.name)
                    tf = false; 
                elseif (~isempty(param.conditions) && isempty(obj.conditions)) ...
                        || ~strcmp(obj.conditions, param.conditions)
                    tf = false; 
                else
                    tf = true;
                end

                if ~tf
                    paramTF = false; return
                else
                    paramTF = true;
                    if isempty(obj.multiplier)
                        curData = [obj.max, obj.typ, obj.min];
                    else
                        if strcmp(class(obj.multiplier), 'cell')
                            curData = [obj.max, obj.typ, obj.min]*obj.componentType.SIprefixes(obj.multiplier{:});
                        else
                            curData = [obj.max, obj.typ, obj.min]*obj.componentType.SIprefixes(obj.multiplier);
                        end
                    end
                    if isempty(obj.multiplier)
                        newData = [param.max, param.typ, param.min];
                    else
                        if strcmp(class(param.multiplier), 'cell')
                            newData = [param.max, param.typ, param.min]*param.componentType.SIprefixes(param.multiplier{:});
                        else
                            newData = [param.max, param.typ, param.min]*param.componentType.SIprefixes(param.multiplier);
                        end
                    end

                    
                    if numel(newData) == numel(curData) && all(newData == curData | (isnan(newData) & isnan(curData)) )
                        paramTF = true;
                    else
                        paramTF = false;
                    end
                                                
                end
            end
        end
    end
end

