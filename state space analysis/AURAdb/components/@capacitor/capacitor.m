classdef capacitor < component
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        partNumber
        manufacturer
        material
        type
        tolerance
    end

    properties (Dependent)
        Capacitance
        MaxDCvoltage
    end
    
%     properties (SetAccess=protected)
%         parameters = componentTableData.empty;
%         graphs = componentPlotData.empty;
%     end
    
    properties (Hidden, Transient, Constant)
    knownParams = {'C','Esr','Esl','Z','Vdc','Vac','Vrms','Irip','Srf','X','Y','T','F','Tc'}
    defaultUnits = {'F','Ohm','H','Ohm','V','V','V','A','Hz','m','m','m','Hz','°C'}
    defaultMultipliers = {'µ','m','n','','','','','','M','m','m','m','',''}
    paramDict = containers.Map({'Cnom','Rs','Ls','Vdcmax','Vacmax','Vrmsmax','Fres','W','L','Z','H','freq','Tj','Temp'},{'Capacitance','Equivalent Series Resistance','Equivalent Series Inductance','Maximum DC Voltage','Maximum AC Voltage','Maximum RMS Voltage','Self-Resonant Frequency','Width','Length','Thickness','Thickness','Frequency','Temperature','Temperature'})
    paramNames = {'Capacitance','Equivalent Series Resistance','Equivalent Series Inductance','Impedance','Maximum DC Voltage','Maximum AC Voltage','Maximum RMS Voltage','Maximum Ripple Current','Self-Resonant Frequency','Width','Length','Thickness','Frequency','Temperature'}

        knownTypes = {'film','ceramic','electrolytic','tantalum'}
        knownMaterials = {'C0G','X7R','X5R','NPO','JB','X7S','X5S','X6S','CH','X7T','NP0','X8R', 'X8L', 'X8M'} %Using temp codes in place of materials for the moment
    end
    
    methods
        function obj = capacitor(partNumber, type, dielectric, varargin)
            %capacitor() Construct an instance of the transistor class
            %   
            %   obj = capacitor(partNumber, type, dielectric)
            
            if nargin == 0
                return
            else
                obj.partNumber = partNumber;
            end
            
            if nargin > 1
                assert(ismember(type,obj.knownTypes), [type ' is not a known type of ' class(obj) '. Valid types are: ' strjoin(obj.knownTypes, ', ')]);
                obj.type = type;
            end
            
            if nargin > 2
                try
                    assert(ismember(dielectric,obj.knownMaterials), [dielectric ' is not a known material of ' class(obj) '. Valid types are: ' strjoin(obj.knownMaterials, ', ')]);
                catch
                    %de-elevated to warning as dielectric list is not
                    %currently comprehensive, and non-ceramics aren't
                    %included
                    warning([dielectric ' is not a known material of ' class(obj) '. Valid types are: ' strjoin(obj.knownMaterials, ', ')]);
                end
                obj.material = dielectric;
            end
            
            if ~isempty(varargin)
                if mod(length(varargin),2) ~= 0 
                    error('parameter-value inputs must be specified in pairs');
                end
                
                for i = 1:2:length(varargin)
                    if isprop(obj, varargin{i})
%                         eval(['obj.' varargin{i} ' = ''' varargin{i+1} ''';']);
                        S = struct();
                        S.type = '.';
                        S.subs = varargin{i};
                        [~] = builtin('subsasgn',obj, S, varargin{i+1});
                    else
                        [name, valid] = isParamOf(obj, varargin{i});
                        if valid
                            if length(varargin{i+1}) == 1 && isnumeric(varargin{i+1})
                                % single value --  typical
                                param = componentTableData(capacitor,name,varargin{i+1},[],[],[]);
                                obj.addParameter(param);
                            elseif length(varargin{i+1}) == 3 && isnumeric(varargin{i+1})
                                % three values == [min typ max]
                                param = componentTableData(transistor,name,varargin{i+1}(2),varargin{i+1}(3),varargin{i+1}(1),[]);
                                obj.addParameter(param);
                            else
                                error('not defined 2');
                            end
                        else
                            warning([varargin{i} ' is not a know parameter of class ' class(obj)]);
                        end
                    end
                end
            end
            obj.upDated = 1;
        end

        function C = get.Capacitance(obj)
            loc = strcmpi({obj.parameters.name}, 'C');
            if any(loc)
                C = approx(obj.parameters(loc));
            else
                C = [];
            end
        end

        function V = get.MaxDCvoltage(obj)
            loc = strcmpi({obj.parameters.name}, 'Vdc');
            if any(loc)
                V = approx(obj.parameters(loc));
            else
                V = [];
            end

        end

        function [C, V] = capacitanceAtV(obj,Vdc)
            graphs = [obj.graphs(:)];
            CplotLoc = find(contains({graphs.yLabel}, 'Capacitance') ...
                & contains({graphs.xLabel}, 'Voltage'),1,'last');

            if isempty(CplotLoc)
                CplotLoc = find(startsWith({graphs.yLabel}, 'C, ')...
                    & startsWith({graphs.yLabel}, 'Vdc, '),1,'last');
            end

            C = [];
            V = [];

            try
                if ~isempty(CplotLoc)
                    if length(graphs(CplotLoc).dataLabels) >1
                        curveLoc = strcmp(graphs(CplotLoc).dataLabels, 'C');
                    else 
                        curveLoc = 1;
                    end
                    if ~isempty(curveLoc)
    
                        CV = graphs(CplotLoc).plotData{curveLoc};
            
                        if nargin == 1
                            C = CV(:,2);
                            V = CV(:,1);
                            return
                        end
    
                        C = interp1(CV(:,1),CV(:,2),abs(Vdc),[],'extrap');
                        V = abs(Vdc);
                        if V > obj.MaxDCvoltage
                            warning('Extrapolating data beyond maximum rated voltage');
                        end
                    end
                end
            catch
                %% Quick fix for bad plots
                warning(['Capacitor ' obj.partNumber ' has bad C-V plot']);
                C = [];
                V = [];
            end

        end

         function [ESR, fo] = ESRatf(obj,f)
            graphs = [obj.graphs(:)];
            ZplotLoc = find(~cellfun(@isempty,strfind({graphs.yLabel}, 'Impedance')),1,'last');
            if isempty(ZplotLoc)
                ZplotLoc = find(startsWith({graphs.yLabel}, 'Z, '),1,'last');
            end

            ESR = [];
            fo = [];

            if ~isempty(ZplotLoc)
                if length(graphs(ZplotLoc).dataLabels) >1
                    curveLoc = strcmp(graphs(ZplotLoc).dataLabels, 'Z');
                else 
                    curveLoc = 1;
                end
                if ~isempty(curveLoc)

                    ZF = graphs(ZplotLoc).plotData{curveLoc};
        
                    if nargin == 1
                        ESR = real(ZF(:,2));
                        fo = ZF(:,1);
                        return
                    end

                    Z = interp1(ZF(:,1),ZF(:,2),abs(f),[],'extrap');
                    ESR = real(Z);
                    fo = abs(f);
                    if fo > max(ZF(:,1)) || fo<min(ZF(:,1))
                        warning('Extrapolating data beyond maximum tested frequency');
                    end
                end
            end

         end

         function cClass = ceramicClass(obj)
            cClass = zeros(numel(obj),1);
            isCeramic = strcmpi({obj.type}, 'ceramic');
            isClassII = startsWith({obj.material}, {'X','Y','Z', 'J'});
            isClassI = startsWith({obj.material}, ...
                        {'C', 'B', 'L', 'A', 'M', 'N', 'P', 'R', 'S', 'T', 'V', 'U'});
            if isCeramic
                if any(isClassII & isCeramic)
                    cClass(isClassII & isCeramic) = 2;
                end
                if any(isClassI & isCeramic)
                    cClass(isClassI & isCeramic) = 1;
                end
            end
         end

        
    end
     


    methods (Hidden)
        function clearUpdated(obj)
            obj.upDated = 0;
        end
    end
    
    
end

