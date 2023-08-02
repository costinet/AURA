classdef capacitor < component
    %CAPACITOR is the class of inductors for AURA database
    %   CAPACITOR class contains the part number, manufacturer, material
    %   and type information as well as the known parameters for inductors
    
    properties
        partNumber
        manufacturer
        material
        type
    end
    
    %     properties (SetAccess=protected)
    %         parameters = componentTableData.empty;
    %         graphs = componentPlotData.empty;
    %     end

    properties (Hidden, Transient, Constant)
        knownParams = {'C','Vdc','Vac','F','Size','Esr','Z','Temp','Tol','Length','Width','Height'}
        defaultUnits = {'F','V','V','Hz','One','Ohm','Ohm','°C','One','m','m','m'}
        defaultMultipliers = {'µ','','','M','','m','','','','m','m','m'}
        paramDict = containers.Map({'Cap','V','fs','f','Cr','Imp','Tr','l','w','h'},{'Capacitance','DC Rated Voltage','Frequency','Frequency','Equivalent Series Resistance','Impedance','Temperature ','Length','Width','Height'})
        paramNames = {'Capacitance','DC Rated Voltage','AC Rated Voltage','Frequency','Size Code','Equivalent Series Resistance','Impedance','Temperature ','Capacitance Tolerance +/-','Length','Width','Height'}



        knownTypes = {'Power', 'Signal'}
        knownMaterials = {'Ceramic', 'Flim', 'Electrolytic'}
    end

    methods
        function obj = capacitor(partNumber, type, material, varargin)
            %INDUCTOR constructs an instance of the transistor class
            %
            %   inductor(partNumber, type, material, varargin)
            %
            %   partNumber is a string that contains the part number
            %
            %   type is the type of the transistor
            %   Currently known types 
            %   {'Power', 'Signal'}
            %
            %   material is the material of the transistor
            %   Currently known materials 
            %   {'Composite', 'Ferrite', 'Powdered Iron', 'Air'}
            % 
            %   varargin currently under work

            
            if nargin == 0
                return
            else
                obj.partNumber = partNumber;
            end
            
            if nargin > 2
                assert(ismember(type,obj.knownTypes), [type ' is not a known type of ' class(obj) '. Valid types are: ' strjoin(obj.knownTypes, ', ')]);
                obj.type = type;
            end
            
            if nargin > 3
                assert(ismember(material,obj.knownMaterials), [material ' is not a known material of ' class(obj) '. Valid types are: ' strjoin(obj.knownMaterials, ', ')]);
                obj.material = material;
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
                                param = componentTableData(inductor,name,varargin{i+1},[],[],[]);
                                obj.addParameter(param);
                            elseif length(varargin{i+1}) == 3 && isnumeric(varargin{i+1})
                                % three values == [min typ max]
                                param = componentTableData(inductor,name,varargin{i+1}(2),varargin{i+1}(3),varargin{i+1}(1),[]);
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
    end
    
    methods (Hidden)
        function clearUpdated(obj)
            obj.upDated = 0;
        end
    end
    
    
end

