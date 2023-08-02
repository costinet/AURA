classdef inductor < component
    %INDUCTOR is the class of inductors for AURA database
    %   INDUCTOR class contains the part number, manufacturer, material
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
        knownParams = {'L','F','I','Rdc','Srf','Isat','Irms','Fspice','R1spice','R2spice','Cspice','K1spice','K2spice','K3spice','K4spice','K5spice','Length','Width','Height'}
        defaultUnits = {'H','Hz','A','Ohm','Hz','A','A','Hz','Ohm','Ohm','F','One','One','One','One','One','m','m','m'}
        defaultMultipliers = {'µ','M','','m','M','','','M','','','p','','','','','','m','m','m'}
        paramDict = containers.Map({'Ls','Lp','fs','f','DCR','SRF','FSPICE','SPICEF','R1SPICE','SPICER1','R2SPICE','SPICER2','CSPICE','SPICEC','K1SPICE','SPICEK1','K2SPICE','SPICEK2','K3SPICE','SPICEK3','K4SPICE','SPICEK4','K5SPICE','SPICEK5','w','h'},{'Inductance','Inductance','Frequency','Frequency','DC Resistance','Self Resonant Frequency','SPICE Frequency Limit ','SPICE Frequency Limit ','SPICE R1','SPICE R1','SPICE R2','SPICE R2','SPICE C ','SPICE C ','SPICE k1','SPICE k1','SPICE k2','SPICE k2','SPICE k3','SPICE k3','SPICE k4','SPICE k4','SPICE k5','SPICE k5','Width','Height'})
        paramNames = {'Inductance','Frequency','Current','DC Resistance','Self Resonant Frequency','Saturation Current','RMS Current','SPICE Frequency Limit ','SPICE R1','SPICE R2','SPICE C ','SPICE k1','SPICE k2','SPICE k3','SPICE k4','SPICE k5','Length','Width','Height'}
        
        
        knownTypes = {'Power', 'Signal'}
        knownMaterials = {'Composite', 'Ferrite', 'Powdered Iron', 'Air'}
    end
    
    methods
        function obj = inductor(partNumber, type, material, varargin)
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

