classdef transistor < component
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
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
        knownParams = {'Vds','Vsd','Vgs','Vth','Ids','Isd','Ig','Idpulse','Ispulse','Igss','Idss','Rds','Rg','Qg','Qgs','Qgd','Qoss','Qrr','Trr','Irr','Coss','Ciss','Crss','Cgd','Cds','Cgs','Eon','Eoff','Eoss','Rthjc','Rthjb','Rthja','Pd','Gfs','Tr','Tdon','Tf','Tdoff','Tj','Tc'}
        defaultUnits = {'V','V','V','V','A','A','A','A','A','A','A','Ohm','Ohm','C','C','C','C','C','s','A','F','F','F','F','F','F','J','J','J','°C/W','°C/W','°C/W','W','A/V','s','s','s','s','°C','°C'}
        defaultMultipliers = {'','','','','','','m','','','µ','µ','m','m','n','n','n','n','n','µ','','p','p','p','p','p','p','µ','µ','µ','','','','','','n','n','n','n','',''}
        paramDict = containers.Map({'Bv','Bvdss','Vdsmax','Vd','Vf','Vgsth','Idspulse','Isdpulse','Ileak','Idsleak','Idleak','Ron','Rdson','Rgint','Qgate','Qm','Qds','Irrm','Irrmax','Irrpk'},{'Drain-Source Voltage','Drain-Source Voltage','Drain-Source Voltage','Source-Drain Voltage','Source-Drain Voltage','Gate Threshold Voltage','Pulsed Drain-Source Current','Pulsed Source-Drain Current','Drain-Source Leakage Current','Drain-Source Leakage Current','Drain-Source Leakage Current','Drain-Source On-State Resistance','Drain-Source On-State Resistance','Internal Gate Resistance','Total Gate Charge','Gate-Drain Charge','Output Charge','Peak Reverse Recovery Current','Peak Reverse Recovery Current','Peak Reverse Recovery Current'})
        paramNames = {'Drain-Source Voltage','Source-Drain Voltage','Gate-Source Voltage','Gate Threshold Voltage','Drain-Source Current','Source-Drain Current','Gate Current','Pulsed Drain-Source Current','Pulsed Source-Drain Current','Gate-Source Leakage Current','Drain-Source Leakage Current','Drain-Source On-State Resistance','Internal Gate Resistance','Total Gate Charge','Gate-Source Charge','Gate-Drain Charge','Output Charge','Reverse Recovery Charge','Reverse Recovery Time','Peak Reverse Recovery Current','Output Capacitance','Input Capacitance','Reverse Transfer Capacitance','Gate-Drain Capacitance','Drain-Source Capacitance','Gate-Source Capacitance','Turn-on Energy','Turn-off Energy','Coss Stored Energy','Junction-Case Thermal Resistance','Junction-Board Thermal Resistance','Junction-Ambient Thermal Resistance','Power Dissipation','Transconductance','Rise Time','Turn-On Delay Time','Fall Time','Turn-Off Delay Time','Junction Temperature','Case Temperature'}'

     
        knownTypes = {'en-nMOS', 'en-pMOS', 'dep-nMOS', 'dep-pMOS', 'en-npMOS', ...
            'HEMT', 'HFET', 'PT-IGBT', 'NPT-IGBT', 'nJFET', 'pJFET', 'NPN', 'PNP'}
        knownMaterials = {'Si', 'GaN', 'SiC', 'GaAs'}
    end
    
    methods
        function obj = transistor(partNumber, type, material, varargin)
            %transistor() Construct an instance of the transistor class
            %   Detailed explanation goes here
            
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
                                param = componentTableData(transistor,name,varargin{i+1},[],[],[]);
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
    end

    methods (Hidden)
        function clearUpdated(obj)
            obj.upDated = 0;
        end
    end
    
    
end

