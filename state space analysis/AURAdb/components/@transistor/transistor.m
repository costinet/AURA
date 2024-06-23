classdef transistor < component
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        partNumber
        manufacturer
        material
        type
        package
    end
    
%     properties (SetAccess=protected)
%         parameters = componentTableData.empty;
%         graphs = componentPlotData.empty;
%     end
    
    properties (Hidden, Transient, Constant)
%         knownParams = {'Vds','Vsd','Vgs','Vth','Ids','Isd','Ig','Idpulse','Ispulse','Igss','Idss','Rds','Rg','Qg','Qgs','Qgd','Qoss','Qrr','Trr','Irr','Coss','Ciss','Crss','Cgd','Cds','Cgs','Eon','Eoff','Eoss','Rthjc','Rthjb','Rthja','Pd','Gfs','Tr','Tdon','Tf','Tdoff','Tj','Tc'}
%         defaultUnits = {'V','V','V','V','A','A','A','A','A','A','A','Ohm','Ohm','C','C','C','C','C','s','A','F','F','F','F','F','F','J','J','J','°C/W','°C/W','°C/W','W','A/V','s','s','s','s','°C','°C'}
%         defaultMultipliers = {'','','','','','','m','','','µ','µ','m','m','n','n','n','n','n','µ','','p','p','p','p','p','p','µ','µ','µ','','','','','','n','n','n','n','',''}
%         paramDict = containers.Map({'Bv','Bvdss','Vdsmax','Vd','Vf','Vgsth','Idspulse','Isdpulse','Ileak','Idsleak','Idleak','Ron','Rdson','Rgint','Qgate','Qm','Qds','Irrm','Irrmax','Irrpk'},{'Drain-Source Voltage','Drain-Source Voltage','Drain-Source Voltage','Source-Drain Voltage','Source-Drain Voltage','Gate Threshold Voltage','Pulsed Drain-Source Current','Pulsed Source-Drain Current','Drain-Source Leakage Current','Drain-Source Leakage Current','Drain-Source Leakage Current','Drain-Source On-State Resistance','Drain-Source On-State Resistance','Internal Gate Resistance','Total Gate Charge','Gate-Drain Charge','Output Charge','Peak Reverse Recovery Current','Peak Reverse Recovery Current','Peak Reverse Recovery Current'})
%         paramNames = {'Drain-Source Voltage','Source-Drain Voltage','Gate-Source Voltage','Gate Threshold Voltage','Drain-Source Current','Source-Drain Current','Gate Current','Pulsed Drain-Source Current','Pulsed Source-Drain Current','Gate-Source Leakage Current','Drain-Source Leakage Current','Drain-Source On-State Resistance','Internal Gate Resistance','Total Gate Charge','Gate-Source Charge','Gate-Drain Charge','Output Charge','Reverse Recovery Charge','Reverse Recovery Time','Peak Reverse Recovery Current','Output Capacitance','Input Capacitance','Reverse Transfer Capacitance','Gate-Drain Capacitance','Drain-Source Capacitance','Gate-Source Capacitance','Turn-on Energy','Turn-off Energy','Coss Stored Energy','Junction-Case Thermal Resistance','Junction-Board Thermal Resistance','Junction-Ambient Thermal Resistance','Power Dissipation','Transconductance','Rise Time','Turn-On Delay Time','Fall Time','Turn-Off Delay Time','Junction Temperature','Case Temperature'}'
       knownParams = {'Vds','Vsd','Vgs','Vth','Ids','Isd','Ig','Idpulse','Ispulse','Igss','Idss','Rds','Rg','Qg','Qgs','Qgd','Qoss','Qrr','Trr','Irr','Coss','Ciss','Crss','Cgd','Cds','Cgs','Eon','Eoff','Eoss','Rthjc','Rthjb','Rthja','Pd','Gfs','Tr','Tdon','Tf','Tdoff','Tj','Tc','X','Y','T'}
        defaultUnits = {'V','V','V','V','A','A','A','A','A','A','A','Ohm','Ohm','C','C','C','C','C','s','A','F','F','F','F','F','F','J','J','J','°C/W','°C/W','°C/W','W','A/V','s','s','s','s','°C','°C','m','m','m'}
        defaultMultipliers = {'','','','','','','m','','','µ','µ','m','m','n','n','n','n','n','µ','','p','p','p','p','p','p','µ','µ','µ','','','','','','n','n','n','n','','','m','m','m'}
        paramDict = containers.Map({'Bv','Bvdss','Vdsmax','Vd','Vf','Vgsth','Idspulse','Isdpulse','Ileak','Idsleak','Idleak','Ron','Rdson','Rgint','Qgate','Qm','Qds','Irrm','Irrmax','Irrpk','Cjo','W','L','Z','H'},{'Drain-Source Voltage','Drain-Source Voltage','Drain-Source Voltage','Source-Drain Voltage','Source-Drain Voltage','Gate Threshold Voltage','Pulsed Drain-Source Current','Pulsed Source-Drain Current','Drain-Source Leakage Current','Drain-Source Leakage Current','Drain-Source Leakage Current','Drain-Source On-State Resistance','Drain-Source On-State Resistance','Internal Gate Resistance','Total Gate Charge','Gate-Drain Charge','Output Charge','Peak Reverse Recovery Current','Peak Reverse Recovery Current','Peak Reverse Recovery Current','Output Capacitance','Width','Length','Thickness','Thickness'})
        paramNames = {'Drain-Source Voltage','Source-Drain Voltage','Gate-Source Voltage','Gate Threshold Voltage','Drain-Source Current','Source-Drain Current','Gate Current','Pulsed Drain-Source Current','Pulsed Source-Drain Current','Gate-Source Leakage Current','Drain-Source Leakage Current','Drain-Source On-State Resistance','Internal Gate Resistance','Total Gate Charge','Gate-Source Charge','Gate-Drain Charge','Output Charge','Reverse Recovery Charge','Reverse Recovery Time','Peak Reverse Recovery Current','Output Capacitance','Input Capacitance','Reverse Transfer Capacitance','Gate-Drain Capacitance','Drain-Source Capacitance','Gate-Source Capacitance','Turn-on Energy','Turn-off Energy','Coss Stored Energy','Junction-Case Thermal Resistance','Junction-Board Thermal Resistance','Junction-Ambient Thermal Resistance','Power Dissipation','Transconductance','Rise Time','Turn-On Delay Time','Fall Time','Turn-Off Delay Time','Junction Temperature','Case Temperature','Width','Length','Thickness'}
        knownTypes = {'en-nMOS', 'en-pMOS', 'dep-nMOS', 'dep-pMOS', 'en-npMOS', ...
            'HEMT', 'HFET', 'PT-IGBT', 'NPT-IGBT', 'nJFET', 'pJFET', 'NPN', 'PNP'}
        knownMaterials = {'Si', 'GaN', 'SiC', 'GaAs'}
    end
    
    methods
        function obj = transistor(partNumber, type, material, varargin)
            %transistor() Construct an instance of the transistor class
            %   
            %   obj = transistor(partNumber, type, material)
            
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

        function Ceq = eqCap(obj, type, Vds)
            %EQCAP Linear-equivalent capacitance to voltage-dependent Coss
            %
            %   [Ceq] = eqCap(obj, type, Vds) finds a scalar
            %   linear-equivalent capacitance Ceq for the nonlinear C-V
            %   relationship contained in a transistors graphs data.
            %       -type is 'Q' for charge equivalent or 'E' for energy
            %       equivalent
            %       -Vds is a maximum voltage considered.  The output is
            %       always an equivalent when charged from 0 to Vds.

            graphs = [obj.graphs(:)];
            CplotLoc = find(strcmp({graphs.yLabel}, 'Capacitance'),1,'last');
            if isempty(CplotLoc)
                CplotLoc = find(startsWith({graphs.yLabel}, 'Coss'),1,'last');
            end
            

            if ~isempty(CplotLoc) 
                if length(graphs(CplotLoc).dataLabels) >1
                    curveLoc = strcmp(graphs(CplotLoc).dataLabels, 'Coss');
                else 
                    curveLoc = 1;
                end
                if ~isempty(curveLoc)

            
                    CossV = graphs(CplotLoc).plotData{curveLoc};
        
                    if nargin == 2 || isempty(Vds) || Vds<0
                        Vds = max(CossV(:,1));
                    end
        
                    Vrange = CossV(:,1) < Vds;
                    if sum(Vrange) < 10 || min(Vrange) > 0 || max(Vrange) < Vds
                        %When operating with sparsely sampeld data over the
                        %Vds range
                        newV = linspace(0,Vds,max(10,sum(Vrange)));
                        CossV = [newV', interp1(CossV(:,1),CossV(:,2),newV,'linear','extrap')'];
                        Vrange = ones(size(newV))==1;
                    end

                    if strcmp(type,'Q')
                        Ceq = 1/Vds*trapz(CossV(Vrange,1), CossV(Vrange,2));
                    elseif strcmp(type,'E')
                        Ceq = 2/Vds^2*trapz(CossV(Vrange,1), CossV(Vrange,1).*CossV(Vrange,2));
                    end
                    return
                end
            end
            warning('Unable to locate Coss-vs-Vds plots in transistor graphs.')
%             Ceq = obj.Coss.approx();
            s.type = '.';
            s.subs = 'Coss';
            tableData = subsref(obj,s);
            Ceq = tableData.approx();
        end

        function [Qg, Qgs1, Qgd, Qgs2, Qsw] = gateCharge(obj, Vgs, Vth)
            %GATECHARGE gate charge values for transistor
            %
            %   [Qg, Qgs1, Qgd, Qgs2, Qsw] = gateCharge(obj, Vgs, Vth)
            %   finds total gate charge at voltage Vgs for the transistor
            %   object obj.  
            %   Additional outputs are partial charges found from the plot
            %   of Vgs-vs-Qg, when available.


            Qg = []; Qgs1 = []; Qgd = []; Qgs2 = []; Qsw = [];
            graphs = [obj.graphs(:)];
            QplotLoc = find(startsWith({graphs.title}, 'Vgs-vs-Qg'),1,'last');
            if isempty(QplotLoc)
                QplotLoc = find(startsWith({graphs.yLabel}, 'Vgs') & ...
                    startsWith({graphs.xLabel}, 'Qg'),1,'last');
            end

            params = obj.parameters;
            if ~exist('Vth','var')
                
                VthLoc = strcmp(params(:).name,'Vth');
    
                if ~any(VthLoc)
                    Vth = 2;
                    warning('No Vth value given or available, guessing Vth=2V')
                else
                    Vth = params(VthLoc).approx();
                end
            end

            if isempty(QplotLoc)
                QgLoc = find(strcmp({params(:).name},'Qg'),1);
                if ~isempty(QgLoc)
                    Qg = params(QgLoc).approx();
                    Qgs1 = [];
                    Qgd = [];
                    Qgs2 = [];
                    Qsw = [];
                else
                    warning('Unable to find any data to base Qg approximation on');
                end
            else
                Qgplot = obj.graphs(QplotLoc);
                if ~exist('Vgs','var')
                    Vgs = max(Qgplot.plotData{1}(:,2));
                end
                totalLoc = find(Qgplot.plotData{1}(:,2) >= Vgs,1);
                if isempty(totalLoc)
                    totalLoc = numel(Qgplot.plotData{1}(:,2));
                    if totalLoc == 0
                        totalLoc = [];
                    end
                end
                QgplotData = Qgplot.plotData{1};
                Qg = QgplotData(totalLoc,1);

                dVdQ = diff(Qgplot.plotData{1},1);
                deldVdQ = diff(dVdQ(:,2)./dVdQ(:,1));
                millerStartLoc = find(deldVdQ == min(deldVdQ),1);
                millerEndLoc = find(deldVdQ == max(deldVdQ),1);
                VthLoc = find(Qgplot.plotData{1}(:,2) >= Vth,1);

                Qs = Qgplot.plotData{1}([VthLoc, millerStartLoc, millerEndLoc, totalLoc],1);

                if numel(Qs) == 4
                    Qgs1 = Qs(2);
                    Qgd = Qs(3)-Qs(2);
                    Qgs2 = Qs(4)-Qs(3);
                    Qsw = Qs(3)-Qs(1);
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

