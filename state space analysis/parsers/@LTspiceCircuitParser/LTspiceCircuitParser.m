classdef LTspiceCircuitParser < circuitParser
    %LTspiceCircuitParser circuitParser handling spice netlist files
    %   LTspiceCircuitParser interfaces spice netlists (LTspice or
    %   otherwise) to parse state space descriptions of circuit topologies.
    %
    %   Currently only V, M, D, L, C, R, I elements are supported.  K
    %   statements (inductor coupling) have limited support.  Parametric
    %   definitions of variables which are available in the MATLAB base
    %   workspace are supported.
    %
    %   Contributed by J. Baxter
    %
    %   Methods replicated from L. Chua and P-M Lin, *Computer Aided
    %   Analysis of Electronic Circuits: Algorithms and Computational
    %   Techniques* and referenced works
    %
    %   See also @PLECScircuitParser, @SMPSim, @circuitParser

    properties 
        %Installation directory of LTspice.  Replace this with your own.
        LTSpiceFolder = ['C:\Users\' getenv('USERNAME') '\Documents\LTspiceXVII'];
        LTSpiceExe = ['C:\Program Files\LTC\LTspiceXVII\' 'XVIIx64.exe'];
        method = 'new'
    end
    
    properties (SetAccess = protected)
        sourcefn
        sourcefdate
        ascfn
        
        topology 
		origComponents
        netListDirectives

        components
		
		 % Symbolic Matrix of ABCD
% % %         Asym
% % %         Bsym
% % %         Csym
% % %         Dsym

        % Numerical Matrix of ABCD
        Anum
        Bnum
        Cnum
        Dnum
        Inum % Only used for PLECS

        eigA % The eiganvalues of Anum
        
        
       

        % Find Switch and Diode Position
        Diodes
        Switches
        FETs
        DMpos % In A matrix

        % Names of State, Input, and Output Variables
        StateNames % All state variables [OutputNames;DependentNames]
        OutputNamesCD % Names of Measurements
        OutputNames % Independent states in circuit
        ConstantNames % Constant Inputs to the system such as Vg or the forward voltage of a diode

        %%% Need to change OutputNamesCD to OutputNames
        %%% Need to change OutputNames to MeasurementNames
        %%% For clarity %%%

        % Netlist files
        NL % Numerical Representation of net list
        NLnets % Cell Representation of net list

       
        Component_Values = {}
        
        index
        
        % Forward votlage values in order of A matrix
        Fwd_Voltage
    end


	
	 properties (Access = private)
        Cutset
        SortedTree_cutloop
        SortedCoTree_cutloop
%         NewNL
%         NewNLnets
        % Htemp for AB and CD matrix for large converters
        % Store indexing values hear as well
        % Diode placement from ABCD matrix
        % Identify topology and converter from name of netlist file
        % Create void fuction that when called takes Htemp and creates
        % numerical ABCD matrix
        % Have a public class and function to input measured states in cell
        % format maybe place in converter class???
        % etc...????

        % Measurement Voltage and Current Nodes
        Meas_Voltage = {}
        Meas_Current = {}
        filename
        
        

        %% Netlist things
        netlistLibraries
        netlistModels
     end

     properties (Hidden)
        defaultRoff = 10e6;

        undefinedExpressions = {}

        %% Depricated (?) properties
%         % Htemp of ABCD (used when unable to solve rref(sym))
%         HtempAB
%         HtempCD
%         dependsAB
%         savedCD
% 
%         InputNames % Such as Vg
%         DependentNames % Dependent states in circuit
%         SwitchNames % Names of Switching elements (FETs and Diodes)
%         StateNumbers % Number of either voltage or current measurement in Y output vector that corresponds to the state variable index
%         StateNumbers_Opposite % Number of either voltage and current source but is the opposite source as StateNumbers
% 
%          NLwhole % Net list file as it is received
%         SortedTree % Elements in Tree
%         SortedCoTree % Elements in CoTree
% 
%         % K value for inductors
%         K
% 
%         % List the chars of the swithces that are on and off for the
%         % state number (column of cell array)
%         ON_States
%         OFF_States
% 
% 
%         % Lists whether a FET or diode is on or off during a state
        ONorOFF
% 
%         % List the refernce postion of state variables to and there
%         % row location in ABCD matricies relative to their postion in
%         % the netlist
%         OrderedNamesnum
%         
%         BD_state % Lists the off diode states index in NL
%         BD_OFF_state % Lists the off diode states index in NL
%         
%         % The Codex
%         Codex
        
     end

   properties (Access = private)
       % From Jared' Topology Class
        order % Need right now becuase it breaks code if not here
        Element_Properties = {}% 1st column is char of all element names 2nd column is the value of that element ** Resistor values for switches must be at the end
        Switch_Resistors  % List of chars that represent the switch names plus '_R'
        Switch_Resistor_Values % [SW_OFF; SW_ON; SW_ON] third one will eventually be diode resistance
        Switch_Sequence % Set of binary on or off values that has the number of colums of swithces and the numer of rows of time intervals
        Switch_Names % The names of the switch resistors

        devType
        isFET
        isDiode
        Vf

        parsedU
    end
    
    methods
        [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrices(obj,circuitPath)
        loadModel(obj, fn, swseq, force)

        readSpiceNetlist(obj,filename)
        evalSpiceParams(obj,param)
        str = spiceNumFormat(obj,str)
        component = parseSpiceComponent(obj, str, type)
        [paramVal, embeddedParams, paramExpr] = parseTwoNetSpiceComponent(obj,str,component)
        [params] = parseSpiceParamList(obj, str, params)
        [components] = deEmbedSpiceParams(obj,component,embeddedParams)
        
        linearizeCircuitModel(obj)
        linearizeCircuitModel2(obj)
        updateComponentValues(obj)
        components = switchLinearSubcircuit(obj, component)
        components = XFdependentSourceSubcircuit(obj, directive, components)

        setSwitchingState(obj, swvec)
        [A,B,C,D,I] = solveStateSpaceRepresentation(obj)
	
		function set.ONorOFF(obj,value)
            obj.ONorOFF = value;
        end
        
        function [value] = get.ONorOFF(obj)
            value = obj.ONorOFF;
        end
            
            
        function initialize(obj,Filename,Voltage,Current,Num_solve)
            % This function checks to see if all variables are of the
            % correct class
            
            if nargin == 5
                if iscellstr(Num_solve(:,1))
                    try
                        numbers = cell2mat(Num_solve(:,2));
                    catch ME
                        error('Ensure the format of the component values array is of the form: component | numeric value')
                    end
                    if isnumeric(numbers)
                        obj.Component_Values = Num_solve;
                    end
                else
                    error('Ensure the format of the component values array is of the form: component | numeric value')
                end
            end
            
            if ~exist(Filename,'file')
                error('File %s not found',Filename)
            end

            if ischar(Filename)
                obj.filename=Filename;
            else
                error('Filename must be of type char \nCurrent class of filename: %s',class(Filename))
            end
            if nargin == 2
                obj.Meas_Voltage = {};
                obj.Meas_Current = {};
                return
            end
            if nargin == 3
                warning('Only voltage measurements detected')
                Current = {};
            end
            if iscell(Voltage)
                if size(Voltage,1)==1||size(Voltage,2)==1||isempty(Voltage)
                    obj.Meas_Voltage = Voltage;
                else
                    error('Voltage is not the correct dimension')
                end
            else
                error('Voltage must be of type cell \nCurrent class of Voltage: %s',class(Voltage))
            end
            if iscell(Current)
                if size(Current,1)==1||size(Current,2)==1||isempty(Current)
                    obj.Meas_Current = Current;
                else
                    error('Current is not the correct dimension')
                end
            else
                error('Current must be of type cell \nCurrent class of Current: %s',class(Current))

            end

            obj.sourcefn = Filename;
            readSpiceNetlist(obj,Filename)
        end

        function [StateNames]=getStateNames(obj)
            StateNames = obj.StateNames;
        end
        
         function obj = LTspiceCircuitParser(topology)
             if nargin ==1 
                assert(isa(topology, 'SMPStopology'), ...
                    'input argument topology must be a handle to an object of class SMPStopology');
                obj.topology = topology;
             else
                 % warning('Incorrect number of inputs supplied')
             end

             if ~exist(obj.LTSpiceFolder,'dir')
                 warning('Value of @LTSpiceCircuitParse.LTSpiceFolder is not the installation directory of LTSpice.  Libarary use may be limited')
             end
         end

    end

    methods(Hidden)
        function [Iname, Vname] = getSwitchMeasSourceNames(obj,sName)
            if ~isa(sName,'cell')
                sName = {sName};
            end

            for i = 1:numel(sName)
                Iname{i} = ['Im_' sName{i}];
                Vname{i} = ['Vm_' sName{i}];
            end
        end

        function setComponents(obj,compList)
            obj.origComponents = compList;
        end

        function storedTopology = saveTopology(obj,name)
            error('Not currently implemented');
        end
        function loadTopology(obj,storedTopolpogy)
            error('Not currently implemented');
        end

    end
end

