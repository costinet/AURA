classdef NetListParse < handle
    %NETLISTPARSE Holds all of the variables need to parse a netlist
    %   Currently unit checks on initialize function
    
    properties
        % Symbolic Matrix of ABCD
        Asym
        Bsym
        Csym
        Dsym

        % Numerical Matrix of ABCD
        Anum
        Bnum
        Cnum
        Dnum
        Inum % Only used for PLECS

        eigA % The eiganvalues of Anum
        
        
        % Htemp of ABCD (used when unable to solve rref(sym))
        HtempAB
        HtempCD
        dependsAB
        savedCD

        % Find Switch and Diode Position
        Diodes
        Switches
        FETs
        DMpos % In A matrix

        % Names of State, Input, and Output Variables
        StateNames % All state variables [OutputNames;DependentNames]
        OutputNamesCD % Names of Measurements
        InputNames % Such as Vg
        DependentNames % Dependent states in circuit
        OutputNames % Independent states in circuit
        SwitchNames % Names of Switching elements (FETs and Diodes)
        StateNumbers % Number of either voltage or current measurement in Y output vector that corresponds to the state variable index
        StateNumbers_Opposite % Number of either voltage and current source but is the opposite source as StateNumbers
        ConstantNames % Constant Inputs to the system such as Vg or the forward voltage of a diode
        
        %%% Need to change OutputNamesCD to OutputNames
        %%% Need to change OutputNames to MeasurementNames
        %%% For clarity %%%

        % Netlist files
        NL % Numerical Representation of net list
        NLnets % Cell Representation of net list
        NLwhole % Net list file as it is received

        SortedTree % Elements in Tree
        SortedCoTree % Elements in CoTree

        % K value for inductors
        K

        % List the chars of the swithces that are on and off for the
        % state number (column of cell array)
        ON_States
        OFF_States


        % Lists whether a FET or diode is on or off during a state
        ONorOFF

        % List the refernce postion of state variables to and there
        % row location in ABCD matricies relative to their postion in
        % the netlist
        
        OrderedNamesnum
        
        BD_state % Lists the off diode states index in NL

        BD_OFF_state % Lists the off diode states index in NL
        
        
        % The Codex
        Codex
        Component_Values
        
        index
        
        % Forward votlage values in order of A matrix
        Fwd_Voltage
        
    end
    properties (Dependent = true)
    %    Component_Values
        
        
    end
    % These are values used for the numerical parser:
    properties (Access = private)
        Cutset
        SortedTree_cutloop
        SortedCoTree_cutloop
        NewNL
        NewNLnets
        end
        
        % Htemp for AB and CD matrix for large converters

        % Store indexing values hear as well

        % Diode placement from ABCD matrix

        % Identify topology and converter from name of netlist file

        % Create void fuction that when called takes Htemp and creates
        % numerical ABCD matrix

        % Have a public class and function to input measured states in cell
        % format maybe place in converter class???

        % etc...????
   



    properties (Access = private)

        % Measurement Voltage and Current Nodes
        Meas_Voltage
        Meas_Current
        filename
        
    end
    methods
        
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
        end

        function [StateNames]=getStateNames(obj)
            StateNames = obj.StateNames;
        end
        
         function obj = LTSpicecircuitParser(topology)
            assert(isa(topology, 'SMPStopology'), ...
                'input argument topology must be a handle to an object of class SMPStopology');
            obj.topology = topology;
         end

        
         %% Getters
         
        % function res = get.Component_Values(obj)
        %     res = obj.topology.Element_Properties;
        % end
         
         
         
         
    end
end
