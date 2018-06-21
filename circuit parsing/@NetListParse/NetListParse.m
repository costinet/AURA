classdef NetListParse < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
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
        
        % Htemp of ABCD (used when unable to solve rref(sym))
        HtempAB
        HtempCD
        dependsAB
        savedCD
        
        % Find Switch and Diode Position
        Diodes
        Switches
        
        % Names of State, Input, and Output Variables
        StateNames % All state variables [OutputNames;DependentNames]
        OutputNamesCD % Names of Measurements
        InputNames % Such as Vg
        DependentNames % Dependent states in circuit
        OutputNames % Independent states in circuit
        
        % Netlist files
        NL % Numerical Representation of net list
        NLnets % Cell Representation of net list
        NLwhole % Net list file as it is received
        
        SortedTree % Elements in Tree
        SortedCoTree % Elements in CoTree
        
        % K value for inductors
        K
        
        % Htemp for AB and CD matrix for large converters
        
        % Store indexing values hear as well
        
        % Diode placement from ABCD matrix
        
        % Identify topology and converter from name of netlist file
        
        % Create void fuction that when called takes Htemp and creates
        % numerical ABCD matrix
        
        % Have a public class and function to input measured states in cell
        % format maybe place in converter class???
        
        % etc...????
    end
    
    
    
    properties (Access = private)
        
        % Measurement Voltage and Current Nodes
        Meas_Voltage
        Meas_Current
        filename
        
    end
    methods
        
        function initialize(obj,Filename,Voltage,Current)
            % This function checks to see if 
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
            
        
    end
end

