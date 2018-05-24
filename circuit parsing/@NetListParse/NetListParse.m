classdef NetListParse < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Symbolic Matrix of ABCD 
        Asym
        Bsym
        Csym
        Dsym
        
        % Htemp of ABCD
        HtempAB
        HtempCD
        
        % Measurement Voltage and Current Nodes
        Meas_Voltage
        Meas_Current
        
        % Find Switch and Diode Position
        Diodes
        Switches
        
        % Names of State, Input, and Output Variables
        StateNames
        OutputNames
        InputNames
        DependentNames
        
        % Netlist files
        NL
        NLnets
        NLwhole
        
        % Filename
        filename
        
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
    
    methods
        
        
        
        
    end
end

