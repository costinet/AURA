classdef SMPSconverter < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Topology

        order % Need right now becuase it breaks code if not here
        Element_Properties % 1st column is char of all element names 2nd column is the value of that element
        Switch_Resistors  % List of chars that represent the switch names plus '_R'
        Switch_Resistor_Values % [SW_OFF; SW_ON; SW_ON] third one will eventually be diode resistance
        Switch_Sequence % Set of binary on or off values that has the number of colums of swithces and the numer of rows of time intervals
        Fwd_Voltage
       
        
        % These are for the active states only 
        ts % time interval length in seconds
        u % input variable (must me in certain order)
        As
        Bs
        Cs
        Ds
        eigA
        
        Simulator % The simulator Class for the Covnerter
        
    end
    
    methods
        
         function Parse_circuit(obj,filename)
        
             % This will delete and set a new Topology class for the
             % converter with the topology given iin the filename
             
             
             top = SMPStopology();
             obj.Topology = top;
             
             parse = NetListParse();
             obj.Topology.Parser = parse;
             
             obj.Topology.Parser.initialize(filename);
             
             
             
             parse = NetListParse();
             parse.initialize(filename);
             parse.ABCD();
             
             top = SMPStopology();
             top.Parser = parse;
             
             conv = SMPSconverter();
             conv.Topology = top;
             
             
         end
         
         function set.Element_Properties(obj,Element_Properties)
             
            if sum(sum(cellfun(@isempty,Element_Properties)))==0
             obj.Element_Properties = Element_Properties;
            else
                error('Not all Element Properties have been declared')
            end
            
         end
         
         function set.ts(obj,ts)
            if sum(ts>0)==length(ts)
                obj.ts = ts;
                
            else
                error('There is a non-postitive time interval length trying to be assigned to the simulation class variable ts')
            end
         
         end
    
    
    end
    
    
end

