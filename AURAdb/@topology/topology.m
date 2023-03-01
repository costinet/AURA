classdef topology < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        Name
        LTspice % binary 1/0 if a LTspice file 
        PLECS % binary 1/0 if a PLECS file
        NetList % Netlist file of the 
        Schematic % .jpeg or png of 
        
        
      
    end
    
    methods
    
        function obj = topology(Name,LTspice,PLECS,NetList,Schematic)
            
            if nargin == 0
                return
            else
                obj.Name = Name;
            end
            
            if nargin > 1
                obj.LTspice = LTspice;
            end
            
            if nargin > 2
                obj.PLECS = PLECS;
            end
            
            if nargin > 3
                obj.NetList = NetList;
            end
            
            if nargin > 4
                obj.Schematic = Schematic;
            end
            
            
        end
        
        
        function merge(obj, newComponent)
           
            
            for i = 1:length(newComponent.parameters)
                obj.addParameter(newComponent.parameters(i))
            end
            
            for i = 1:length(newComponent.graphs)
                obj.addGraph(newComponent.graphs(i))
            end
            
            %            obj.parameters = [obj.parameters, newComponent.parameters];
            %            obj.graphs = [obj.graphs, newComponent.graphs];
            
            if isempty(obj.manufacturer)
                obj.manufacturer = newComponent.manufacturer;
            end
            if isempty(obj.type)
                obj.type = newComponent.type;
            end
            if isempty(obj.material)
                obj.material = newComponent.material;
            end
            
        end
        
        
    end
    
end

