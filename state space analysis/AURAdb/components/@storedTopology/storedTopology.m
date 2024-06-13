classdef storedTopology < handle
    %storedTopology minimal implementation of circuit parser for storing
    %topology description
    %   Detailed explanation goes here
    
    properties
        name
        description
    end

    properties (Hidden)
        components
        props
        switches
        schemPositions
        schemWires

        %Only needed for modular Parsers
        baseParser
        moduleParser
        sheetConnectorOutSym 
        sheetConnectorInSym
        sheetConnectorOutFormat
        sheetConnectorInFormat
    end

    properties (Hidden, SetAccess=protected)
        upDated
    end
    

    
    methods
        function obj = storedTopology(name, description, parser)
            %storedTopology() Construct an instance of the storedTopology class
            %   
            %   obj = storedTopology(name, description, varargin)
            arguments
                name {mustBeText} = ''
                description {mustBeText} = ''
                parser circuitParser = NetlistCircuitParser();
            end

            obj.name = name;
            obj.description = description;

            

            if ~isempty(parser)
                if isa(parser,"modularCircuitParser")
                    obj.baseParser = storedTopology(name, 'base' , parser.baseParser);
                    obj.moduleParser = storedTopology(name, 'module' , parser.moduleParser);
                    obj.sheetConnectorOutSym = parser.sheetConnectorOutSym;
                    obj.sheetConnectorInSym = parser.sheetConnectorInSym;
                    obj.sheetConnectorOutFormat = parser.sheetConnectorOutFormat;
                    obj.sheetConnectorInFormat = parser.sheetConnectorInFormat;
                else
                    obj.components = parser.origComponents;
                    obj.props = parser.netListDirectives;
                    if ~isempty(parser.topology)
                        obj.switches = parser.topology.switchLabels;
                    end
                    obj.schemPositions = parser.schemPositions;
                    obj.schemWires = parser.schemWires;
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

