classdef circuitParser < handle
    %CIRCUITPARSER abstract class defining requirements of a circuitParser
    %   
    %   See also @LTspiceCircuitParser, @PLECScircuitParser
    
    properties (Abstract, SetAccess = protected)
        sourcefn
        sourcefdate
        
        topology 
    end
    
    methods (Abstract)
        [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrices(obj,circuitPath)
        loadModel(obj, fn, swseq, force)
        updateComponentValues(obj)
        storedTopology = saveTopology(obj,name)
        loadTopology(obj,storedTopolpogy)
    end

    methods (Abstract, Hidden)
        [Iname, Vname] = getSwitchMeasSourceNames(obj,sName)
    end

    methods(Hidden)
        function linkTopology(obj,top)
            obj.topology = top;
        end
    end
    
%     methods
%         %%constructor
%         function obj = circuitParser(topology)
%             assert(isa(topology, 'SMPStopology'), ...
%                 'input argument topology must be a handle to an object of class SMPStopology');
%             obj.topology = topology;
%         end
%     end
end

