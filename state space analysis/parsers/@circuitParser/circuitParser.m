classdef circuitParser < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, SetAccess = private)
        sourcefn
        sourcefdate
        
        topology 
    end
    
    methods (Abstract)
        [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrices(obj,circuitPath)
        loadModel(obj, fn, swseq, force)
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

