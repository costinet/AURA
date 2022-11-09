classdef PLECScircuitParser < circuitParser
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
   properties (SetAccess = private)
        sourcefn
        sourcefdate
        topology
   end
   
   properties (Hidden, SetAccess = private)
        devType
        isFET
        isDiode
        Vf
        allSwitchSignals

        constraintsWarning = 1
   end
    
   properties (Hidden)
       diodeCurrentThreshold = 1e-2;
       bodyDiodeVoltage = 0.5;
       bodyDiodeVoltageHysteresis = 0.25;
       removeConstraintMeters = 1;
   end
    
    methods
        [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrices(obj,circuitPath, forceRefresh)
        loadModel(obj, fn, swseq, force)
        
         function obj = PLECScircuitParser(topology)
            assert(isa(topology, 'SMPStopology'), ...
                'input argument topology must be a handle to an object of class SMPStopology');
            obj.topology = topology;
        end
    end
end

