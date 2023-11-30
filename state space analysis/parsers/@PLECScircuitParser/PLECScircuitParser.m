classdef PLECScircuitParser < circuitParser
    %PLECScircuitParser Interface to PLECS to extract state space
    %descriptions from simulink-embedded PLECS circuit schematics
    %   
    %   PLECScircuitParser implements a @circuitParser class that uses the
    %   built-in command line interface with PLECS (commercial tool) to
    %   parse state space descriptions.
    %
    %   See also @LTspiceCircuitParser, @SMPSim, @circuitParser
    
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

        missingParams = {}
   end
    
   properties (Hidden)
       diodeCurrentThreshold = 1e-1;
       bodyDiodeVoltage = 0.5;
       bodyDiodeVoltageHysteresis = 0.25;
       removeConstraintMeters = 1;
       VFhystScaleFactor = 0.9;
       VfhystMinEquivResistance = 1;
       diodeMaxExpectedCurrent = 10;
   end
    
    methods
        [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrices(obj,circuitPath, forceRefresh)
        loadModel(obj, fn, swseq, force)
        inspect(obj)
        
         function obj = PLECScircuitParser(topology)
            assert(isa(topology, 'SMPStopology'), ...
                'input argument topology must be a handle to an object of class SMPStopology');
            obj.topology = topology;
         end

        function readOnLoadError(obj, msg, fn)
            obj.sourcefn = fn;
            
            [obj.missingParams.type,~] = regexp(msg, "(?<=Error evaluating parameter ')[\w\s-]+(?=')" , 'match', 'tokens');
            [obj.missingParams.component,~] = regexp(msg, ['(?<=' fn '/)[\w]+(?=\:)' ] , 'match', 'tokens');
            [obj.missingParams.val,~] = regexp(msg, "(?<=Unrecognized function or variable ')[\w\s-]+(?=')" , 'match', 'tokens');
                
        end
    end
end

