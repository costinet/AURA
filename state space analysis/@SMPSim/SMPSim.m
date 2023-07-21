classdef SMPSim < handle
    %Simulator object for use with AURA
    % Methods include
%     
    properties
        Xs
        converter
    end

    properties (Dependent = true)
        As
        Bs
        Cs
        Ds
        Is
      
        topology
        stateNames
        outputNames
        switchNames
        inputNames

        ts
        u
    end
    

    
    properties (Hidden)
        % debugging variables
        debug = 0;
        debug2 = 0;
        debugFigNo = 47;
    
        timeSteppingInit = 0;
        finalRunMethod = 0;

        maxItns = 100;

        IHC = [];

        tryOpt = 0;
        condThreshold = 1e9;        
    end
    
    methods (Access = private)
        %% Private Methods from external files
        [ Xs] = SS_Soln(obj, Xi, Bi) 
        [fresp, intEAt] = forcedResponse(obj, A, expA, B, u, t, storeResult) 
        [ Xs] = AugmentedSteadyState(obj, dts)

        [J, J2, XssF, XssB, X0, dt] = discreteJacobian(obj, order)
        [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
        [JoutStart,JoutEnd] = discreteJacobianConstraint(obj)

        [violateMargin,targetVal] = checkStateValidity(obj, X, u, swind)
        [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = checkDiscreteErr(obj)
        [tLocs,insertAt,adjType] = findRequiredUncontrolledSwitching(obj,violateMarginStart,violateMarginEnd)

        [Xf,ts,swinds] = timeSteppingPeriod(obj, Xs, ts, origSwind )
        [newts,newswinds] = format1DtimingVector(obj,ts,swinds)
        [weightTotalErr] = getWeightedTotalError(obj, errBefore,errAfter)

%         [ valid, newt, dist ] = validByInterval(obj, si, Xs) % DEPRICATED
%         [ x, xdot ] = stateValue_at_t(obj, x0, t, si)% DEPRICATED
%         [ Xs] = perturbedSteadyState(obj, dts) %DEPRICATED
%         [ Xs] = SS_Soln2(obj, Xi, Bi) % DEPRICATED
%         [ output_args ] = regulate( obj )% DEPRICATED
%         [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xi, Si, Vmax, Vmin, progBar)% DEPRICATED
        %[margins(before,after,with &w/o hysteresit)] = checkDiscreteErrors
        % checkContinuousError
        %[altered] = updateForUncontrolledSwitching
        %[deltaTs] = switchingErrorGradientDescent
        
        
    end
    
    methods
        %% Methods from external files
        
        [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
        [ avgXs, avgYs ] = ssAvgs(obj, Xss)
        plotWaveforms(obj, type, fn, oSelect, subplots)     
        niter = findValidSteadyState(obj)

        %% Debugging (Verbose) helper functions
        describeDiscreteErrors(obj)
        describeInsertedIntervals(obj, allChanges)
        [T] = describeSwitchState(obj)
        describeAlteredTimes(obj,oldts)

        
        %% Constructors
        function obj = SMPSim(conv)
            if nargin == 1
                obj.converter = conv;
            elseif nargin == 0
                conv = SMPSconverter();
                top = SMPStopology();
                conv.topology = top;
                obj.converter = conv;
            end
        end
        
        %% Getters
        function res = get.As(obj)
            res = obj.converter.As;
        end
        
        function res = get.Bs(obj)
            res = obj.converter.Bs;
        end
        
        function res = get.Cs(obj)
            res = obj.converter.Cs;
        end
        
        function res = get.Ds(obj)
            res = obj.converter.Ds;
        end
        
        function res = get.Is(obj)
            res = obj.converter.Is;
        end
        
        function res = get.ts(obj)
            res = obj.converter.ts;
        end
        
        function res = get.u(obj)
            res = obj.converter.u;
        end
        
        function res = get.topology(obj)
            res = obj.converter.topology;
        end
        
        
        function res = get.stateNames(obj)
            res = obj.converter.topology.stateLabels;
        end
        
        function res = get.outputNames(obj)
            res = obj.converter.topology.outputLabels;
        end
        
        function res = get.inputNames(obj)
            res = obj.converter.topology.inputLabels;
        end
        
        function res = get.switchNames(obj)
            res = obj.converter.topology.switchLabels;
        end
            
        
        %% Setters
        function set.converter(obj, conv)
            obj.converter = conv;
            obj.Xs = [];
        end
        
        function set.u(obj,newU)
            obj.converter.u = newU;
        end
        
        function set.ts(varargin)
             error('Setting ts is not recommended for class SMPSsim.  Use methods in SMPSconverter');
        end
        
        
        %% Locally-defined methods
        function Xss = steadyState(obj, dts)
%             try
            if nargin > 1 
                Xss = obj.AugmentedSteadyState(dts);
            else
                Xss = obj.AugmentedSteadyState();
            end

%             catch e 
%                 Xss = SS_Soln(obj);
%             end
        end
  

        function plotAllStates(obj, fn, oSelect, subplots)
            if(nargin <= 2)
                subplots = 1;
                oSelect = 1:size(obj.Xs,1); 
            elseif(nargin <=3)
                subplots = 1;
            end
            plotWaveforms(obj, 1, fn, oSelect, subplots);
        end

        function plotAllOutputs(obj, fn, oSelect, subplots)
            if(nargin <= 2)
                subplots = 1;
                oSelect = 1:size(obj.Cs,1); 
            elseif(nargin <=3)
                subplots = 1;
            end
            plotWaveforms(obj, 2, fn, oSelect, subplots);
        end
    

        
        function loc = sigLoc(obj, name, type)
            if nargin == 3
                if strcmp(type, 'x') %state
                    loc = find(strcmp(obj.topology.stateLabels, name));
                elseif strcmp(type, 'u') %input
                    loc = find(strcmp(obj.topology.inputLabels, name));
                elseif strcmp(type, 'y') %output 
                    loc = find(strcmp(obj.topology.outputLabels, name));
                elseif strcmp(type, 'sw') %switch
                    loc = find(strcmp(obj.topology.switchLabels, name));
                else
                   error('Incompatable type.  Variable type must be x, u, y, or s'); 
                end
            elseif nargin == 2
                allLabels = [obj.topology.stateLabels; 
                    obj.topology.inputLabels; 
                    obj.topology.outputLabels; 
                    obj.topology.switchLabels'];
                lengths = [0; cumsum([length(obj.topology.stateLabels); 
                    length(obj.topology.inputLabels); 
                    length(obj.topology.outputLabels)])];
               
                loc = find(strcmp(allLabels, name));
                if ~isempty(loc)
                    loc = loc - max(lengths(lengths<loc));
                end
            else
                loc = [];
            end
            if isempty(loc)
                error('cannot find specified signal in topology');
            end
        end
              
    end
    
end

