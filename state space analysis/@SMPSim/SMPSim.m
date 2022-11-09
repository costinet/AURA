classdef SMPSim < handle
    %Simulator object for use with AURA
    % Methods include
%         [ Xs] = SS_Soln(obj, Xi, Bi) 
%         [ Xs] = SS_Soln2(obj, Xi, Bi) 
%         [ output_args ] = regulate( obj )
%         [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
%         [ avgXs, avgYs ] = ssAvgs(obj, Xss)
%         plotAllStates(obj, fn, subplots)
%         plotAllOutputs(obj, fn, subplots)
%         
%         [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xi, Si, Vmax, Vmin, progBar)
%         [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
%     
    properties (Access = private)
%         LSQoptions = optimoptions('lsqlin','algorithm','trust-region-reflective','Display','none');
        tryOpt = 0;
        condThreshold = 1e9;
        gmin = 1/100e6;
        
        % speedup varaibles -> solution memory
        % These are out-of-date (if used)
        oldAs
        oldts
        oldIntEAt
        
        % debugging variables
        debugFigNo = 47;
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
    end
    
    properties (Dependent = true)
        ts
        u
    end
    
    properties
        Xs
        converter
    end
    
    properties (Hidden)
        caching = 0
    end
    
    methods (Access = private)
        %% Private Methods from external files
        [fresp, intEAt] = forcedResponse(obj, A, expA, B, u, t, storeResult)
        [ valid, newt, dist ] = validByInterval(obj, si, Xs) % DEPRICATED
        [ x, xdot ] = stateValue_at_t(obj, x0, t, si)% DEPRICATED
        %AdjustDiodeConduction
        
    end
    
    methods
        %% Methods from external files
        [ Xs] = SS_Soln(obj, Xi, Bi) 
        [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
        [ avgXs, avgYs ] = ssAvgs(obj, Xss)
        
        plotAllStates(obj, fn, oSelect, subplots)
        plotAllOutputs(obj, fn, oSelect, subplots)
        
        [ Xs] = perturbedSteadyState(obj, dts) 
        [J, J2, XssF, XssB, X0, dt] = discreteJacobian(obj, order)
        [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
        [JoutStart,JoutEnd] = discreteJacobianConstraint(obj)
        
        [ output_args ] = regulate( obj )% DEPRICATED
        [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xi, Si, Vmax, Vmin, progBar)% DEPRICATED
        [ Xs] = SS_Soln2(obj, Xi, Bi) % DEPRICATED
        
        %% 
        [violateMargin,targetVal] = checkStateValidity(obj, X, u, swind)
        [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = checkDiscreteErr(obj)
        [tLocs,insertAt,adjType] = findRequiredUncontrolledSwitching(obj,violateMarginStart,violateMarginEnd)
        %[margins(before,after,with &w/o hysteresit)] = checkDiscreteErrors
        % checkContinuousError
        %[altered] = updateForUncontrolledSwitching
        %[deltaTs] = switchingErrorGradientDescent
        
        [Xf,ts,swinds] = timeSteppingPeriod(obj, Xs, ts, origSwind )
        [newts,newswinds] = format1DtimingVector(obj,ts,swinds)

        [weightTotalErr] = getWeightedTotalError(obj, errBefore,errAfter)

        %% Debugging (Verbose) helper functions
        describeDiscreteErrors(obj)
        describeInsertedIntervals(obj, allChanges)
        [T] = describeSwitchState(obj)

        
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
        
        function res = get.oldAs(obj)
            if isempty(obj.oldAs)
                res = zeros(size(obj.As));
            else
                res = obj.oldAs;
            end
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
            obj.clearStoredResults();
        end
        
        function set.u(obj,newU)
            obj.converter.u = newU;
        end
        
        function set.ts(obj,newT)
             error('Setting ts is not recommended for class SMPSsim.  Use methods in SMPSconverter');
%             obj.converter.ts = newT;
        end
        
        
        %% Locally-defined methods
%         function settopology(obj, As, Bs, Cs, Ds)
%            obj.As = As;
%            obj.Bs = Bs;
%            obj.Cs = Cs;
%            obj.Ds = Ds;
%            
%            obj.oldAs = zeros(size(As));
%            obj.oldIntEAt = zeros(size(As));
%            
%            obj.Xs = [];
%         end
%         
%         function setmodulation(obj, ts)
%             obj.ts = ts;
%             obj.oldts = zeros(size(ts));
% 
%             obj.Xs = [];
%         end
%         
%         function setinputs(obj, u)
%             obj.u = u;
%             
%             obj.Xs = [];
%         end

        function Xss = steadyState(obj)
            Xss = SS_Soln(obj);
        end

        function clearStoredResults(obj)
            if obj.caching
                obj.oldAs = zeros(size(obj.As));
                obj.oldIntEAt = zeros(size(obj.As));
                obj.oldts = zeros(size(obj.ts));
            end
            obj.Xs = [];     
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
        
        
        %% Test functions        
        function loadTestConverter(obj,dotmatfile)
            try
                load(dotmatfile, 'conv');
                obj.converter = conv;
%                 params = load(matfile);
            catch err
                ME = MException('AURA:IncompleteConverter', ...
                       'Error: test converter file does not contain all requred variables. Required variables are As, Bs, Cs, Ds, ts, and u');
                throw(ME);
            end
            
%             obj.converter = conv;
%             obj.settopology(conv.topology.As, conv.topology.Bs, conv.topology.Cs, conv.topology.Ds);
%             obj.setmodulation(conv.ts);
%             obj.setinputs(conv.u);
            
            obj.Xs = [];
        end
       
    end
    
end

