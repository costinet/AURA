classdef SMPSim < handle
    %Simulator object for use with AURA
    %   SMPSim works in conjunction with SMPSTopology and SMPSConverter to
    %   simulate steady-state performance of Switched Mode Power Supplies.
    %
    %   Constructors:
    %   obj = SMPSim() returns an empty class object
    %   obj = SMPSim(circuitPath) circuitPath is a character vector path to an 
    %       LTspice or PLECS (simulink-embedded) file containing a description of a converter.  
    %   obj = SMPSim(circuitPath, swvec, us, ts) will load the specified
    %   modulation pattern and set the input vector.  Inputs can be
    %   ommitted by replacing them with []
    %
    %   See also SMPSconverter, SMPStopology, AURAdb
   
    properties
        Xs % Current steady-state solution e.g. as solved by SteadyState
        converter % link to a SMPSconverter class object specifying the converter
    end

    properties (Dependent = true)
        As % 3-dimensional matrix of values for Ai, where As(:,:,i) is the 2D square matrix Ai during the ith interval.
        Bs % 3-dimensional matrix of values for Bi, where Bs(:,:,i) is the 2D matrix/vector Bi during the ith interval.
        Cs % 3-dimensional matrix of values for Ci, where Cs(:,:,i) is the 2D square matrix Ci during the ith interval.  
        Ds % 3-dimensional matrix of values for Di, where Ds(:,:,i) is the 2D vector Di during the ith interval.
        Is % 3-dimensional matrix of values for Ii, where Is(:,:,i) is the 2D vector Ii during the ith interval.
      
        topology % link to a SMPStopology class object specifying the converter topology
        parser % link to a circuitParser class interfancing the topology description
        stateNames % Cell array with names of states ordered to correspond to their order in Xs
        outputNames % Cell array with names of outputs ordered to correspond to their order in Y
        switchNames % Cell array with names of switches ordered to correspond to their order in the SMPSconverter swvec
        inputNames  % % Cell array with names of inputs ordered to correspond to their order in u

        ts % ts is a vector of the time durations of each inverval
        u % u is the (assumed constant within each interval) independent input vector
        swvec %swvec is a matrix of binary switch states.  swvec(i,j) is the on/off status of switch i during time interval j

        Ys %Outputs at current steady-state solution
        YsEnd % Outputs at the end of each subinterval for current steady-state solution
    end

    properties (Hidden, Dependent)
        fullu
    end
    

    
    properties (Hidden)
        % debugging variables
        debug = 0;
        debug2 = 0;
        debugFigNo = 47;
    
        suppressIterationOutput = 1; %will display weighted error once per iteration in findValidSteadyState 

        timeSteppingInit = 0;   
        finalRunMethod = 0; %0= eig2tis, 1=timestepping

        maxItns = 100;

        allowHalfCycleReduction = 0;
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

        [Xf,ts,swinds,us] = timeSteppingPeriod(obj, Xs, ts, origSwind, us)
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

    methods (Hidden)
        plotWaveforms(obj, type, fn, oSelect, subplots)    

        %% Debugging (Verbose) helper functions
        describeDiscreteErrors(obj)
        describeInsertedIntervals(obj, allChanges)
        [T] = describeSwitchState(obj)
        describeAlteredTimes(obj,oldts)

        [ avgXs, avgYs ] = ssAvgsInvert(obj, Xss)
    end
    
    methods
        %% Methods from external files
        [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
        [ avgXs, avgYs, Ints ] = ssAvgs(obj, Xss)
        niter = findValidSteadyState(obj)
        Gz = findSSTF(obj, tp, oi, tm, ut)
        Gs = findAvgModelTF(obj, tp, oi)
        [ xs, t, ys, xd ] = simTransient(obj, nPeriods, Xss, tsteps)

        %% Locally-defined methods

         %% Constructors
        function obj = SMPSim(circuitPath, swvec, us, ts)
        %construct a new object of class SMPSim
        % 
        %   obj = SMPSim(circuitPath, swvec, us, ts) all inputs are
        %   optional. 

            conv = SMPSconverter();
            top = SMPStopology();
            conv.topology = top;
            obj.converter = conv;

            if nargin == 0
                return
            elseif nargin == 1
                initialize(obj, circuitPath);
            elseif nargin == 2
                initialize(obj, circuitPath, swvec);
            elseif nargin == 3
                initialize(obj, circuitPath, swvec, us);
            elseif nargin == 4
                initialize(obj, circuitPath, swvec, us, ts);
            end
            
            
        end
        
        function [top, conv] = initialize(obj, circuitPath, swvec, us, ts)
        % load circuit file and necessary parameters for simulation
        %
        %   initialize combines calls to the topology's loadCircuit and the
        %   converter's setSwitchingPattern functions, along with some
        %   format checks.
        %
        %   [top, conv] = initialize(obj, circuitPath, swvec, us, ts)
        %   circuitPath is the path to a parseable circuit definition file
        %   swvec is a 2D binary vector describing the switching states
        %   us is the input vector
        %   ts is the time vector
        %   top and conv are the object's reference topology and converter
        %  
        %   See also SMPSconverter.setSwitchingPattern, SMPStopology.loadCircuit
            if nargin == 2
                swvec = [];
                us = [];
                ts = [];
            elseif nargin == 3
                us = [];
                ts = [];
            elseif nargin == 4
                ts = [];
            end

            obj.IHC = [];           %Need to drop any half-cycle identification 
            %Specific issue if you findValidSteadyState, get IHC, then re-run 
            %initialize, then plot a steadyState -without- another
            %findValidSteadyState, you get multiple inverted periods

            top = obj.topology;
            conv = obj.converter;

            if ~isempty(swvec)
                top.loadCircuit(circuitPath,swvec,1);
            else
                top.loadCircuit(circuitPath,[],1);
                try 
                    swvec = evalin('base','swvec');
                catch e
                    return
                end
            end

            if ~isempty(us)
                obj.u = us;
            else
                try
                    obj.u = evalin('base', 'us');
                catch e
                    return
                end
            end

            if size(obj.u,1) ~= size(obj.Bs,2)
                if(size(obj.u,2) ~= size(obj.Bs,2)) && size(obj.u,1) == 1
                    obj.u = reshape(obj.u, [size(obj.Bs,2), 1, size(obj.u,3)]);
                else
                    error('Invalid size of u specified.')
                end
            end
            
             if ~isempty(ts)
                conv.setSwitchingPattern(swvec, ts);
             else
                try
                    ts = evalin('base','ts');
                catch e
                    return
                end
                try
                    conv.setSwitchingPattern(swvec, ts);
                catch e
                    if nargin < 5
                        % They didn't supply all parameters, and finding
                        % them on base didn't work.
                        warning('Not all parameters defined, SMPSim object will need to be updated with ts, us, and swvec before running.')
                        return
                    else
                        rethrow(e);
                    end
                end
            end
            
        end


        function Xss = steadyState(obj, dts)
        %solve steady-state solution of converter at current switching pattern and timing
        %
        %   [ Xss] = steadyState(obj) finds the state values Xss in
        %   steady-state for the switched system described by the SMPSim object obj
        %
        %   [ Xss] = steadyState(obj,dts) finds the steady-state solution
        %   with perturbations dts to the timing of each interval.
        %   
        %   Note that steadyState does not consider state-dependent
        %   swithing.  For that use FindValidSteadyState
        %
        %   See Also SMPSim, SMPSim.FindValidSteadyState,
        %   SMPSim.AugmentedSteadyState, SMPSim.SS_Soln
            assert(~isempty(obj.u), 'No input vector u specified.')
            assert(~isempty(obj.ts), 'No timing vector ts specified.')
            assert(size(obj.Bs,2) == size(obj.u,1), 'Invalid input vector u specificed')

            if nargin > 1 
                Xss = obj.AugmentedSteadyState(dts);
            else
                Xss = obj.AugmentedSteadyState();
            end
        end
  

        function plotAllStates(obj, fn, oSelect, subplots)
        %plot all states over one period from the current steady-state solution
        %
        %   plotAllStates(obj) will plot all states as seperate subplots of
        %   a single, new figure
        %
        %   plotAllStates(obj,fn) will use the figure number fn
        %
        %   plotAllStates(obj, fn, oSelect) will only plot state numbers
        %   specified by the vector oSelect
        %
        %   plotAllStates(obj, fn, oSelect, subplots) uses the boolean
        %   variable subplots to determine whether to plot each state on
        %   its own subplot or plot them all on a single plot.
        %
        %   See Also SMPSim, SMPSim.SS_WF_Reconstruct, SMPSim.plotWaveforms
            if(nargin <= 1)
                f = figure;
                fn = f.Number;
                subplots = 1;
                oSelect = 1:size(obj.Xs,1); 
            elseif(nargin <= 2)
                subplots = 1;
                oSelect = 1:size(obj.Xs,1); 
            elseif(nargin <=3)
                subplots = 1;
            end
            plotWaveforms(obj, 1, fn, oSelect, subplots);
        end

        function plotAllOutputs(obj, fn, oSelect, subplots)
        %plot all outputs over one period from the current steady-state solution
        %
        %   plotAllOutputs(obj) will plot all outputs as seperate subplots of
        %   a single, new figure
        %
        %   plotAllOutputs(obj,fn) will use the figure number fn
        %
        %   plotAllOutputs(obj, fn, oSelect) will only plot outputs numbers
        %   specified by the vector oSelect.  If oSelect is omitted or left
        %   blank, all original circuit outputs will be included.  If
        %   oSelect is -1 or any vector, the outputs resulting from
        %   state-dependent switch and other measurements will be included.
        %
        %   plotAllOutputs(obj, fn, oSelect, subplots) uses the boolean
        %   variable subplots to determine whether to plot each outputs on
        %   its own subplot or plot them all on a single plot.
        %
        %   See Also SMPSim, SMPSim.SS_WF_Reconstruct,  SMPSim.plotWaveforms   

            no = size(obj.Cs,1);
            [Iname, Vname] = obj.parser.getSwitchMeasSourceNames(obj.switchNames);
            [~,origOuts] = setdiff(obj.outputNames, [Vname, Iname]');
            if(nargin <= 1)
                f = figure;
                fn = f.Number;
                subplots = 1;
                oSelect = origOuts; 
            elseif(nargin <= 2)
                subplots = 1;
                oSelect = origOuts; 
            elseif(nargin <=3)
                subplots = 1;
                if numel(oSelect) == 1 || isempty(oSelect)
                    if isempty(oSelect) ||  oSelect == 0
                        oSelect = origOuts;
                    elseif oSelect == -1
                        oSelect = 1:no;
                    end
                end
            end
            plotWaveforms(obj, 2, fn, oSelect, subplots);
        end
    

        
        function [loc, type] = sigLoc(obj, name, type)
        %Locate specific signal in state, input, output, or switch vectors
        %
        %   loc = sigLoc(obj, name) returns the numerical index loc
        %   into the appropriate vector for the state, input, output, or
        %   switch named name
        %
        %   [loc, type] = sigLoc(obj, name, type) uses additional variable type to
        %   select which vector to search where 
        %       type == 'x' searches the states
        %       type == 'u' searches the inputs
        %       type == 'y' searches the outputs
        %       tpye == 'sw' searches the switches
        %
        %   See Also SMPSim
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
                    obj.topology.switchLabels];
                lengths = [0; cumsum([length(obj.topology.stateLabels); 
                    length(obj.topology.inputLabels); 
                    length(obj.topology.outputLabels)])];
               
                loc = find(strcmp(allLabels, name));
                if ~isempty(loc)
                    vec = find(lengths<loc,1,'last');
                    loc = loc - max(lengths(lengths<loc));
                    switch vec
                        case 1 
                            type = 'x';
                        case 2
                            type = 'u';
                        case 3
                            type = 'y';
                        case 4
                            type = 'sw';
                        otherwise
                            type = [];
                    end
                end
            else
                loc = [];
            end
            if isempty(loc)
                warning('cannot find specified signal in topology');
                type = [];
            end
        end

        function PlossSW = switchLosses(obj, tsteps)
        %Calculate averaged power loss for each switch in steady-state
        %
        %   PlossFET = switchLosses(obj) returns a vector of average power
        %   losses where PlossFET(i) is the loss of the FET corresponding
        %   to obj.switchNamse(i).  The calculated losses are total losses,
        %   including swithcing and conduction as specified by the models
        %   used.
        %
        %   PlossFET = switchLosses(obj, tsteps) also allows the user to
        %   specify the number of time steps in the expansion fo the
        %   waveforms.
        %
        %   See Also SMPSim
            arguments
                obj SMPSim
                tsteps = []
            end
            if ~isempty(tsteps)
                [ ~, t, ys ] = obj.SS_WF_Reconstruct(tsteps);
            else
                [ ~, t, ys ] = obj.SS_WF_Reconstruct();
            end
            for j = 1:numel(obj.switchNames)
                [Iname, Vname] = obj.parser.getSwitchMeasSourceNames(obj.switchNames{j});
                vloc = obj.sigLoc(Vname{1},'y');
                iloc = obj.sigLoc(Iname{1},'y');
                PlossSW(j) = trapz(t,ys(vloc,:).*ys(iloc,:))/sum(obj.ts);
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

        function res = get.swvec(obj)
            res = obj.converter.swvec;
        end

        function res = get.fullu(obj)
            if size(obj.converter.u,3) == 1
                res = repmat(obj.u,1,1,(length(obj.converter.ts)));
            elseif size(obj.converter.u,3) == length(obj.converter.ts)
                res = obj.converter.u;
            elseif size(obj.converter.u,3) == size(obj.converter.fullts,2)
                [~, ints, ~] = obj.converter.getIntervalts();
                res = obj.converter.u(:,:,ints);
            else
                error('invalid input vector u specificed');
            end
        end
        
        function res = get.topology(obj)
            res = obj.converter.topology;
        end

        function res = get.parser(obj)
            res = obj.converter.topology.circuitParser;
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

        function res = get.Ys(obj)
            Xs3 = reshape(obj.Xs,[size(obj.Xs,1), 1, size(obj.Xs,2)]);
            res = squeeze(pagemtimes(obj.Cs,Xs3(:,:,1:end-1)) + pagemtimes(obj.Ds,obj.fullu));
        end

        function res = get.YsEnd(obj)
            Xs3 = reshape(obj.Xs,[size(obj.Xs,1), 1, size(obj.Xs,2)]);
            res = squeeze(pagemtimes(circshift(obj.Cs,1,3),Xs3(:,:,1:end-1)) + pagemtimes(circshift(obj.Ds,1,3),obj.fullu));
            res = circshift(res,-1,2);
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

        function set.swvec(varargin)
            error('Setting swvec is not recommended for class SMPSsim.  Use methods in SMPSconverter');
        end

        function set.parser(obj, newParser)
            obj.converter.topology.circuitParser = newParser;
        end
              
    end

    methods (Hidden)
        %Extending superclass methods and adding Hidden atribute
        % Done so they won't appear and muddy the auto-generated documentation
        function h = findobj(varargin)
            h = findobj@handle(varargin);
        end
        function el = addlistener(varargin)
            el = addlistener@handle(varargin);
        end
        function delete(H)
            delete@handle(H);
        end
        function TF = eq(H1, H2)
            TF = eq@handle(H1,H2);
        end
        function TF = ge(H1, H2)
            TF = ge@handle(H1, H2);
        end
        function TF = gt(H1, H2)
            TF = gt@handle(H1, H2);
        end
        function TF = le(H1, H2)
            TF = le@handle(H1, H2);
        end
        function TF = lt(H1, H2)
            TF = lt@handle(H1, H2);
        end
        function TF = ne(H1, H2)
            TF = ne@handle(H1, H2);
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin);
        end
        function el = listener(varargin)
            el = listener@handle(varargin);
        end
        function notify(varargin)
            notify@handle(varargin);
        end

    end
    
end

