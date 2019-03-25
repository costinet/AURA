classdef SMPSim < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        LSQoptions = optimoptions('lsqlin','algorithm','trust-region-reflective','Display','none');
        tryOpt = 1;
        condThreshold = 1e9;
        gmin = 1/100e6;
        
        % speedup varaibles -> solution memory
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
    end
    
    properties (Dependent = true)
        ts
        u
    end
    
    properties
        Xs
        converter
    end
    
    methods (Access = private)
        %% Private Methods from external files
        [fresp, intEAt] = forcedResponse(obj, A, expA, B, u, t, storeResult)
        [ valid, newt, dist ] = validByInterval(obj, si, Xs)
        [ x, xdot ] = stateValue_at_t(obj, x0, t, si)
        %AdjustDiodeConduction
        
    end
    
    methods
        %% Methods from external files
        [ Xs] = SS_Soln(obj, Xi, Bi) 
        [ output_args ] = regulate( obj )
        [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
        [ avgXs, avgYs ] = ssAvgs(obj, Xss)
        plotAllStates(obj, fn, subplots)
        plotAllOutputs(obj, fn, subplots)
        
        [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xi, Si, Vmax, Vmin, progBar)
        [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
        
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
            warning('Setting ts is not recommended for class SMPSsim.  Use methods in SMPSconverter');
            obj.converter.ts = newT;
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

        function clearStoredResults(obj)
            obj.oldAs = zeros(size(obj.As));
            obj.oldIntEAt = zeros(size(obj.As));
            obj.oldts = zeros(size(obj.ts));
            obj.Xs = [];
        end    
        
        function sn = getstatenames(obj)
            sn = obj.converter.topology.stateLabels;
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

