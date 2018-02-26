classdef SMPSim < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        LSQoptions = optimoptions('lsqlin','algorithm','trust-region-reflective','Display','none');
        tryOpt = 1;
        condThreshold = 1e9;
        gmin = 1/100e6;
        
        converter;
        
        % speedup varaibles -> solution memory
        oldAs;
        oldts;
        oldIntEAt;
    end
    
    properties
        As
        Bs
        Cs
        Ds
        ts
        u
        
        Xs
    end
    
    methods (Access = private)
        %% Private Methods from external files
        [fresp, intEAt] = forcedResponse(obj, A, expA, B, u, t, storeResult)
        %AdjustDiodeConduction
    end
    
    methods
        function obj = SMPSim()
            obj.Xs = [];
        end
        
        %% Methods from external files
        [ Xs] = SS_Soln(obj, Xi, Bi)    
        [ xs, t, ys ] = SS_WF_Reconstruct(obj, tsteps)
        [ avgXs, avgYs ] = ssAvgs(obj, Xss)
        plotAllStates(obj, fn)
        
        [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xi, Si, Vmax, Vmin, progBar)
        [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
        
        
        %% Locally-defined methods
        function settopology(obj, As, Bs, Cs, Ds)
           obj.As = As;
           obj.Bs = Bs;
           obj.Cs = Cs;
           obj.Ds = Ds;
           
           obj.oldAs = zeros(size(As));
           obj.oldIntEAt = zeros(size(As));
           
           obj.Xs = [];
        end
        
        function setmodulation(obj, ts)
            obj.ts = ts;
            obj.oldts = zeros(size(ts));

            obj.Xs = [];
        end
        
        function setinputs(obj, u)
            obj.u = u;
            
            obj.Xs = [];
        end
        
        function sn = getstatenames(obj)
            sn = obj.converter.topology.stateLabels;
        end
        
        
        %% Test functions        
        function loadTestConverter(obj,dotmatfile)
            try
                load(dotmatfile, 'conv');
                obj.converter = conv;
%                 params = load(matfile);
            catch err
                ME = MException('resultisNaN:noSuchVariable', ...
                       'Error: test converter file does not contain all requred variables. Required variables are As, Bs, Cs, Ds, ts, and u');
                throw(ME);
            end
            
            obj.settopology(conv.topology.As, conv.topology.Bs, conv.topology.Cs, conv.topology.Ds);
            obj.setmodulation(conv.ts);
            obj.setinputs(conv.u);
            
            obj.Xs = [];
        end
    end
    
end

