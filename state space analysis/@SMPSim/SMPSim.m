classdef SMPSim < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        LSQoptions = optimoptions('lsqlin','algorithm','trust-region-reflective','Display','none');
        tryOpt = 1;
        condThreshold = 1e11;
        
        % speedup varaibles -> solution memory
        oldAs;
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
            
            obj.Xs = [];
        end
        
        function setinputs(obj, u)
            obj.u = u;
            
            obj.Xs = [];
        end
        
        
        %% Test functions        
        function loadTestConverter(obj,dotmatfile)
            try
                load(dotmatfile, 'As','Bs','Cs','Ds','u','ts');
%                 params = load(matfile);
            catch err
                ME = MException('resultisNaN:noSuchVariable', ...
                       'Error: test converter file does not contain all requred variables. Required variables are As, Bs, Cs, Ds, ts, and u');
                throw(ME);
            end
            
            obj.settopology(As, Bs, Cs, Ds);
            obj.setmodulation(ts);
            obj.setinputs(u);
            
            obj.Xs = [];
        end
    end
    
end

