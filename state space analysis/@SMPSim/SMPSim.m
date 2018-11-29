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
        Converter % Converter
        order
        
        
        %         Aw
        %         Bw
        %         Cw
        %         Dw
        
        Xs
        Xs_circuit
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
        [check] = VfwdIrev(obj)
        [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xi, Si, Vmax, Vmin, progBar)
        [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
        [t] = deadtimecalc(obj,start_value,end_value,state_pos,time_pos)
        
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
            sn = obj.Converter.Topology.stateLabels;
        end
        
        
        %% Test functions
        function loadTestConverter(obj,dotmatfile)
            try
                load(dotmatfile, 'conv');
                obj.Converter = conv;
                %                 params = load(matfile);
            catch err
                ME = MException('resultisNaN:noSuchVariable', ...
                    'Error: test converter file does not contain all requred variables. Required variables are As, Bs, Cs, Ds, ts, and u');
                throw(ME);
            end
            
            obj.settopology(conv.Topology.As, conv.Topology.Bs, conv.Topology.Cs, conv.Topology.Ds);
            obj.setmodulation(conv.ts);
            obj.setinputs(conv.u);
            
            obj.Xs = [];
        end
        
        %% Test functions
        % for TEST_PARSE_SOLN_DIODE
        function loadTestConverter2(obj,conv)
            %{
 try
                load(dotmatfile, 'conv');
                obj.converter = conv;
                %                 params = load(matfile);
            catch err
                ME = MException('resultisNaN:noSuchVariable', ...
                    'Error: test converter file does not contain all requred variables. Required variables are As, Bs, Cs, Ds, ts, and u');
                throw(ME);
            end
            %}
            
            % To ensure matrix is probably statistically not likely to not be
            % invertable (remove dependent states before solve) and puts
            % states in order
            
%             n = size(conv.topology.Parse.OutputNames,1);
%             A = conv.topology.Parse.Anum(:,:,:);
%             B = conv.topology.Parse.Bnum(:,:,:);
%             C = conv.topology.Parse.Cnum(:,:,:);
%             D = conv.topology.Parse.Dnum(:,:,:);
            obj.order = conv.order;
            switchorder = obj.order;
            
            obj.As = [];
            obj.Bs = [];
            obj.Cs = [];
            obj.Ds = [];
            
            obj.As = conv.Topology.Parser.Anum(:,:,switchorder);
            obj.Bs = conv.Topology.Parser.Bnum(:,:,switchorder);
            obj.Cs = conv.Topology.Parser.Cnum(:,:,switchorder);
            obj.Ds = conv.Topology.Parser.Dnum(:,:,switchorder);
            
            %obj.settopology(Ass,Bss,Css,Dss);
            obj.ts = conv.ts;
            obj.u = conv.u;
            obj.oldAs = zeros(size(obj.As));
            obj.oldIntEAt = zeros(size(obj.As));
            obj.Xs = [];
            obj.Converter = conv;
            
            
            %             Ax = conv.topology.Parse.Anum;
            %             Bx = conv.topology.Parse.Bnum;
            %             Cx = conv.topology.Parse.Cnum;
            %             Dx = conv.topology.Parse.Dnum;
            %
            %             obj.Aw = cat(3,Ax(:,:,2),Ax(:,:,1),Ax(:,:,3),Ax(:,:,1));
            %             obj.Bw = cat(3,Bx(:,:,2),Bx(:,:,1),Bx(:,:,3),Bx(:,:,1));
            %             obj.Cw = cat(3,Cx(:,:,2),Cx(:,:,1),Cx(:,:,3),Cx(:,:,1));
            %             obj.Dw = cat(3,Dx(:,:,2),Dx(:,:,1),Dx(:,:,3),Dx(:,:,1));
            
        end
        
    end
    
end

