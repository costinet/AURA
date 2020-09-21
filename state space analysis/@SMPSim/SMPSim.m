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
        
        eigA
        
        % All of the OG Stuff will be pushed to the converter class
        As_OG
        Bs_OG
        Cs_OG
        Ds_OG
        u_OG
        eigA_OG
        ONorOFF_OG
        ts_OG
        
        ts_history
        Xs_history
        %         Aw
        %         Bw
        %         Cw
        %         Dw
        
        As_saved
        Bs_saved
        Cs_saved
        Ds_saved
        u_saved
        eigA_saved
        ONorOFF_saved
        ts_saved
        Xs_saved
        
        
        
        
        
        Xs
        Xs_circuit % In current implementation this should be Xs if correctXs function is run to ensure the values do not violate KVL and KCL
        
        dead_time_intervals
        dead_time_states
        dead_time_goals
        power_state
        power_goal
        power_
        Vo_index
        Vo_ideal_value
        Perturb1_index
        Perturb2_index
        binary % This is the ON and OFF Position of switches for the converter to be simulated
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
        [Xs] = SS_Soln(obj,keep_SS, As,Bs,ts,u)
        [ xs, t, ys, interval_end ] = SS_WF_Reconstruct(obj, tsteps)
        [ avgXs, avgYs ] = ssAvgs(obj, Xss)
        plotAllStates(obj, fn)
        [check] = VfwdIrev(obj)
        [ ts, dxsdt, hardSwNecessary, multcross, overresonant] = adjustDiodeConduction(obj, Xs, Xi, Si, Vmax, Vmin, progBar)
        [ dXs ] = StateSensitivity(obj, varToPerturb, pI, dX, cI)
        [t] = deadtimecalc(obj,start_value,end_value,state_pos,time_pos)
        [m]=binary_search(obj,A,n,T)
        []=Y_Power(obj)
        [X] = SS_Soln_Aug(obj,keep_SS,As,Bs,ts,u)
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
        
        
        function set.ts(obj,ts)
            if sum(ts>0)==length(ts)
                obj.ts = ts;
            else
                error('There is a non-postitive time interval length trying to be assigned to the simulation class variable ts')
            end
            
            
        end
        
        function setts(obj,ts)
            if sum(ts>0)==length(ts)
                obj.ts = ts;
            else
                error('There is a non-postitive time interval length trying to be assigned to the simulation class variable ts')
            end
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
        
        function sn = getstatenames_Opp(obj)
            sn = obj.Converter.Topology.stateLabels_Opp;
        end
        
        function sn = getstatenames2(obj)
            sn = obj.Converter.Topology.stateLabels;
        end
        
        % History Functions
        
        function ts_hist(obj,ts)
        
           if isempty(obj.ts_history)
               obj.ts_history = ts;
           
           elseif size(obj.ts_history,2) ~= size(ts,2)
                obj.ts_history = [];
               obj.ts_history = ts;
           
           else
               obj.ts_history(end+1,:) = ts;
           end
            
        end    
        
        function Xs_hist(obj,Xs)
        
            if isempty(obj.Xs_history)
               obj.Xs_history = Xs;
           
           elseif size(obj.Xs_history,2) ~= size(Xs,2)
                obj.Xs_history = [];
               obj.Xs_history = Xs;
           
           else
               obj.Xs_history(:,:,end+1) = Xs;
           end
            
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

            obj.order = conv.order;
            switchorder = obj.order;
            
            if max(conv.order)>size(conv.Topology.Parser.Anum,3)
                error('Value of states order given is invalid (exceeds size of A matrix)\n')
            end
            obj.As = [];
            obj.Bs = [];
            obj.Cs = [];
            obj.Ds = [];
            
            obj.As = conv.Topology.Parser.Anum(:,:,switchorder);
            obj.Bs = conv.Topology.Parser.Bnum(:,:,switchorder);
            obj.Cs = conv.Topology.Parser.Cnum(:,:,switchorder);
            obj.Ds = conv.Topology.Parser.Dnum(:,:,switchorder);
            
            if length(conv.ts)~=length(switchorder)
                error('The length of the states and time intervals do not match\n')
            end
            %obj.settopology(Ass,Bss,Css,Dss);
            obj.ts = conv.ts;
            obj.u = conv.u;
            obj.oldAs = zeros(size(obj.As));
            obj.oldIntEAt = zeros(size(obj.As));
            obj.Xs = [];
            obj.Converter = conv;
            
            if size(obj.u,1)~=size(obj.Bs,2)
                error('The size of B and u do not allow matrix multiplication\n')
            end
            
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
        
        
        function updateTestConverter(obj)

            switchorder = obj.order;
            
            obj.As = [];
            obj.Bs = [];
            obj.Cs = [];
            obj.Ds = [];
            
            obj.As = obj.Converter.Topology.Parser.Anum(:,:,switchorder);
            obj.Bs = obj.Converter.Topology.Parser.Bnum(:,:,switchorder);
            obj.Cs = obj.Converter.Topology.Parser.Cnum(:,:,switchorder);
            obj.Ds = obj.Converter.Topology.Parser.Dnum(:,:,switchorder);
            
            obj.oldAs = zeros(size(obj.As));
            obj.oldIntEAt = zeros(size(obj.As));
            
            
            if size(obj.u,1)~=size(obj.Bs,2)
                error('The size of B and u do not allow matrix multiplication\n')
            end
            
            
        end
    end
    
end

