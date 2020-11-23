classdef SMPStopology < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        As
        Bs
        Cs
        Ds
        Is % Matrices denoting any dependent states
        K % Matrix of state component values (L's and C's)
        
        % Constraints
        Cbnd
        Dbnd
        
        stateLabels
        outputLabels
        inputLabels
        switchLabels
        
        swseq
        
        constraints
    end
    
    properties (SetAccess = private)
        expAs
    end
    
    methods
        function setSS(obj, As, Bs, Cs, Ds, K)
            obj.As = As;
            obj.Bs = Bs;
            
            if(nargin>3)
                obj.Cs = Cs;
            else
                obj.Cs = zeros(1,size(As,1));
            end
            
            if(nargin>4)
                obj.Ds = Ds;
            else
                obj.Ds = zeros(1,size(Bs,1));
            end
            
            if(nargin>5)
                obj.K = K;
            else
                obj.K = eye(length(obj.As));
            end
        end
        
        function set.As(obj,newAs)
            obj.As = newAs;
            obj.refreshCache;
            % need to notify up the ladder???
            % obj.converter.refreshcache
        end
        
        function refreshCache(obj)
            for i = 1:size(obj.As,3)
                obj.expAs(:,:,i) = expm(obj.As(:,:,i));
            end 
        end
        
        function expAs = get.expAs(obj)
            expAs = obj.expAs;
        end
        
        function setConstraints(obj, constraints)
            if ~isa(constraints,'constraints')
                ME = MException('resultisNaN:noSuchVariable', ...
                       'Error: constraints must be an object of constraints class');
                throw(ME);
            else
                obj.constraints = constraints;
            end
        end
        
        
        
        function loadPLECsModel(obj, fn, swseq)
            try
                plecs('set',fn,'EnableStateSpaceSplitting', 'off');
            catch
                warning('Unable to turn of State Space Splitting in PLECS.  Resulting matrices may not be complete')
            end
            
            open_system(fn(1:strfind(fn,'/')-1),'loadonly');
            ssOrder = plecs('get', fn, 'StateSpaceOrder');

            obj.switchLabels = cellfun(@(x) x(strfind(x,'FET')+3:end), ssOrder.Switches, 'un', 0)';
            obj.stateLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.States, 'un', 0);
            obj.outputLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Outputs, 'un', 0);
            obj.inputLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Inputs, 'un', 0);
            
            if nargin > 2
                obj.swseq = swseq;
                for i = 1:size(swseq,1)
                    plecs('set', fn, 'SwitchVector', swseq(i,:));
                    names = plecs('get', fn, 'Topology');
                    obj.As(:,:,i) = names.A;
                    obj.Bs(:,:,i) = names.B;
                    obj.Cs(:,:,i) = names.C;
                    obj.Ds(:,:,i) = names.D;
                    obj.Is(:,:,i) = names.I;
                end
            end
            
            obj.K = eye(length(obj.stateLabels));
            
            for i = 1:length(obj.stateLabels)
                cloc = strfind(obj.stateLabels{i}, ':');
                if isempty(cloc)
                    element = plecs('get',[fn, '/', obj.stateLabels{i}]);
                else
                    element = plecs('get',[fn, '/', obj.stateLabels{i}(1:cloc-1)]);
                end
                if strcmp(element.Type, 'Capacitor')
                    obj.K(i,i) = evalin('base',element.C);
                elseif strcmp(element.Type, 'Inductor')
                    obj.K(i,i)  = evalin('base',element.L);
                elseif strcmp(element.Type, 'Transformer')
                    obj.K(i,i)  = evalin('base',element.Lm);
                end
            end

        end

              
    end
    
end

