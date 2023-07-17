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
        
        swseq
        
        circuitParser
    end
    
    properties (Dependent)
        constraints
        labels
        
        sourcefn
    end
    
    properties (Hidden)
       % Constraints
        Cbnd
        Dbnd
        bndHyst     %[mean, std]
        switchRef   %[SwitchNo, State its in]
        
        stateLabels
        outputLabels
        inputLabels
        switchLabels 
    end
    
%     properties (SetAccess = private)
%         plecsfn
%         plecsfdate
%     end
    
    properties (Hidden, SetAccess = private)
        expAs
        eigAs
    end
    
    methods
        % From external files
        function loadCircuit(obj, fn, swvec, force)
            
%             if nargin > 2
%                 obj.swseq = [obj.swseq; setdiff(swvec, obj.swseq, 'rows')];
%             end
            if nargin < 4
                force = 0;
            end
            
            if isempty(fn) && ~isempty(obj.circuitParser) && ~isempty(obj.sourcefn)
                fn = obj.sourcefn;
            end
            
            if isempty(obj.circuitParser)
                try %Test if the file is a valid plecs circuit
                    plecs('get', fn, 'StateSpaceOrder');
                    obj.circuitParser = PLECScircuitParser(obj);
                catch e
                end

                if isempty(obj.circuitParser)
                    try
                        %Other Circuit Parsers
                    catch
                    end
                    warning('Circuit is not a PLECS file and alternatives arent implemented yet');
                    rethrow(e);
                end
            end
            loadModel(obj.circuitParser, fn, swvec, force); 
            obj.circuitParser.getConstraintMatrices;
        end
        
        function setSS(obj, As, Bs, Cs, Ds, K)
            obj.As = As;
            obj.Bs = Bs;
            
            if(nargin>3)
                obj.Cs = Cs;
            else
                obj.Cs = zeros(1,size(As,1), size(obj.As,3));
            end
            
            if(nargin>4)
                obj.Ds = Ds;
            else
                obj.Ds = zeros(1, 1, size(obj.As,3));
            end
            
            if(nargin>5)
                obj.K = K;
            else
                obj.K = eye(length(obj.As));
            end
        end
        
        function refreshCache(obj)
            for i = 1:size(obj.As,3)
                obj.expAs(:,:,i) = expm(obj.As(:,:,i));
            end 
        end
        
        %% Setters
        
        function set.As(obj,newAs)
            obj.As = newAs;
            if isempty(newAs)
                obj.eigAs = [];
            else
                for i = 1:size(newAs,3)
                    if isempty(obj.eigAs)
                        obj.eigAs = eigs(newAs(:,:,i));
                    else
                        obj.eigAs(:,:,i) = eigs(newAs(:,:,i));
                    end
    %                 obj.expAs(:,:,i) = expm(newAs(:,:,i));
                end
            end
%             obj.refreshCache;
            % need to notify up the ladder???
            % obj.converter.refreshcache
        end
        
        %% Getters
        
        function expAs = get.expAs(obj)
            expAs = obj.expAs;
        end
        
        function res = get.constraints(obj)
            res.Cbnd = obj.Cbnd;
            res.Dbnd = obj.Dbnd;
            res.bndHys = obj.bndHyst;
            res.switchRef = obj.switchRef;
        end
        
        function res = get.labels(obj)
            res.stateLabels = obj.stateLabels;
            res.outputLabels = obj.outputLabels;
            res.inputLabels = obj.inputLabels;
            res.switchLabels = obj.switchLabels;
        end
        
        function res = get.sourcefn(obj)
            res = obj.circuitParser.sourcefn;
        end

              
    end
    
end

