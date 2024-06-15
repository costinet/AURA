classdef SMPStopology < handle
    %Topology object for use with AURA
    %   SMPStopology is contained within the SMPSConverter class.  SMPStopology
    %   contains methods and handles information related to topological 
    %   description of the circuit.  In most cases, SMPStopology is
    %   automatically constructed when SMPSim() is constructed.
    %
    %   See also SMPSim, SMPSconverter, AURAdb, @circuitParser
    
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
    
   
    properties (Hidden, SetAccess = private)
        expAs
        eigAs

        converter
    end
    
    methods
        function result = loadCircuit(obj, fn, swvec, force)
        %LOADCIRCUIT loads plecs or LTspice circuit 
        %
        %   loadCircuit(obj, fn)
        %   
        %   loadCircuit(obj, fn, swvec)
        %
        %   loadCircuit(obj, [], swvec)
        %   if a topology has previously been loaded, the filename may be
        %   omitted and the function will use the SMPStopology's stored
        %   file as reference for a new swvec
        %
        %   loadCircuit(___, force)
        %   setting force=1 will drop all previously parsed switching
        %   states and repopulate from the supplied swvec

            if nargin < 4
                force = 0;
            end

            if nargin == 2
                swvec = [];
            end
            
            if isempty(fn) && ~isempty(obj.circuitParser) && ~isempty(obj.sourcefn)
                fn = obj.sourcefn;
            end

            if isempty(fn)
                error('loadCircuit requires a file defining the circuit structure');
            end

            % Throw error if file doesn't exist
            if ~exist(fn,'file')
                % if simulink, check if referencing subcircuit in model
                subfn = split(fn,'/');
                subfn = [subfn{1:end-1}];
                if ~exist(subfn,'file')
                    error(['Supplied filename ' fn ' does not exist, is not on the path, or is not readable ']);
                end
            end
            
            if isempty(obj.circuitParser)
                try %Test if the file is a valid PLECS circuit
                    % If it is a simulink-embedded plecs circuit
                    parts = split(fn, '/');
                    modelFile = join(parts{1:end-1}, '/');
                    if exist(modelFile,'file') == 4
                        plecs('get', fn, 'StateSpaceOrder');
                        obj.circuitParser = PLECScircuitParser(obj);
                    end
                catch e 
                    % allow to continue trying other parsers
                    if startsWith(e.message, "Undefined function 'plecs' ") 
                        e = MException('LOADCIRCUIT:plecsNotInstalled','supplied file appears to be a PLECS circuit but plecs is not installed or is not on the path ');
                        % e2 = addCause(e,addlInfo);
						throw(e)
                    elseif startsWith(e.message, 'Error evaluating parameter') 
                        obj.circuitParser = PLECScircuitParser(obj);
                        obj.circuitParser.readOnLoadError(e.message, fn);
                        obj.circuitParser.inspect();
                        obj.circuitParser = [];
                        result = 0;
                        return
                    elseif startsWith(e.message, 'State/source dependence in components')
                        addlInfo = MException('LOADCIRCUIT:badCircuitConfiguration','Circuit has a C-V loop or LI-net');
                        e2 = addCause(e,addlInfo);
						throw(e2)
                    elseif startsWith(e.message, 'Invalid component path')
                         addlInfo = MException('LOADCIRCUIT:cannotOpenFile','Unable to access PLECS file.  Try opening the PLECs model directly in your browser before calling loadCircuit');
                        e2 = addCause(e,addlInfo);
						throw(e2)
                    end
                end

                if isempty(obj.circuitParser)
                    try
                        [~,~,EXT] = fileparts(fn);
                        if strcmp(EXT,'.net') && exist(fn,'file')
						    % obj.circuitParser = LTspiceCircuitParser(obj);
                            obj.circuitParser = NetlistCircuitParser(obj);
                        elseif strcmp(EXT,'.asc') && exist(fn,'file')
                            % obj.circuitParser = LTspiceCircuitParser(obj);
                            obj.circuitParser = NetlistCircuitParser(obj);
                        else
                            error('Circuit cannot be parsed as either a PLECS or LTSpice file')
                        end
                    catch e
                        addlInfo = MException('LOADCIRCUIT:BadInputFile','Circuit cannot be parsed as either a PLECS or LTSpice file');
                        e2 = addCause(e,addlInfo);
						throw(e2)
                    end
                    
                end
            end
            loadModel(obj.circuitParser, fn, swvec, force); 
            obj.circuitParser.getConstraintMatrices;
            result = 1;
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
%             if isempty(newAs)
            obj.eigAs = [];
%             else
%                 for i = 1:size(newAs,3)
%                     if isempty(obj.eigAs)
%                         obj.eigAs = eigs(newAs(:,:,i));
%                     else
%                         obj.eigAs(:,:,i) = eigs(newAs(:,:,i));
%                     end
%     %                 obj.expAs(:,:,i) = expm(newAs(:,:,i));
%                 end
%             end
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

        function res = get.swseq(obj)
            if isempty(obj.swseq)
                try
                    numSw = startsWith({obj.circuitParser.origComponents.Type}, {'D','M'});
                catch 
                    numSw = [];
                end
                if ~isempty(numSw)
                    res = double.empty(0,sum(numSw));
                else
                    res = [];
                end
            else
                res = obj.swseq;
            end
        end

        function res = get.eigAs(obj)
            if ~isempty(obj.eigAs)
                res = obj.eigAs;
            else
                res = zeros(size(obj.As,1),1,size(obj.As,3));
                 for i = 1:size(obj.As,3)
                    evs = eigs(obj.As(:,:,i));
                    res(1:length(evs),:,i) = eigs(obj.As(:,:,i));
                 end
                obj.eigAs = res;
            end
        end

        %% Setters
        function set.circuitParser(obj,parser)
            arguments
                obj SMPStopology
                parser
            end
            obj.circuitParser = parser;
            parser.linkTopology(obj);
        end
              
    end

    methods (Hidden)
        function linkConverter(obj, conv)
            % for making link bidrectional when converter adds topology.
            obj.converter = conv;
        end
    end
    
end

