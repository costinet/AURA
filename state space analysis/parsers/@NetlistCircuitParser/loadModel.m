function loadModel(obj, fn, swseq, force)
% loadModel(obj, fn, swseq, force) (Netlist Implementation)
% fn = path to netlist file (*.asc or *.net), or a string netlist
% swseq = binary vector of switching states to be parsed
% force = binary scalar.  If force == 0, the action won't
% repeat on subsequent calls unless the file has been updated
% or a new swseq is included.
% Required behavior:
% populates fields As, Bs, Cs, Ds, Is, K
%     switchLabels, stateLabels, outputLabels, inputLabels

arguments
    obj circuitParser
    fn {mustBeText} = '';
    swseq {mustBeNumericOrLogical} = []
    force {mustBeNumericOrLogical} = 0
end
% obj.initialize(fn,{},{});

if nargin == 1
    fn = obj.sourcefn;
    force = 0;
    swseq = [];
end

if ~isempty(fn)
    [fp,fn2,EXT] = fileparts(fn);
    if strcmp(EXT,'.asc')
        assert(exist(fn,'file'), 'Supplied filename does not exist or is not on the current path');
        try
            obj.ascfn = fn;
            fn = fullfile(fp,[fn2 '.net']);
            if ~exist(fn,'file') || datetime(dir(fn).date) < datetime(dir(obj.ascfn).date)  %if no netlist already in existence, or the asc file has been modified recently
                if isempty(fileparts(which(obj.ascfn)))
                    status = system(['"' obj.LTSpiceExe '" -netlist "' obj.ascfn '"' ]);
                else
                    status = system(['"' obj.LTSpiceExe '" -netlist "' fileparts(which(obj.ascfn)) filesep obj.ascfn '"' ]);
                end
                if status == -1 || ~exist(fn,'file')
                    msg = ['Unable to execute LTSpice from the command line. This could be a file path or permissions issue.  ' ...
                        'Make sure NetlistCircuitParser().LTSpiceExe is the correct path to the LTSpice executable' ...
                        ' e.g. ' newline newline 'NP = NetlistCircuitParser(); NP.LTSpiceExe = ''C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe'''];
                    error(msg);
                end
            end
        catch e
            
            msg = ['Unable to call LTSpice to generate netlist from asc file. '...
                ' Check <a href="matlab:matlab.desktop.editor.openAndGoToLine(''',...
                which('NetlistCircuitParser'),''',22)">NetlistCircuitParser</a> parameter LTSpiceExe or supply a .net file directly.'];
            % errID = 'LOADMODEL:BadIndex';
            % disp(msg);
            % addlException = MException(errID,msg);
            % e = addCause(e,addlException);
            % throw(e)
            error(msg);
        end
    end

    if isempty(obj.topology.As) || force || ~strcmp(obj.sourcefn,fn)
        if ~strcmp(obj.sourcefn,fn) || isempty(obj.origComponents) || force
            obj.readSpiceNetlist(fn)       
        end
    end
end


if isempty(obj.topology.As) || force || ~strcmp(obj.sourcefn,fn)
    if ~strcmp(obj.sourcefn,fn) || isempty(obj.components) || force
        obj.linearizeCircuitModel()
    end
end


if ~exist('swseq','var') || isempty(swseq)
    try
        swseq = evalin('base', 'swvec');
    catch 
        % None defined, just use all FETS off
        numSw = sum(strcmp({obj.origComponents.Type}, 'M') + strcmp({obj.origComponents.Type}, 'D'));
        swseq = zeros(1,numSw);
    end       
end


numSw = sum(strcmp({obj.origComponents.Type}, 'M') + strcmp({obj.origComponents.Type}, 'D'));
 if size(swseq,2) < numSw
     warning('found a switching sequence with fewer switches than the topology contains.  Remaining switches will be set to OFF.')
     swseq = [swseq, zeros(size(swseq,1), numSw-size(swseq,2))];
 end





if ~force
    if nargin > 2
        if ~isempty(obj.topology.swseq) && all(ismember(swseq, obj.topology.swseq, 'rows'))
            return
        else
            %% only new swseqs, just add those
            newSwSeq = setdiff(swseq, obj.topology.swseq, 'rows');
            startLoc = size(obj.topology.As,3);
            obj.topology.swseq  = [obj.topology.swseq; newSwSeq];
        end
    elseif nargin <= 2
        obj.topology.swseq = swseq;
        newSwSeq = swseq;
        startLoc = 0;
    end
else
    obj.topology.swseq = swseq;
    newSwSeq = swseq;
    startLoc = 0;
end

everythingDefined = 1;
if ~isempty(obj.undefinedExpressions)
    try
        for i = 1:size(obj.undefinedExpressions,1)
            evalin('base',obj.undefinedExpressions(i,3))
            obj.undefinedExpressions(i,:) = {};
        end
    catch
        warning(['Undefined value ' obj.undefinedExpressions{i,3} ' for component ' obj.undefinedExpressions{i,1}])
        everythingDefined = 0;
    end
    if everythingDefined == 1
        obj.updateComponentValues();
    end
end

if everythingDefined == 1


    for k = 1:1:size(newSwSeq,1)

        obj.setSwitchingState(newSwSeq(k,:));
        [A,B,C,D,I] = solveStateSpaceRepresentation(obj);

        if ~all(size(I) == size(A))
            % when there is a state-source dependence, the above function
            % will return I and ns x (ns+ni) to capture the source
            % dependence.  Split that up into state dependence I and source
            % dependence BI
            BI = I(1:size(A,1),size(A,2)+1:end);
            I = I(1:size(A,1),1:size(A,1));
        else
            BI = zeros(size(B));
        end
        
        if k+startLoc == 1 && isempty(obj.topology.As)
            obj.topology.As = A;
            obj.topology.Bs = B;
            obj.topology.Cs = C;
            obj.topology.Ds = D;
            obj.topology.Is = I; 
            obj.BInum = BI; %BIs stays here for now.  
            % Should eventually move to topology, but will have to address
            % PLECs compatability.
        else
            obj.topology.As(:,:,k+startLoc) = A;
            obj.topology.Bs(:,:,k+startLoc) = B;
            obj.topology.Cs(:,:,k+startLoc) = C;
            obj.topology.Ds(:,:,k+startLoc) = D;
            obj.topology.Is(:,:,k+startLoc) = I; 
            obj.BInum(:,:,k+startLoc) = BI; %BIs stays here for now.  
            % Should eventually move to topology, but will have to address
            % PLECs compatability.
        end

    end

    obj.getConstraintMatrices;
end


try
    if isempty(obj.topology.converter.u)
        parsedU = zeros(length(obj.topology.inputLabels),1);
        [~,IA, IB] = intersect(obj.topology.inputLabels,{obj.components.Name});
        parsedU(IA) = [obj.components(IB).paramVals];

        % assignin('base','us', parsedU);   %Gradually moving away from base assignment.
        obj.topology.converter.u = parsedU;
    end
    %Note, this won't find Vfs, only the sources
catch e
    %If we can't, the input isn't defined (ok).
end

end



