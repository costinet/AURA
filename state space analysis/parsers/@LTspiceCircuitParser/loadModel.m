function loadModel(obj, fn, swseq, force)
% loadModel(obj, fn, swseq, force) (LTSpice Implementation)
% fn = path to ltspice model
% swseq = binary vector of switching states to be parsed
% force = binary scalar.  If force == 0, the action won't
% repeat on subsequent calls unless the file has been updated
% or a new swseq is included.
% Required behavior:
% populates fields As, Bs, Cs, Ds, Is, K
%     switchLabels, stateLabels, outputLabels, inputLabels



% %conv = obj.converter;
% top = obj;
% %ts = conv.ts;
% %u = conv.u ;
% Order = obj.order ;
% Numerical_Components = obj.Element_Properties;
% Switch_Resistors = obj.Switch_Resistors;
% SW = obj.Switch_Resistor_Values;
% Switch_Sequence = obj.Switch_Sequence;
% Diode_Forward_Voltage = obj.Fwd_Voltage;
% Binary_for_DAB = Switch_Sequence;
% obj.Component_Values = Numerical_Components;
%%  The above all replaced by LTSpice parsing

% obj.initialize(fn,{},{});
arguments 
    obj LTspiceCircuitParser
    fn {mustBeText} = obj.sourcefn
    swseq {mustBeNumeric} = [];
    force {mustBeScalarOrEmpty} = 0;
end

% if nargin == 1
%     fn = obj.sourcefn;
%     force = 0;
%     swseq = [];
% end

if ~isempty(fn)
    [fp,fn2,EXT] = fileparts(fn);
    if strcmp(EXT,'.asc')
        try
            obj.ascfn = fn;
            
            fn = fullfile(fp,[fn2 '.net']);
            if ~exist(fn,'file') %if no netlist already in existence
                status = system(['"' obj.LTSpiceExe '" -netlist "' fileparts(which(obj.ascfn)) filesep obj.ascfn '"' ]);
                if status == -1
                    msg = 'Unable to execute LTSpice from the command line. This could be a file path or permissions issue';
                    error(msg);
                end
            end
        catch e
            errID = 'LOADMODEL:BadIndex';
            msg = ['Unable to call LTSpice to generate netlist from asc file. '...
                ' Check <a href="matlab:matlab.desktop.editor.openAndGoToLine(''',...
                which('LTspiceCircuitParser'),''',22)">LTspiceCircuitParser</a> parameter LTSpiceExe or supply a .net file directly.'];
            addlException = MException(errID,msg);
            e = addCause(e,addlException);
            throw(e)
        end
    end

    if isempty(obj.Anum) || force || ~strcmp(obj.sourcefn,fn)
        if ~strcmp(obj.sourcefn,fn) || isempty(obj.origComponents) || force
            obj.sourcefn = fn;
            obj.readSpiceNetlist(fn)       
        end
    end
end


if isempty(obj.Anum) || force || ~strcmp(obj.sourcefn,fn)
    if(strcmp(obj.method, 'new'))
        if ~strcmp(obj.sourcefn,fn) || isempty(obj.components) || force
            obj.linearizeCircuitModel2()
        end
    elseif(strcmp(obj.method, 'old'))
        obj.linearizeCircuitModel()
    
        obj.read_file_num();
        obj.addmeasure;
        [switches]=obj.findDM;
        obj.ON_States = cell(length(switches),1);
        obj.OFF_States = cell(length(switches),1);
        [obj.NL,obj.NLnets,~]=obj.Single_states_D(1,1,switches);
        obj.cutset_loop_num();
        obj.Component_Values = obj.Element_Properties;
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

if(strcmp(obj.method, 'old'))
    if size(swseq,2) < length(obj.Switch_Resistors)
        swseq = [swseq, zeros(size(swseq,1), length(obj.Switch_Resistors)-size(swseq,2))];
    end
elseif(strcmp(obj.method, 'new'))
    numSw = sum(strcmp({obj.origComponents.Type}, 'M') + strcmp({obj.origComponents.Type}, 'D'));
     if size(swseq,2) < numSw
         warning('found a switching sequence with fewer switches than the topology contains.  Remaining switches will be set to OFF.')
         swseq = [swseq, zeros(size(swseq,1), numSw-size(swseq,2))];
     end
end

% What is this?  It seems errant -- length(Diodes_POS) is alwasys
% length(Switch_Resistors)
if(strcmp(obj.method, 'old'))
    Diodes_POS = contains(obj.Switch_Resistors,'D')';
    Switch_L = length(Diodes_POS);
end


if ~force
    if nargin > 2
        if ~isempty(obj.topology.swseq) && all(ismember(swseq, obj.topology.swseq, 'rows'))
            return
        else
            %% only new swseqs, just add those
            newSwSeq = setdiff(swseq, obj.topology.swseq, 'rows');
            startLoc = size(obj.Anum,3);
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
        if(strcmp(obj.method, 'new'))
            obj.setSwitchingState(newSwSeq(k,:));
            [A,B,C,D,I] = solveStateSpaceRepresentation(obj);
        elseif(strcmp(obj.method, 'old'))
            [A,B,C,D,I] = obj.ABCD_num(obj.Switch_Resistors,obj.Switch_Resistor_Values,newSwSeq(k,:));
        end
        
        if(strcmp(obj.method, 'old'))
            % Only allows diodes that are on to have non zeros in Bs and
            % Ds
            B_R = size(B,1);
            B_C = size(B,2);
            Diode_adjust_B = ones(B_R,B_C);
            D_Key = newSwSeq(k,:).* Diodes_POS;
            Diode_adjust_B(:,B_C-Switch_L+1:end) = repmat(D_Key,[B_R,1]);
            B = B.*Diode_adjust_B;
            
            D_R = size(D,1);
            D_C = size(D,2);
            Diode_adjust_D = ones(D_R,D_C);
            D_Key = swseq(k,:).* Diodes_POS;
            Diode_adjust_D(:,D_C-Switch_L+1:end) = repmat(D_Key,[D_R,1]);
            D = D.*Diode_adjust_D;
        end
        
        
        
        obj.Anum(:,:,k+startLoc) = A;
        obj.Bnum(:,:,k+startLoc) = B;
        obj.Cnum(:,:,k+startLoc) = C;
        obj.Dnum(:,:,k+startLoc) = D;
        obj.Inum(:,:,k+startLoc) = I;
        [obj.eigA(:,k+startLoc)] = eig(A);
    
        
        
    end

    


    %%% This to account for there being to diodes on at the beginning of
    %%% the run
    B = obj.Bnum;
    B(:,contains(obj.ConstantNames,'C'),:)=0;
    obj.Bnum = B;
    D = obj.Dnum;
    D(:,contains(obj.ConstantNames,'C'),:)=0;
    obj.Dnum = D;
    
    
    obj.topology.As = obj.Anum;
    obj.topology.Bs = obj.Bnum;
    obj.topology.Cs = obj.Cnum;
    obj.topology.Ds = obj.Dnum;
    obj.topology.Is = obj.Inum; %repmat(eye(size(obj.Anum,1)),[1 1 size(obj.Anum,3)]);
end

obj.topology.switchLabels = obj.Switch_Names;
obj.topology.stateLabels = obj.StateNames;
obj.topology.outputLabels = obj.OutputNamesCD;
obj.topology.inputLabels = obj.ConstantNames;

obj.parsedU = zeros(length(obj.topology.inputLabels),1);
[~,IA, IB] = intersect(obj.topology.inputLabels,{obj.components.Name});
obj.parsedU(IA) = [obj.components(IB).paramVals];
%Note, this won't find Vfs, only the sources

end



