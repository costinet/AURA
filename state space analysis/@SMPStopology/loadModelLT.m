function loadModelLT(obj, fn, swseq, force)
% LOADMODELlt (obj, fn, swseq, force) (LTSpice Implementation)
% fn = path to ltspice model (not required here resolved earlier)
% swseq = binary vector of switching states to be parsed
% force = binary scalar.  If force == 0, the action won't
% repeat on subsequent calls unless the file has been updated
% or a new swseq is included.
% Required behavior:
% populates fields As, Bs, Cs, Ds, Is, K
%     switchLabels, stateLabels, outputLabels, inputLabels

%conv = obj.converter;
top = obj;
%ts = conv.ts;
%u = conv.u ;

% Call class variables needed
Order = top.order ;
Numerical_Components = top.Element_Properties;
Switch_Resistors = top.Switch_Resistors;
SW = top.Switch_Resistor_Values;
Switch_Sequence = top.Switch_Sequence;
Diode_Forward_Voltage = top.Fwd_Voltage;
Binary_for_DAB = Switch_Sequence;
obj.circuitParser.Component_Values = Numerical_Components;

Didoes_POS = contains(Switch_Resistors,'D')';
Switch_L = length(Didoes_POS);

if ~force % only add new switching sequences
    
    if nargin > 2
        if all(ismember(swseq, obj.swseq, 'rows'))
            return
        else
            % only new swseqs, just add those
            newSwSeq = setdiff(swseq, obj.swseq, 'rows');
            
            for k = 1:1:size(newSwSeq,1)
                
                [A,B,C,D] = obj.circuitParser.ABCD_num(Switch_Resistors,SW,newSwSeq(k,:));
                
                % Only allows diodes that are on to have non zeros in Bs and
                % Ds
                B_R = size(B,1);
                B_C = size(B,2);
                Diode_adjust_B = ones(B_R,B_C);
                Diode_adjust_B2 = ones(B_R,B_C);
                B_Key = newSwSeq(k,:).* Didoes_POS;
                B_Key = B_Key(Didoes_POS);
                Diode_adjust_B(:,contains(obj.circuitParser.ConstantNames,{'D'})) = repmat(B_Key,[B_R,1]);
                Diode_adjust_B2(:,contains(obj.circuitParser.ConstantNames,{'M'})) = zeros(size(Diode_adjust_B2(:,contains(obj.circuitParser.ConstantNames,{'M'}))));
                B = B.*Diode_adjust_B.*Diode_adjust_B2;

                D_R = size(D,1);
                D_C = size(D,2);
                Diode_adjust_D = ones(D_R,D_C);
                Diode_adjust_D2 = ones(D_R,D_C);
                D_Key = newSwSeq(k,:).* Didoes_POS;
                D_Key = D_Key(Didoes_POS);
                Diode_adjust_D(:,contains(obj.circuitParser.ConstantNames,{'D'})) = repmat(D_Key,[D_R,1]);
                Diode_adjust_D2(:,contains(obj.circuitParser.ConstantNames,{'M'})) = zeros(size(Diode_adjust_D2(:,contains(obj.circuitParser.ConstantNames,{'M'}))));
                D = D.*Diode_adjust_D.*Diode_adjust_D2;



                obj.circuitParser.Anum(:,:,end+1) = A;
                obj.circuitParser.Bnum(:,:,end+1) = B;
                obj.circuitParser.Cnum(:,:,end+1) = C;
                obj.circuitParser.Dnum(:,:,end+1) = D;
                [obj.circuitParser.eigA(:,end+1)] = eig(A);
                obj.swseq  = [obj.swseq; newSwSeq(k,:)];
                
                
            end
        end
    end
end

if force % add everything even old swseq
    
    obj.swseq = swseq;

    % This is sort of a patch there was an issue where forece was casuing
    % the switching vector number to reset in the Anum database. The isse
    % was in the SC Fib in buck/boost mode it would find M7 and M10 on as
    % the 9 and 10 new parsed matrices but would then insert those matrices
    % into at 9 and 10 when moving on to SC mode when it should have been
    % M1 and M5 turning on. 
%     obj.circuitParser.Anum = [];
%     obj.circuitParser.Bnum = [];
%     obj.circuitParser.Cnum = [];
%     obj.circuitParser.Dnum = [];
%     [obj.circuitParser.eigA] = [];


    for k = 1:1:size(swseq,1)
        
        
        % Find state-space 
        [A,B,C,D] = obj.circuitParser.ABCD_num(Switch_Resistors,SW,swseq(k,:));
        
        % Only allows diodes that are on to have non zeros in Bs and
        % Ds

        % Adjusted April 2023: There was a bug where a constant current
        % source would cause wrong columns of B to be zeroed out to account
        % for diode forward voltage

        B_R = size(B,1);
        B_C = size(B,2);
        Diode_adjust_B = ones(B_R,B_C);
        D_Key = swseq(k,:).* Didoes_POS;
        Diode_adjust_B(:,contains(obj.circuitParser.ConstantNames,{'M','D'})) = repmat(D_Key,[B_R,1]);
        B = B.*Diode_adjust_B;
        
        D_R = size(D,1);
        D_C = size(D,2);
        Diode_adjust_D = ones(D_R,D_C);
        D_Key = swseq(k,:).* Didoes_POS;
        Diode_adjust_D(:,contains(obj.circuitParser.ConstantNames,{'M','D'}))  = repmat(D_Key,[D_R,1]);
        D = D.*Diode_adjust_D;
        
        
        obj.circuitParser.Anum(:,:,k) = A;
        obj.circuitParser.Bnum(:,:,k) = B;
        obj.circuitParser.Cnum(:,:,k) = C;
        obj.circuitParser.Dnum(:,:,k) = D;
        [obj.circuitParser.eigA(:,k)] = eig(A);
        
    end
end

%%%%%%

% It should no longer be needed to delete 

%%% This to account for there being to diodes on at the beginning of
%%% the run 
B = obj.circuitParser.Bnum;
B(:,contains(obj.circuitParser.ConstantNames,'C'),:)=0;
obj.circuitParser.Bnum = B;
D = obj.circuitParser.Dnum;
D(:,contains(obj.circuitParser.ConstantNames,'C'),:)=0;
obj.circuitParser.Dnum = D;



% Assign state-space variables
obj.As = obj.circuitParser.Anum;
obj.Bs = obj.circuitParser.Bnum;
obj.Cs = obj.circuitParser.Cnum;
obj.Ds = obj.circuitParser.Dnum;

% Assign Label Names
obj.switchLabels = obj.Switch_Names;
obj.stateLabels = obj.circuitParser.StateNames;
obj.outputLabels = obj.circuitParser.OutputNamesCD;
obj.inputLabels = obj.circuitParser.ConstantNames;



end



%{
    
    % Old code for plecs solver
    try
                plecs('set',fn,'EnableStateSpaceSplitting', 'off');
            catch
                warning('Unable to turn of State Space Splitting in PLECS.  Resulting matrices may not be complete')
            end
            
            if nargin == 3
                force = 0;
            end
            
            file = dir([fn(1:strfind(fn,'/')-1) '.slx']);
            
            
            
            if ~force
                if ~isempty(obj.sourcefn)
                    if strcmp(obj.sourcefn, fn)
                        if strcmp(obj.sourcefdate, file.date)
                            if nargin > 2
                                if all(ismember(swseq, obj.topology.swseq, 'rows'))
                                    return
                                else
                                    %% only new swseqs, just add those
                                    newSwSeq = setdiff(swseq, obj.topology.swseq, 'rows');
                                    for i = 1:size(newSwSeq,1)
                                        plecs('set', fn, 'SwitchVector', newSwSeq(i,:));
                                        names = plecs('get', fn, 'Topology');
                                        obj.topology.As(:,:,end+1) = names.A;
                                        obj.topology.Bs(:,:,end+1) = names.B;
                                        obj.topology.Cs(:,:,end+1) = names.C;
                                        obj.topology.Ds(:,:,end+1) = names.D;
                                        obj.topology.Is(:,:,end+1) = names.I;
                                        obj.topology.swseq = [obj.topology.swseq; newSwSeq(i,:)];
                                    end
                                    return
                                    
                                end
                            else
                                return
                            end
                        end
                    end
                end
            end
            
            open_system(fn(1:strfind(fn,'/')-1),'loadonly');
            ssOrder = plecs('get', fn, 'StateSpaceOrder');

            obj.topology.switchLabels = cellfun(@(x) x(strfind(x,'FET')+3:end), ssOrder.Switches, 'un', 0)';
            obj.topology.stateLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.States, 'un', 0);
            obj.topology.outputLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Outputs, 'un', 0);
            obj.topology.inputLabels = cellfun(@(x) x(strfind(x,'/')+1:end), ssOrder.Inputs, 'un', 0);
            
            if nargin > 2
                obj.topology.swseq = swseq;
                for i = 1:size(obj.topology.swseq,1)
                    plecs('set', fn, 'SwitchVector', swseq(i,:));
                    names = plecs('get', fn, 'Topology');
                    obj.topology.As(:,:,i) = names.A;
                    obj.topology.Bs(:,:,i) = names.B;
                    obj.topology.Cs(:,:,i) = names.C;
                    obj.topology.Ds(:,:,i) = names.D;
                    obj.topology.Is(:,:,i) = names.I;
                end
            end
            
            obj.topology.K = eye(length(obj.topology.stateLabels));
            
            for i = 1:length(obj.topology.stateLabels)
                cloc = strfind(obj.topology.stateLabels{i}, ':');
                if isempty(cloc)
                    element = plecs('get',[fn, '/', obj.topology.stateLabels{i}]);
                else
                    element = plecs('get',[fn, '/', obj.topology.stateLabels{i}(1:cloc-1)]);
                end
                if strcmp(element.Type, 'Capacitor')
                    obj.topology.K(i,i) = evalin('base',element.C);
                elseif strcmp(element.Type, 'Inductor')
                    obj.topology.K(i,i)  = evalin('base',element.L);
                elseif strcmp(element.Type, 'Transformer')
                    obj.topology.K(i,i)  = evalin('base',element.Lm);
                end
            end
            
            obj.sourcefn = fn;
            obj.sourcefdate = file.date;

        end
%}