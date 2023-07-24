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

%%  The above all needs to be determined by parsing or eval in base workspace

obj.initialize(fn,{},{});
obj.cutset_loop_num();


Didoes_POS = contains(Switch_Resistors,'D')';
Switch_L = length(Didoes_POS);

if ~force
    
    if nargin > 2
        if all(ismember(swseq, obj.swseq, 'rows'))
            return
        else
            %% only new swseqs, just add those
            newSwSeq = setdiff(swseq, obj.topology.swseq, 'rows');
            
            for k = 1:1:size(newSwSeq,1)
                
                [A,B,C,D] = obj.ABCD_num(Switch_Resistors,SW,newSwSeq(k,:));
                
                % Only allows diodes that are on to have non zeros in Bs and
                % Ds
                B_R = size(B,1);
                B_C = size(B,2);
                Diode_adjust_B = ones(B_R,B_C);
                D_Key = newSwSeq(k,:).* Didoes_POS;
                Diode_adjust_B(:,B_C-Switch_L+1:end) = repmat(D_Key,[B_R,1]);
                B = B.*Diode_adjust_B;
                
                D_R = size(D,1);
                D_C = size(D,2);
                Diode_adjust_D = ones(D_R,D_C);
                D_Key = swseq(k,:).* Didoes_POS;
                Diode_adjust_D(:,D_C-Switch_L+1:end) = repmat(D_Key,[D_R,1]);
                D = D.*Diode_adjust_D;
                
                
                
                obj.Anum(:,:,end+1) = A;
                obj.Bnum(:,:,end+1) = B;
                obj.Cnum(:,:,end+1) = C;
                obj.Dnum(:,:,end+1) = D;
                [obj.eigA(:,end+1)] = eig(A);
                obj.swseq  = [obj.swseq; newSwSeq(k,:)];
                
                
            end
        end
    end
end

if force
    
    obj.topology.swseq = swseq;
    for k = 1:1:size(swseq,1)
        
        
        
        [A,B,C,D] = obj.ABCD_num(Switch_Resistors,SW,swseq(k,:));
        
        % Only allows diodes that are on to have non zeros in Bs and
        % Ds
        B_R = size(B,1);
        B_C = size(B,2);
        Diode_adjust_B = ones(B_R,B_C);
        D_Key = swseq(k,:).* Didoes_POS;
        Diode_adjust_B(:,B_C-Switch_L+1:end) = repmat(D_Key,[B_R,1]);
        B = B.*Diode_adjust_B;
        
        D_R = size(D,1);
        D_C = size(D,2);
        Diode_adjust_D = ones(D_R,D_C);
        D_Key = swseq(k,:).* Didoes_POS;
        Diode_adjust_D(:,D_C-Switch_L+1:end) = repmat(D_Key,[D_R,1]);
        D = D.*Diode_adjust_D;
        
        
        obj.Anum(:,:,k) = A;
        obj.Bnum(:,:,k) = B;
        obj.Cnum(:,:,k) = C;
        obj.Dnum(:,:,k) = D;
        [obj.eigA(:,k)] = eig(A);
        
    end
end

%%%%%%



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


obj.topology.switchLabels = obj.Switch_Names;
obj.topology.stateLabels = obj.StateNames;
obj.topology.outputLabels = obj.OutputNamesCD;
obj.topology.inputLabels = obj.ConstantNames;



end



