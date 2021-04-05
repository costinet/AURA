function EVAL_SS_Converter_Fib(conv)



ts = conv.ts;
u = conv.u ;
Order = conv.order ;
Numerical_Components = conv.Element_Properties;
Switch_Resistors = conv.Switch_Resistors;
SW = conv.Switch_Resistor_Values;
Switch_Sequence = conv.Switch_Sequence;
Diode_Forward_Voltage = conv.Fwd_Voltage;
parse = conv.Topology.Parser;
top = conv.Topology;
simulator = SMPSim();
Binary_for_DAB = Switch_Sequence;

%{
% Switch Resistors Contains the same order of resistors as was given
% in the Binary matrix provided by the user to determine when states
% are on or off
for i = 1:1:size(Switch_Resistors,1)
    eval([Switch_Resistors{i} '=' 'out{i}']);
end
%}


for i = 1:1:size(Numerical_Components,1)
    eval([Numerical_Components{i,1} '=' 'Numerical_Components{i,2};'])
end


A = parse.Asym;
B = parse.Bsym;
C = parse.Csym;
D = parse.Dsym;

SortedTree = parse.SortedTree;
SortedCoTree = parse.SortedCoTree;

if isempty(A)
    for k = 1:1:size(Binary_for_DAB,1)
        
        out = {};
        
        for j = 1:1:size(Binary_for_DAB,2)
            out{j} = SW(Binary_for_DAB(k,j)+1,j);
        end
        for i = 1:1:size(Switch_Resistors,1)
            eval([Switch_Resistors{i} '=' 'out{i};'])
        end
        
        HtempAB(:,:,k) = eval(parse.HtempAB(:,:,1));
        HtempCD(:,:,k) = eval(parse.HtempCD(:,:,1));
        
        if ~isempty(parse.dependsAB)
            dependsAB(:,:,k) = eval(parse.dependsAB(:,:,1));
            savedCD(:,:,k) = eval(parse.savedCD(:,:,1));
            for j = 1:1:size(parse.DependentNames(:,1),1)
                DependentNames(j,k) = eval(parse.DependentNames{j,1});
            end
        else
            dependsAB = [];
            savedCD = [];
            DependentNames = [];
        end
        for j = 1:1:size(parse.OutputNames(:,1),1)
            OutputNames(j,k) = eval(parse.OutputNames{j,1});
        end
    end
    for k = 1:1:size(HtempAB,3)
        [A,B,C,D] = parse.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
        [C,D] = parse.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),SortedTree(:,:,1),SortedCoTree(:,:,1));
        parse.Anum(:,:,k)=A;
        parse.Bnum(:,:,k)=B;
        parse.Cnum(:,:,k)=C;
        parse.Dnum(:,:,k)=D;
        [parse.eigA(:,k)] = eig(parse.Anum(:,:,k));
    end
    
    % Set all diode forward voltages to be off
    B = parse.Bnum;
    B(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Bnum = B;
    D = parse.Dnum;
    D(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Dnum = D;
    
%elseif ~isempty(parse.Anum)
%    fprintf('Confirm Plecs used');
else
    for k = 1:1:size(Binary_for_DAB,1)
        
        out = {};
        
        for j = 1:1:size(Binary_for_DAB,2)
            out{j} = SW(Binary_for_DAB(k,j)+1,j);
        end
        
        for i = 1:1:size(Switch_Resistors,1)
            eval([Switch_Resistors{i} '=' 'out{i};']);
        end
        
        parse.Anum(:,:,k) = eval(A(:,:,1));
        parse.Bnum(:,:,k) = eval(B(:,:,1));
        parse.Cnum(:,:,k) = eval(C(:,:,1));
        parse.Dnum(:,:,k) = eval(D(:,:,1));
        [parse.eigA(:,k)] = eig(parse.Anum(:,:,k));
    end
    % Set all diode forward voltages to be off
    B = parse.Bnum;
    B(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Bnum = B;
    D = parse.Dnum;
    D(:,contains(parse.ConstantNames,'VF'),:)=0;
    parse.Dnum = D;
end


end