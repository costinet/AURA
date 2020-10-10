function  [simulator] = Run_SS_Converter(conv)
% This fucntion find the steady-state of a converter and outputs the
% simulation parameters

%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

simulator.loadTestConverter2(conv);
simulator.eigA = parse.eigA;
simulator.binary = Binary_for_DAB;
Xss = simulator.SS_Soln();
Xss_Aug=simulator.SS_Soln_Aug();

%% Reconstruction of Dependent variables
parse.StateVarIndex();
simulator.CorrectXs();


for i = 1:1:length(parse.StateNumbers)
    if strcmp(parse.OutputNamesCD{parse.StateNumbers(i),1}(1,end),'A')
        top.stateLabels(end+1,1) = strcat('I_{', parse.StateNames(i,1),'} (A)');
        top.stateLabels_Opp(end+1,1) = strcat('V_{', parse.StateNames(i,1),'} (V)');
    else
        top.stateLabels(end+1,1) = strcat(parse.OutputNamesCD{parse.StateNumbers(i),1}(1,end),'_{', parse.StateNames(i,1),'} (V)');
        top.stateLabels_Opp(end+1,1) = strcat('I_{', parse.StateNames(i,1),'} (A)');
    end
end


parse.find_diode_new(Order,Binary_for_DAB,Diode_Forward_Voltage);
iterations = 50;
%Optimization = SMPSOptim;
%Optimization.Simulator = simulator;
%Optimization.opptimization_loop;
%{
cycle = 0;
toc
fail = simulator.Three_tier_diode_correct(iterations,0,0);
while fail
    
fail = simulator.Three_tier_diode_correct(iterations,0,1);
cycle = cycle+1;
end
%}

toc


StateNumbers = simulator.Converter.Topology.Parser.StateNumbers;
StateNumbers_Opp = simulator.Converter.Topology.Parser.StateNumbers_Opposite;

%% Adjust Power times

    simulator.As_OG = simulator.As;
    simulator.Bs_OG = simulator.Bs;
    simulator.Cs_OG = simulator.Cs;
    simulator.Ds_OG = simulator.Ds;
    simulator.u_OG = simulator.u;
    simulator.eigA_OG = simulator.eigA;
    simulator.ONorOFF_OG = simulator.Converter.Topology.Parser.ONorOFF;
    simulator.ts_OG = simulator.ts;
    
onward = true;

started=tic;
while onward
    onward = false;
    round=tic;
    simulator.Three_tier_diode_correct(50,0,0);
    toc(round)

    
    %%  All of this is for the buck boot converter
    ts = simulator.ts;
    Xs = simulator.Xs;
    As = simulator.As;
    Bs = simulator.Bs;
    u = simulator.u;
    ONorOFF = simulator.Converter.Topology.Parser.ONorOFF;
    
    
%     cumulativesum=cumsum(ts)./sum(ts);
%     findcumulative = find(abs(cumulativesum-0.5)<eps*100);
%     
    half_way_ONofOFF  =  [2 ;    2 ;   -1;     2  ;  -1  ;   0  ;  -1   ;  0   ;  0];
    
    
    ONorOFFtofind =  [ -1 ;    2  ;   2  ;  -1  ;  -1   ;  0  ;   2  ;   0  ;   0];
    
    indextochange= [];
    for i = 1:size(ONorOFF,2)
        if ONorOFF(:,i)== ONorOFFtofind
            indextochange(end+1) = i;
        end
        
        if ONorOFF(:,i)== half_way_ONofOFF
            findcumulative = i-1;
        end
            
    end
    
    Xsindex = 6; %This is to find the inductor row in the steady state solution
    GoaliL = 4.9; % This is the goal minimum inductor current
    Sir=Xsindex;
    

    
    for z = 1:2
        
        if z == 1
            Xic = findcumulative+1;
        elseif z==2
            Xic  = size(Xs,2);
        end
        
        
        
        if (abs(Xs(Xsindex,Xic)-GoaliL)>0.1) 
            onward = true;
            progBar = 0.001;
            maxStep = sum(ts)/100/(.5+1.5*progBar);
            
            
            
            
            Ti = indextochange(z);
            

            
            
            delta_DTs = max(min(ts)/1000, sum(ts)/100000);
            [dXs,delta_DTs] = simulator.Baxter_StateSensitivity2(0, 'ts', Ti, delta_DTs);
            dxsdt = (dXs-Xs)/delta_DTs;
            
            tdelta = (GoaliL-Xs(Sir,Xic))/(dxsdt(Sir,Xic));
            % tdelta = min(max(tdelta, -ts(Ti) + delta_DTs), tsmax(Ti) - ts(Ti));
            tdelta = sign(tdelta)*min(abs(tdelta), maxStep);
            ts(Ti) = ts(Ti) + tdelta;
            if(ts(Ti)<0)
                change = abs(ts(Ti))*0.1;
                ts(Ti) = change;
            end
            
            if z ==1
                simulator.ts_OG(5) = simulator.ts_OG(5)+tdelta;
            end
            
            
            if z ==2
                simulator.ts_OG(11) = simulator.ts_OG(11)+tdelta;
            end
            break
        end
    end
    
    simulator.setts(ts);
    simulator.SS_Soln();
    simulator.CorrectXs();

end
toc(started)


toc
simulator.Three_tier_diode_correct(50,0,0);

return




end