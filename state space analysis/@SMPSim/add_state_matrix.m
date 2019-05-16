function [A,B,C,D] = add_state_matrix(obj,new_states,time_pos)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% DAB Converter
Vg = 48;
L3 = 1e-6;
L2 = 1296e-6;
L1 = 76e-6;
C1 = (4*100+680)*10^-6; %Cout
C2 = C1;
fs = 0.2e6;
Ts = 1/fs;
V = 1.2;
Io = 10; % was 1
M1_C = 110e-12; % CHS
M2_C = 110e-12; % LHS
M3_C = 110e-12; % CHS
M4_C = 110e-12; % LHS
M5_C = 1620e-12; % CHS
M6_C = 1620e-12; % LHS
M7_C = 1620e-12; % CHS
M8_C = 1620e-12; % LHS
D1_C = 1620e-12; % LHS


R2 = 0.01;
R3 = 0.01;

R1 =  0.16; %the out put resistor
dt = Ts/1000;%5e-10;
Vdr = 5;
M1_R_ON = 0.05; % ronHS
M2_R_ON = 0.05; % ronLS
M3_R_ON = 0.05; % ronHS
M4_R_ON = 0.05; % ronLS
M5_R_ON = 0.0015; % ronHS
M6_R_ON = .0015; % ronLS
M7_R_ON = .0015; % ronHS
M8_R_ON = .0015; % ronLS

M1_R_D = 0.05; % ronHS
M2_R_D = 0.05; % ronLS
M3_R_D = 0.05; % ronHS
M4_R_D = 0.05; % ronLS
M5_R_D = .0015; % ronHS
M6_R_D = .0015; % ronLS
M7_R_D = .0015; % ronHS
M8_R_D = .0015; % ronLS



[D1_R_OFF,M1_R_OFF,M2_R_OFF,M3_R_OFF,M4_R_OFF,M5_R_OFF,M6_R_OFF,M7_R_OFF,M8_R_OFF] = deal(100000000); % ronLS

Order = [1 2 3 4 5 6 7 8]; % The order that the states must go in after being parsed
% ts = [Ts*.2376 Ts*.0277 Ts*.23 Ts*.00023 Ts*.2376 Ts*.02775 Ts*.23 Ts*.00023];



%[D1_R,M1_R,M2_R,M3_R,M4_R,M5_R,M6_R,M7_R,M8_R] = deal(100000000); % ronLS



ts = [Ts*.32 Ts*.01 Ts*.16 Ts*.001 Ts*.32 Ts*.01 Ts*.16 Ts*.001];

u = [Vg 1 1 1 1 1 1 1 1]';


TestparseWaveform = false;

ON = 1;
OFF = 0;


SW_OFF = ones(1,8).*10000000;

SW_ON = [0.05 0.05 0.05 0.05 0.0015 0.0015 0.0015 0.0015];

SW_Diode = SW_ON;

SW = [SW_OFF;SW_Diode;SW_ON];


out = {};
orig_new_state = new_states;

k = 1;

% Index fixing to get off = 1 diode = 2 and ON = 3

new_state = new_states(:,time_pos);

new_state(new_state==2)=3; % Turn logic 0 (OFF) to a -1
new_state(new_state==1)=2; % Turn logic 1 (ON) to a 2
new_state(new_state==-1)=1; % Turn logic 0 (OFF) to a -1


The_Codex = obj.Converter.Topology.Parser.Codex;


for i = 1:1:size(The_Codex,2)
     Binary(i,:) = new_state(The_Codex(:,i)); % Assign correct row
end

        for j = 1:1:size(Binary,1)
            out{j} = SW(Binary(j),j);
        end
        [M1_R, M2_R, M3_R, M4_R, M5_R, M6_R, M7_R, M8_R] = deal(out{:});
        
        
        
        HtempAB(:,:,k) = eval(obj.Converter.Topology.Parser.HtempAB(:,:,1));
        HtempCD(:,:,k) = eval(obj.Converter.Topology.Parser.HtempCD(:,:,1));
        dependsAB(:,:,k) = eval(obj.Converter.Topology.Parser.dependsAB(:,:,1));
        savedCD(:,:,k) = eval(obj.Converter.Topology.Parser.savedCD(:,:,1));
        for j = 1:1:size(obj.Converter.Topology.Parser.DependentNames(:,1),1)
            DependentNames(j,k) = eval(obj.Converter.Topology.Parser.DependentNames{j,1});
        end
        for j = 1:1:size(obj.Converter.Topology.Parser.OutputNames(:,1),1)
            OutputNames(j,k) = eval(obj.Converter.Topology.Parser.OutputNames{j,1});
        end
    
    for k = 1:1:size(HtempAB,3)
        [A,B,C,D] = obj.Converter.Topology.Parser.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
        [C,D] = obj.Converter.Topology.Parser.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),obj.Converter.Topology.Parser.SortedTree(:,:,1),obj.Converter.Topology.Parser.SortedCoTree(:,:,1));
    end

  
    B(:,contains(obj.Converter.Topology.Parser.ConstantNames,'VF'))=B(:,contains(obj.Converter.Topology.Parser.ConstantNames,'VF')).*(Binary==2)';
   
    
    

end

