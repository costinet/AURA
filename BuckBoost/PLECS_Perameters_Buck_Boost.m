%{

This code forms the parameters for the buck Boost Converter that
Quillen uses

%}

clear

tic
Vg = 5;
Pout = 25;
Dbuck = 0.92;
Dboost = 0.3118;
M = Dbuck/(1-Dboost);
Vout = M*Vg;
Iout = Pout/Vout;
Rout = 0.73+0.73;
L = 200e-9;
Cout = 23.5e-6;
Cfly = 4.4e-6;
RL = 2.15e-3;
Rg = 0.001;
Vf = 1;


% EPC 2024
Ron = 1.5e-3;
Cds = 1620e-12;


% Duty cycle is in pu
% Phase delay is in seconds
fs = 2000e3;
fs_buck = fs;
Ts = 1/fs;
deadtime = 5e-9;
M1_d = Dbuck-deadtime/Ts;
M1_pd = 0;
M2_d = 1-Dbuck-deadtime/Ts;
M2_pd = Dbuck*Ts;


fs = 1000e3;
fs_boost = fs;
Ts = 1/fs;
deadtime = 5e-9;

M3_d = 1-Dboost-deadtime/Ts;
M3_pd = Dboost*Ts;

% NORMMALLY ONNNNNNNN So its the opposite of the others M4 that is
M4_d = Dboost+deadtime/Ts;
M4_pd = Ts/2-deadtime;

M5_d = Dboost-deadtime/Ts;
M5_pd = Ts/2;

M6_d = Dboost-deadtime/Ts;
M6_pd = 0;




M1_V = 0;
M2_V = 0;
M3_V = 0;
M4_V = 0;
M5_V = 0;
M6_V = 0;
M7_V = 0;
M8_V = 0;
L1_I = 0;
C1_V = 0;
C2_V = 0;

M1_V_data = 0;
M2_V_data = 0;
M3_V_data = 0;
M4_V_data = 0;
M5_V_data = 0;
M6_V_data = 0;
M7_V_data = 0;
M8_V_data = 0;
L1_I_data = 0;
C1_V_data = 0;
C2_V_data = 0;
        
        
        
        
        
        
        sim('BuckBoost'); % This is like pressing play in Simulink
        toc
        C1_V_data(end+1) = C1sim.data(end);
        L1_I_data(end+1) = L1sim.data(end);
        M1_V_data(end+1) = M1sim.data(end);
        M2_V_data(end+1) = M2sim.data(end);
        M3_V_data(end+1) = M3sim.data(end);
        M4_V_data(end+1) = M4sim.data(end);
        M5_V_data(end+1) = M5sim.data(end);
        M6_V_data(end+1) = M6sim.data(end);
        C2_V_data(end+1) = C2sim.data(end);

        
        
        C1_V = C1sim.data(end);
        L1_I = L1sim.data(end);
        M1_V = M1sim.data(end);
        M2_V = M2sim.data(end);
        M3_V = M3sim.data(end);
        M4_V = M4sim.data(end);
        M5_V = M5sim.data(end);
        M6_V = M6sim.data(end);
        C2_V = C2sim.data(end);

        
        
        while abs(C1_V_data(end)-C1_V_data(end-1))>1e-6 && abs(L1_I_data(end)-L1_I_data(end-1))>1e-6
            
            sim('BuckBoost'); % This is like pressing play in Simulink
            
        C1_V_data(end+1) = C1sim.data(end);
        L1_I_data(end+1) = L1sim.data(end);
        M1_V_data(end+1) = M1sim.data(end);
        M2_V_data(end+1) = M2sim.data(end);
        M3_V_data(end+1) = M3sim.data(end);
        M4_V_data(end+1) = M4sim.data(end);
        M5_V_data(end+1) = M5sim.data(end);
        M6_V_data(end+1) = M6sim.data(end);
        C2_V_data(end+1) = C2sim.data(end);

        
        
        C1_V = C1sim.data(end);
        L1_I = L1sim.data(end);
        M1_V = M1sim.data(end);
        M2_V = M2sim.data(end);
        M3_V = M3sim.data(end);
        M4_V = M4sim.data(end);
        M5_V = M5sim.data(end);
        M6_V = M6sim.data(end);
        C2_V = C2sim.data(end);
        end



toc




% 
% ModelPath = 'BuckBoost/Circuit';
% 
% swVec(1,:) = [1 0 0 1 0 0];
% 
% ssOrder = plecs('get', ModelPath, 'StateSpaceOrder');
% plecs('set', ModelPath, 'SwitchVector', swVec(1,:));
% names = plecs('get', ModelPath, 'Topology');
% As(:,:,1) = names.A;
% Bs(:,:,1) = names.B;
% Cs(:,:,1) = names.C;
% Ds(:,:,1) = names.D;
% Is(:,:,1) = names.I;
% 
% eigA = eig(As(:,:,1));




