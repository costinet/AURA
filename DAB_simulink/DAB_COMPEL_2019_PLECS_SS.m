% This is used to run the DAB simulink file and extract the
% eigenvalues

% Jared Baxter


%     _   _   _  ____    _    
%    / \ | | | |/ _  |  / \   
%   / _ \| | | | (_| | / _ \  
%  / ___ | |_| |> _  |/ ___ \ 
% /_/   \_\___//_/ |_/_/   \_\

clear

Cp = 110e-12;
Cs = 1620e-12;
Ls = 76e-6;
Ron_s = 0.0015;
Ron_p = 0.050;
Vg = 48;
Vout = 1.2;
np = 36;
ns = 1;
I_out = 3;
Cout = (4*100+680)*10^-6;
R_out = 0.16;
R_in = 0.005;
Cin = Cout;
L_par_p = 1e-12;
L_par_s = 1e-12;
V_f_p = 1;
V_f_s = 1;
V_f = 1;


C_DC = 100e-6;
R_L = 4.57534;
R_out = 0.141;


order1 = [
    0 1 1 0 0 1 1 0  % REVERSE POWER
    0 0 0 0 0 1 1 0  % primary sw
    1 0 0 1 0 1 1 0  % phase shift
    1 0 0 1 0 0 0 0  % secondary sw
    1 0 0 1 1 0 0 1  % POWER
    0 0 0 0 1 0 0 1  % primary sw
    0 1 1 0 1 0 0 1  % phase shift
    0 1 1 0 0 0 0 0];  % secondary sw

fs = 0.2e6;
Ts = 1/fs;

u = [Vg];

% ts = [Ts*.32 Ts*.01 Ts*.16 Ts*.001 Ts*.32 Ts*.01 Ts*.16 Ts*.001 ];

% original used for eff before conference
% ts = [(Ts/2)-50e-9-400e-9-13.33e-9 50e-9 400e-9 13.33e-9 (Ts/2)-50e-9-400e-9-13.33e-9 50e-9 400e-9 13.33e-9];

P_out = 9;

switch P_out
    case 12
        % 12 W
        pri_sw = 50e-9;
        ps = 556.66e-9;
        sec_sw = 13.33e-9;
        power = Ts/2-pri_sw-ps-sec_sw;
        R_out = 0.141;
        % R_out = 0.12
        % R_out ideal = 0.12
        
    case 9
        % 9 W
        
        pri_sw = 86.66e-9;
        ps = 410e-9;
        sec_sw = 13.33e-9;
        power = Ts/2-pri_sw-ps-sec_sw;
        R_out = 0.16;
         % R_out ideal = 0.16
        
    case 6
        % 6 W
        
        pri_sw = 100e-9;
        ps = 203.33e-9;
        sec_sw = 13.33e-9;
        power = Ts/2-pri_sw-ps-sec_sw;
        R_out = 0.22;
        % R_out ideal = 0.24
        
    case 3
        % 3 W
        
        pri_sw = 136.66e-9;
        ps = 126.6e-9;
        sec_sw = 13.33e-9;
        power = Ts/2-pri_sw-ps-sec_sw;
        R_out = 0.24;
        % R_out ideal = 0.48
end

ts = [power pri_sw ps sec_sw power pri_sw ps sec_sw];


fs_adj = 1/sum(ts);



M14_delay = sum(ts(1:2));
M58_delay = sum(ts(1:4));
M67_delay = 0;
M23_delay = ts(1); %%%% Normally ON *******

M14_ON = sum(ts(3:5))/sum(ts);
M58_ON = sum(ts(5:7))/sum(ts);
M67_ON = sum(ts(1:3))/sum(ts);
M23_ON = sum(ts(2:6))/sum(ts); %%%% Normally ON *******



ModelPath = 'DAB_COMPEL_2019/Circuit8';


%             M1 D1 D2 D3 D4 M2 M3 M4 D5 D6 D7 D8 M5 M6 M7 M8

swVec(1,:) = [1  0  0  0  0  0  0  1  0  0  0  0  1  0  0  1];

ssOrder = plecs('get', ModelPath, 'StateSpaceOrder');
plecs('set', ModelPath, 'SwitchVector', swVec(1,:));
names = plecs('get', ModelPath, 'Topology');
Ass(:,:,1) = names.A;
Bss(:,:,1) = names.B;
Css(:,:,1) = names.C;
Dss(:,:,1) = names.D;
Is(:,:,1) = names.I;

eigA = eig(Ass(:,:,1));

%{

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


% This is like pressing play in Simulink
for i = 1:10
SimOut=sim('DAB_COMPEL_2019');


C1_V = C1sim.data;
C2_V = C2sim.data;
L1_I = L1sim.data;
M1_V = M1sim.data;
M2_V = M2sim.data;
M3_V = M3sim.data;
M4_V = M4sim.data;
M5_V = M5sim.data;
M6_V = M6sim.data;
M7_V = M7sim.data;
M8_V = M8sim.data;

C1_V_data(end+1) = C1sim.data;
C2_V_data(end+1) = C2sim.data;
L1_I_data(end+1) = L1sim.data;
M1_V_data(end+1) = M1sim.data;
M2_V_data(end+1) = M2sim.data;
M3_V_data(end+1) = M3sim.data;
M4_V_data(end+1) = M4sim.data;
M5_V_data(end+1) = M5sim.data;
M6_V_data(end+1) = M6sim.data;
M7_V_data(end+1) = M7sim.data;
M8_V_data(end+1) = M8sim.data;

end


% The one that uses PLECS SS solver (uses setting set in GUI from Simulink)
% tic
% plsteadystate('DAB_COMPEL_2019/Steady-State Analysis');
% toc;
%}
return




