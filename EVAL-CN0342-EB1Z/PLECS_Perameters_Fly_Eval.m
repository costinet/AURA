%{

This code forms the parameters for the Active Clamp flyback model
found in the reference listed below:

R. Perrin, N. Quentin, B. Allard, C. Martin and M. Ali,
"High-Temperature GaN Active-Clamp Flyback Converter With Resonant
Operation Mode," in IEEE Journal of Emerging and Selected Topics in
Power Electronics, vol. 4, no. 3, pp. 1077-1085, Sept. 2016.

Modifications have been made to account for a 5 to 20 V 3W power
supply for a gate driver based on the topology.

%}

clear


ShortedL = 0.47382e-6;
ShortedR = 326.9e-3;

Ll = 1.0393e-6-ShortedL;
Lm = 16.228e-6-ShortedL;


Co = (1+10+100+47)*1e-6;
Rds1 = 0.054;

Drec_Ron = 0.0305;
Drec_Vf = 0.2;

R3 = 0.63521-ShortedR;


Dsnub_Ron = 0.2834;
Dsnub_Vf = 0.2;


Rsnub = 390;
Csnub = 220e-9;



Cds1 = 300e-12;
Cds2 = 271.9951e-12;
Cds3 = 28.47e-12;
duty = .52;
np = 1;
ns = 1;

Ro = 17;
Vg = 5;
fs = 200e3;
deadtime = 0;

Ts = 1/fs;

Ns = 1000;
omega = 2*pi()*1.1765e+07;
tic

M1_V = 0;
M2_V = 0;
M3_V = 0;
L1_I = 0;
L2_I = 0;
C1_V = 0;
C2_V = 0;

M1_V_data = 0;
M2_V_data = 0;
M3_V_data = 0;
L1_I_data = 0;
L2_I_data = 0;
C1_V_data = 0;
C2_V_data = 0;

sim('Flyback_PLECS_Model_EVAL'); % This is like pressing play in Simulink
toc
C1_V_data(end+1) = C1sim.data(end);
L1_I_data(end+1) = L1sim.data(end);
L2_I_data(end+1) = L2sim.data(end);
M1_V_data(end+1) = M1sim.data(end);
M2_V_data(end+1) = M2sim.data(end);
M3_V_data(end+1) = M3sim.data(end);
C2_V_data(end+1) = C2sim.data(end);



C1_V = C1sim.data(end);
L1_I = L1sim.data(end);
M1_V = M1sim.data(end);
M2_V = M2sim.data(end);
M3_V = M3sim.data(end);
L2_I = L2sim.data(end);
C2_V = C2sim.data(end);



while abs(C1_V_data(end)-C1_V_data(end-1))>1e-6 && abs(L1_I_data(end)-L1_I_data(end-1))>1e-6
    
    sim('Flyback_PLECS_Model_EVAL'); % This is like pressing play in Simulink
    
    C1_V_data(end+1) = C1sim.data(end);
    L1_I_data(end+1) = L1sim.data(end);
    L2_I_data(end+1) = L2sim.data(end);
    M1_V_data(end+1) = M1sim.data(end);
    M2_V_data(end+1) = M2sim.data(end);
    M3_V_data(end+1) = M3sim.data(end);
    C2_V_data(end+1) = C2sim.data(end);
    
    
    
    C1_V = C1sim.data(end);
    L1_I = L1sim.data(end);
    M1_V = M1sim.data(end);
    M2_V = M2sim.data(end);
    M3_V = M3sim.data(end);
    L2_I = L2sim.data(end);
    C2_V = C2sim.data(end);
end



sim('Flyback_PLECS_Model_EVAL');
toc

% tic
% plsteadystate('SS_Flyback_PLECS_Model_EVAL2/Steady-State Analysis');
% toc;

return

ModelPath = 'Flyback_PLECS_Model_EVAL/Circuit';

% Switch vector is in the order of [Rec_D Snub_D FET_D Fet]
swVec(1,:) = [1 0 0 0];

ssOrder = plecs('get', ModelPath, 'StateSpaceOrder');
plecs('set', ModelPath, 'SwitchVector', swVec(1,:));
names = plecs('get', ModelPath, 'Topology');
As(:,:,1) = names.A;
Bs(:,:,1) = names.B;
Cs(:,:,1) = names.C;
Ds(:,:,1) = names.D;
Is(:,:,1) = names.I;

eigA = eig(As(:,:,1));



