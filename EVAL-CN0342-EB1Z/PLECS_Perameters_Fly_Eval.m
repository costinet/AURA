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


Dsnub_Ron = 0.2834;
Dsnub_Vf = 0.2;


Rsnub = 390;
Csnub = 220e-9;


Rsub = 


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



