

Vg = 24;
V = 5;
Io = 5;

fs = .5e6;
Ts = 1/fs;

L = 143e-9;
Cout = 100e-6;
C1 =4.7e-6;
C2 = C1;
C3 = C1;
ESR1 = 1.2e-3 + 4e-3;
ESR2 = 1.2e-3 +15e-3;
ESR3 = 1.2e-3 + 4e-3;
ron = 1.2e-3 + 10e-3*(.5*3/4*5/6);
Rl = 5e-3;
Coss = 2.5e-9;
Qg = 18e-9;


%% State Space Construction
% x = [V1 V2 V3 Vout IL]
u = [Vg Io]';

ic3 = [1, -1, -1, 0, 3*ron + ESR1 + ESR2]/(5*ron + ESR1 + ESR2 + ESR3);
ic1 = -ic3 + [0 0 0 0 1]; %iL-ic3
ic2 = -ic1;

A1a = cat(1, ...
ic1, ...
ic2, ...
ic3, ...
[0, 0, 0, 0, 1], ...
-ic3*(2*ron + ESR3) - [0 0 1 1 Rl] );

B1a = cat(1, ...
[-1/(5*ron + ESR1 + ESR2 + ESR3), 0], ...
[1/(5*ron + ESR1 + ESR2 + ESR3), 0], ...
[1/(5*ron + ESR1 + ESR2 + ESR3), 0], ...
[0, -1], ...
[-1/(5*ron + ESR1 + ESR2 + ESR3)*(2*ron + ESR3) + 1, 0]);


Cig1a = ic3;
Dig1a = [1/(5*ron + ESR1 + ESR2 + ESR3), 0];

ic3 = [1, 1, -1, 0, -2*ron + ESR1]/(5*ron + ESR1 + ESR2 + ESR3);
ic1 = -ic3 + [0 0 0 0 -1]; 
ic2 = -ic3;

A2a = cat(1, ...
ic1, ...
ic2, ...
ic3, ...
[0, 0, 0, 0, 1], ...
ic1*(2*ron + ESR1) + [1 0 0 -1 -Rl] );

B2a = cat(1, ...
[0, 0], ...
[0, 0], ...
[0, 0], ...
[0, -1], ...
[0, 0]);


Cig2a = [0 0 0 0 0];
Dig2a = [0 0];


A1b = [0, 0, 0, 0, 1;
  0, 0, 0, 0, -1;
  0, 0, 0, 0, 0;
  0, 0, 0, 0, 1;
  -1, 1, 0, -1, (-3*ron-Rl-ESR1-ESR2)];

B1b = [0, 0;
0, 0;
0, 0;
0, -1;
0, 0];

Cig1b = [0 0 0 0 0];
Dig1b = [0 0];

A2b = [0, 0, 0, 0, 0;
  0, 0, 0, 0, 1;
  0, 0, 0, 0, -1;
  0, 0, 0, 0, 1;
  0, -1, 1, -1, (-3*ron-Rl-ESR3-ESR2)];

B2b = [0, 0;
0, 0;
0, 0;
0, -1;
0, 0];

Cig2b = [0 0 0 0 0];
Dig2b = [0 0];

A3 = [0, 0, 0, 0, 0;
  0, 0, 0, 0, 0;
  0, 0, 0, 0, 0;
  0, 0, 0, 0, 1;
  0, 0, 0, -1, -(ron)-Rl];

B3 = [0, 0;
0, 0;
0, 0;
0, -1;
0, 0];

Cig3 = [0, 0,0, 0,0];
Dig3 = [0, 0];

stateElem = [C1 C2 C3 Cout L];
K = diag(stateElem);


M = 5/6;
phi = 0;%-M*0*Ts;
th = 0;%-M*0.05*Ts;

Da = 3/4;

ts = [Da*Ts/2*M-phi+th, (1-Da)*Ts/2*M+phi+th, Ts/2*(1-M), Da*Ts/2*M-phi-th, (1-Da)*Ts/2*M+phi-th, Ts/2*(1-M)];


As = cat(3, A1a, A1b, A3, A2a, A2b, A3);
Bs = cat(3, B1a, B1b, B3, B2a, B2b, B3);
Cs = cat(3, Cig1a, Cig1b, Cig3, Cig2a, Cig2b, Cig3);
Ds = cat(3, Dig1a, Dig1b, Dig3, Dig2a, Dig2b, Dig3);


for i = 1:size(As,3)
    As(:,:,i) = K^-1*As(:,:,i);
    Bs(:,:,i) = K^-1*Bs(:,:,i);
end

Xi = [Vg/4, Vg/2, 3*Vg/4, 6*M Io];
Bi = [.9 .9 .9 .9 .9];

top = SMPStopology();
top.setSS(As, Bs, Cs, Ds, K);
top.Xi = Xi;
top.Bi = Bi;
top.stateLabels =  {'V_{C1}', 'V_{C2}', 'V_{C3}', 'V_{out}', 'I_L'};
top.outputLabels = {'i_g'};

conv = SMPSconverter();
conv.topology = top;
conv.ts = ts;
conv.u = u;



% save('HybridDickson.mat', 'As','Bs','Cs','Ds','u','ts', 'Xi', 'Bi');
save('HybridDickson.mat', 'conv');
    