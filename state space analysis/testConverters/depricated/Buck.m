Vg = 5;
L = 230e-9;
Cout = 4040e-9;
fs = 2e6;
Ts = 1/fs;
V = 1.8;
Io = 1;
CHS = 3.4874e-10;
CLS = 3.4874e-10;

Rl =  .01;
dt = Ts/4;%50e-9;
Vdr = 5;
ronHS = .05;
ronLS = .05;

ts = [Ts*.5-dt dt Ts*.5-dt dt];


%% Waveform reconstruction
% x = [Vhs Vls Il Vo]
u = [Vg Io]';

A1 = [-1/ronHS, 0, 0, 0;
    0, -1/ronHS, -1, 0;
    0, 1, -Rl, -1;
    0, 0, 1, 0];
B1 = [0, 1/ronHS, 0, 0;
    0, 0, 0, -1]';
C1 = [-1/ronHS, -1/ronHS, 0, 0];
D1 = [1/ronHS, 0];

A2 = [0, 0, 1/2, 0;
    0, 0, -1/2, 0;
    0, 1, -Rl, -1;
    0, 0, 1, 0];
B2 = [0, 0, 0, 0;
    0, 0, 0, -1]';
C2 = [0, 0, 1/2, 0];
D2 = [0, 0];

A3 = [-1/ronLS, 0, 1, 0;
    0, -1/ronLS, 0, 0;
    0, 1, -Rl, -1;
    0, 0, 1, 0];
B3 = [1/ronLS, 0, 0, 0;
    0, 0, 0, -1]';
C3 = [-1/ronLS, 0, 1, 0];
D3 = [1/ronLS, 0];




K = [CHS 0 0 0; 0 CLS 0 0; 0 0 L 0 ; 0 0 0 Cout]; 

As = cat(3, A1, A2, A3, A2);
Bs = cat(3, B1, B2, B3, B2);
Cs = cat(3, C1, C2, C3, C2);
Ds = cat(3, D1, D2, D3, D2);
     
for i = 1:size(As,3)
    As(:,:,i) = K^-1*As(:,:,i);
    Bs(:,:,i) = K^-1*Bs(:,:,i);
end

Xi = [Vg, 0, Io, V];
Bi = [.9 .9 .9 .9];

% save('Buck.mat', 'As','Bs','Cs','Ds','u','ts', 'Xi', 'Bi');

% x = [Vhs Vls Il Vo]
Clim = [1, 0, 0, 0;
        0, 1, 0, 0];
Dmax = [-1, 0;
        -1, 0];
Dmin = [0, 0;
        0, 0];
    
Maxlim = [1, 1]';
Minlim = [-1, -1]';


Creg = [0 0 0 1];
Dreg = [0 0];
regtarget = [V];

const = constraints();
const.Cmax = Clim;
const.Dmax = Dmax;
const.Maxlim = Maxlim;
const.Cmin = Clim;
const.Dmin = Dmin;
const.Minlim = Minlim;
const.Creg = Creg;
const.Dreg = Dreg;
const.regtarget = regtarget;

top = SMPStopology();
top.setSS(As, Bs, Cs, Ds);
top.Xi = Xi;
top.Bi = Bi;
top.stateLabels =  {'V_{Chs}', 'V_{Cls}', 'i_L', 'V_{out}'};
top.outputLabels = {'i_g'};
top.setConstraints(const)

conv = SMPSconverter();
conv.topology = top;
conv.ts = ts;
conv.tsmax = [Ts dt Ts dt];
conv.tcomp = [3 3 1 1];
conv.u = u;

save('Buck.mat', 'conv');