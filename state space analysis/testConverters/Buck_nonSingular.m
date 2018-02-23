Vg = 5;
L = 230e-9;
Cout = 4040e-9;
fs = 2e6;
Ts = 1/fs;
V = 1.8;
Io = 1;
Cp = 3.4874e-10;

Rl =  .01;
dt = Ts/1000;%50e-9;
Vdr = 5;
Rshunt = 10*270e3;
ron = .05;

ts = [Ts*.5-dt dt Ts*.5-dt dt];

%% Waveform reconstruction - Nonsingular Case
% x = [Vp Il Vo]
u = [Vg Io]';


A1 = [-1/ron -1 0; 1 -Rl -1; 0 1 0];
A2 = [0 -1 0; 1 -Rl -1; 0 1 -1/Rshunt];
A3 = [-1/ron -1 0; 1 -Rl -1; 0 1 0];
A4 = [-1/Rshunt 0 0; 1 -Rl -1; 0 1 -1/Rshunt];

B1 = [1/ron 0; 0 0; 0 -1];
B2 = [0 0; 0 0; 0 -1];
B3 = [0 0; 0 0; 0 -1];
B4 = [0 0; 0 0; 0 -1];

K = [Cp 0 0 ; 0 L 0 ; 0 0 Cout];

C1 = [0 1 0];
D1 = [0 0];
C2 = [0 0 0];
D2 = [0 0];
C3 = C2;
D3 = D2;

As = cat(3, A1, A2, A3, A2);
Bs = cat(3,B1, B2, B3, B2);

Cs = cat(3, C1, C2, C3, C2);
Ds = cat(3, D1, D2, D3, D2);
     
for i = 1:size(As,3)
    As(:,:,i) = K^-1*As(:,:,i);
    Bs(:,:,i) = K^-1*Bs(:,:,i);
end

Xi = [Vg, 0, Io, V];
Bi = [.9 .9 .9 .9];

save('Buck_nonSingular.mat', 'As','Bs','Cs','Ds','u','ts', 'Xi', 'Bi');