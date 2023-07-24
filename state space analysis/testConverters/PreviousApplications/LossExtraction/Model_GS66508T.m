load 'results_300V_1.5A_D_step.mat';
% I_L V_DS V_GS V_OUT

useCw = 0;

Vg = 280;
L = 23e-6;
Cout = 440e-9;
fs = 2.76e6;
Ts = 1/fs;
V = 77.5;
Io = 1.65;

Rl =  5;
Rl_400k = .05; 
Rl_800k = .05;
dt = 60e-9;
Vdr = 5;
Rshunt = 270e3;

Rw = 1000;
Cw = 1;

if(useCw == 0)
    Cw = [];
end

D = [.5 .75];

addpath('GS66508T');
GS66508T_digitized;

load 'intVds_off'
col = find(xq(1,:) > Vg,1,'first');
vtest = xq(:,col);
itest = yq(:,col);
lambda_off = zq(:,col);

load 'intVds_on'
col = find(xq(1,:) > Vg,1,'first');
vtest = xq(:,col);
itest = yq(:,col);
lambda_on = zq(:,col);


lambda = 7.2470e-07;
tswoff = lambda*2/Vg;


ron = .05;
qg = 5.8e-9;
Rja = 24;
IsdVsd = Isd_Vsd_Vgs0V;
Vbr = 650;


%% Waveform reconstruction
% x = [Vp Il Vo]
C = Cout;

%% Nominal
% A1 = [-1/ron -1 0; 1 0 -1; 0 1 0];
% A2 = [0 -1 0; 1 0 -1; 0 1 -1/Rshunt];
% A3 = [-1/ron -1 0; 1 0 -1; 0 1 0];
% A4 = A2;
% 
% B1 = [1/ron 0; 0 0; 0 -1];
% B2 = [0 0; 0 0; 0 -1];
% B3 = [0 0; 0 0; 0 -1];
% B4 = B2;

% With RL
if(useCw == 0)
    A1 = [-1/ron -1 0; 1 -Rl -1; 0 1 0];
    A2 = [0 -1 0; 1 -Rl -1; 0 1 -1/Rshunt];
    A3 = [-1/ron -1 0; 1 -Rl -1; 0 1 0];
%     A4 = A2;
    A4 = [-1/Rshunt 0 0; 1 -Rl -1; 0 1 -1/Rshunt];

    B1 = [1/ron 0; 0 0; 0 -1];
    B2 = [0 0; 0 0; 0 -1];
    B3 = [0 0; 0 0; 0 -1];
%     B4 = B2;
    B4 = [1/tswoff 0; 0 0; 0 -1];
    
else
%% With RL, Cw, Rw
    A1 = [(-1/ron-1/Rw) -1 1/Rw 1/Rw; 
           1 -Rl -1 0;
           1/Rw 1 (-1/Rw-1/Rshunt) -1/Rw;
           1/Rw 0 -1/Rw -1/Rw];
    A2 = [(-1/Rw) -1 1/Rw 1/Rw; 
           1 -Rl -1 0;
           1/Rw 1 (-1/Rw-1/Rshunt) -1/Rw;
           1/Rw 0 -1/Rw -1/Rw];
    A3 = A1; 
    A4 = A2;

    B1 = [1/ron 0; 0 0; 0 -1; 0 0];
    B2 = [0 0; 0 0; 0 -1; 0 0];
    B3 = [0 0; 0 0; 0 -1; 0 0];
    B4 = B2;
end


% ronS = 5*ron;
% A1(1,1) = (-1/(ronS) - 1/Rw);
% B1(1,1) = 1/(ronS);

As = cat(3, A1, A2, A3, A4);
Bs = cat(3,B1, B2, B3, B4);

 
HS_MOS = MOSFET(ron, qg, Rja, CossVds, IsdVsd, RonT, Vbr);
LS_MOS = HS_MOS;
Buck = BuckConverter(V, Cout, L, Rl, Rl_400k, Rl_800k, fs, Ts, dt, Vdr, HS_MOS, LS_MOS, As, Bs, Cw);
[Xss, ys, t, ts] = simulate(Buck, Io, Vg);

TestPd = 2.5002e6:2.5002e6+1740;
TestTime = linspace(0,Ts,length(TestPd));
ILskew = -30;

subplot(2,1,1)
hold off;
plot(TestTime,V_DS(TestPd), 'linewidth', 3);
hold on;
plot(TestTime,V_OUT(TestPd), 'linewidth', 3);

plot(t, ys(1,:), 'linewidth', 5);
plot(t, ys(3,:), 'linewidth', 5);

ylims = ylim;
for tx = cumsum(ts);
    plot([tx tx], [ylims(1) ylims(2)], ':r');
end


subplot(2,1,2)
hold off;
plot(TestTime,smooth(I_L(TestPd + ILskew),100), 'linewidth', 3);
hold on;
if(useCw)
    plot(t, ys(2,:) + (ys(1,:)-ys(3,:)-ys(4,:))/Rw, 'linewidth', 5);
else
    plot(t, ys(2,:), 'linewidth', 5);
end

ylim([.75 2.5])

ylims = ylim;
for tx = cumsum(ts);
    plot([tx tx], [ylims(1) ylims(2)], ':r');
end

