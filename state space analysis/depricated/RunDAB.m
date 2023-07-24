ronp =1e-3;
rons = 1e-3;
Rp = 1e-3;
Rs = 1e-3;
Lp = 10e-6;
np = 48;
ns = 1;
Lm = 100e-6;

Cp = 50e-12;
Cs = 50e-12;

Vg = 48;
Vout = 1;

fs = 1e6;
Ts = 1/fs;

sim = SMPSim();

conv = sim.converter;
top = conv.topology;

addpath('TestConverters');
fn = 'PLECSConverters/DAB';

swseq = [1 0 0 1 1 0 0 1;
    0 0 0 0 1 0 0 1;
    0 1 1 0 1 0 0 1;
    0 1 1 0 0 0 0 0;
    0 1 1 0 0 1 1 0;
    0 0 0 0 0 1 1 0;
    1 0 0 1 0 1 1 0;
    1 0 0 1 0 0 0 0];


top.loadPLECsModel(fn, swseq)
sim.ts = [Ts/2*.9, Ts/2*.01, Ts*.3, Ts/2*.0001, Ts/2*.9, Ts/2*.01, Ts*.1, Ts/2*.0001];
sim.u = [Vg, Vout]';

K = diag([Cp, Cs, Lm, Lp]);
top.K = K;

Xs = sim.SS_Soln();
sim.plotAllStates(1);

