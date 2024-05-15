%% General Demos
[ts1,swvec1] = dutyMod(.5,1e-6, 'dt', 100e-9, 'phase', -.25e-6, 'phaseUnits', 'time');
plotModWf(ts1,swvec1,1);

[ts2,swvec2] = dutyMod(.5,1e-6, 'dt', 300e-9, 'phase', .05e-6, 'phaseUnits', 'time');
plotModWf(ts2,swvec2,2);

switchCol1 = (1:2);
switchCol2 = (3:4);

[ts,swvec] = combineMods(ts1,swvec1,switchCol1,ts2,swvec2,switchCol2);
plotModWf(ts,swvec,3);

[ts3,swvec3] = dutyMod(.1,.25e-6, 'dt', 5e-9, 'phase', -.25e-6, 'phaseUnits', 'time');
[ts4,swvec4] = dutyMod(.9,.75e-6, 'dt', 10e-9, 'phase', -.25e-6, 'phaseUnits', 'time');
[ts5,swvec5] = concatMods(ts3,swvec3,ts4,swvec4);
plotModWf(ts5,swvec5,4);

[ts,swvec] = combineMods(ts,swvec,(1:4),ts5,swvec5,(5:6));
plotModWf(ts,swvec,6);

%% DAB PSM
Ts = 1e-6;
phi = pi/4;
[ts1,swvec1] = dutyMod(.5,Ts, 'dt', 20e-9);
[ts2,swvec2] = dutyMod(.5,Ts, 'dt', 20e-9, 'phase', phi, 'phaseUnits', 'rad');

[ts,swvec] = combineMods(ts1,swvec1,[1 2; 4 3],ts2,swvec2,[5 6; 8 7]);
plotModWf(ts,swvec,8);

%% Alternative with phaseshift
[ts1,swvec1] = dutyMod(.5,Ts, 'dt', 20e-9);
[ts2,swvec2] = phaseShiftMod(ts1,swvec1,phi, 'phaseUnits', 'rad');

[ts,swvec] = combineMods(ts1,swvec1,[1 2; 4 3],ts2,swvec2,[5 6; 8 7]);
hold(gcf().Children,"on")
plotModWf(ts,swvec,8, ':');