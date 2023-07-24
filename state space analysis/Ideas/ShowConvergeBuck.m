
sdir = mfilename('fullpath');
sdir = sdir(1:find(sdir=='\',1,'last')-1);
addpath(sdir);

modelfile = 'MRBuckDiodes'; PLECsModel = 'MRBuck';

        
% find_system(modelfile,'SearchDepth',1, 'IncludeCommented', 'on')
open_system(modelfile,'loadonly');
circuitPath = [modelfile '/' PLECsModel];
set_param(circuitPath,'Commented','on');
clear sim;
simout = sim(modelfile,eps);

for i = 1:length(simout.properties)
    assignin('base',simout.properties{i},eval(['simout.' simout.properties{i}]));
end

set_param(circuitPath,'Commented','off');

Lr = 0;
Lf = 0.25e-6;
swvec = [0 0 0 0; 0 1 0 0; 0 0 0 0; 1 0 0 0];
Ts = .2e-6;
ts = [.35 1.2 .25 .14];
ts = ts/sum(ts)*Ts;
% swvec = [0 1 ;  1 0];
% Ts = [5e-6 .1e-6 1e-6];
% Ts = Ts(startPoint);
% ts = [ .9  .1]*Ts;

us(1)=0;
Rload = 1;

sim = SMPSim();

Coss = 1e-9;
CdsL = Coss;
CdsH = Coss;
conv = sim.converter;
top = sim.topology;

top.loadCircuit(circuitPath,swvec,1);
sim.u = us';

conv.setSwitchingPattern(swvec, ts)

sim.steadyState;
sim.plotAllStates(1);

hold(gcf().Children,'on')

us(2)=48;

% ts = [.175 1.02 .075 .1];
ts = ts/sum(ts)*Ts;
top.loadCircuit(circuitPath,swvec,1);
sim.u = us';
conv.setSwitchingPattern(swvec, ts)

% sim.converter.timingThreshold = sim.converter.timingThreshold*1000;

tic
sim.findValidSteadyState;
toc 

sim.plotAllStates(1);

ylim(gcf().Children(end), [2.5 5.5])
ylim(gcf().Children(end-1), [.5 7.2])
ylim(gcf().Children(end-2), [-40 50])

fig = gcf;
fig.Units = 'Inches';
fig.Position = [0 .4 3.5 5];
fig.Renderer = 'painters';
