clear all; 
clc;

testDir = [erase(mfilename('fullpath'), mfilename) 'testConverters'];
addpath(testDir);

circuit = 'syncbuck';
circuitPath = ['PLECSConverters/', circuit];

if strcmp(circuit, 'syncbuck')
    if(0) % MRBuck
%         Rin = 0;
%         Io = 15/1.5;
%         Vg = 48;
%         CdsL = 1500e-12 + 550*1500e-12;
%         CdsH = 0;
%         Lr = 25e-9;
%         Lf = 100e-6;
%         Co = 1e-6;
%         ron = 1.2e-3;
% 
%         Ts = 1e-6;
%         % ts = Ts/4*ones(1,4);
%         ts = [1/48,44/48];
%         ts = [ts, 1-sum(ts)]*Ts;
%         u = [Io, Vg]';
% 
%         Vf = 1;
% 
%         Cbnd = diag([0, 0, 1, 1, 0]);
%         Dbnd = [1, 1, Vf, Vf, 1]';
%         
%         swseq = [1 1; 1 0; 0 1; 0 0];
%         
%         Iorange = 1;
    elseif(1) % Normal Buck
        Rin = 1e-6;
        Vg = 48;
        CdsL = 1500e-12;
        CdsH = 550e-12;
%         Lr = 0;
        Lf = 1e-6;
        Co = 1e-6;
        ron = 1.2e-3;

        Ts = .1e-6;
        % ts = Ts/4*ones(1,4);
        dt = 50/100;
        ts = [1/48, dt, 47/48, dt];
        ts = ts/sum(ts)*Ts;
        

        Vf = 1;

        Cbnd = diag([0, 0, 1, 1]);
        Dbnd = [1, 1, Vf, Vf]';
        
        swseq = [0 1; 0 0 ; 1 0];
        
        Iorange = 10;
        Io = Iorange(1);
        u = [Io, Vg]';
    end
end
        
sim = SMPSim();
conv = sim.converter;
top = sim.topology;

top.Cbnd = Cbnd; top.Dbnd = Dbnd;




top.loadPLECsModel(circuitPath,swseq);
conv.ts = ts;
conv.swseq = [1 2 3 2];

Igloc = find(strcmp(top.outputLabels, 'Ig'));
Vgloc = find(strcmp(top.inputLabels, 'Vg'));

Ioloc = find(strcmp(top.inputLabels, 'Io'));
Voloc = find(strcmp(top.outputLabels, 'Vo'));

% if strcmp(circuit, 'HybridDickson')
%     Ioloc = 2;
%     Vgloc = 1;
%     Iorange = logspace(log10(.5), log10(18), 25);
%     Voutloc = 4;
% elseif strcmp(circuit, 'Buck')
%     Ioloc = 2;
%     Vgloc = 1;
%     Iorange = logspace(log10(.5), log10(18), 25);
%     Voutloc = 4;
% elseif strcmp(circuit, 'Buck_nonSingular')
%     Ioloc = 2;
%     Vgloc = 1;
%     Iorange = logspace(log10(.5), log10(5), 25);
%     Voutloc = 3;
% end

zeroAlloc = zeros(size(Iorange));
Pout = zeroAlloc;
Pin = zeroAlloc;
eta = zeroAlloc;
PoutWF = zeroAlloc;
PinWF = zeroAlloc;
etaWF = zeroAlloc;

tic
% for td1 = linspace(0, Ts/20,100);
%     for td2 = linspace(0, Ts/20,100);
for j = 1:length(Iorange)
    Io = Iorange(j);
    for i = 1:length(top.inputLabels)
        u(i) = eval(top.inputLabels{i});
    end
    sim.u = u;
    sim.ts = ts;
    Xss = sim.SS_Soln();
%     sim.regulate();
    

    [xs, t, y] = sim.SS_WF_Reconstruct();
    PoutWF(j) = trapz(t,xs(Voloc,:).*Iorange(j))/max(t);
    PinWF(j) = trapz(t,y(Igloc,:)*sim.u(Vgloc))/max(t);
    etaWF(j) = PoutWF(j)/PinWF(j);
    
    sim.plotAllStates(1);
    
    err = (Cbnd*xs(:,end/2) + Dbnd);
    errH = err

     
     [ avgXs, avgYs ] = sim.ssAvgs(Xss);
     Pout(j) = avgXs(Voloc)*Iorange(j);
     Pin(j) = avgYs(Igloc)*sim.u(Vgloc);
     eta(j) = Pout(j)/Pin(j);

end
toc

return
figure(2);
plot(Pout, eta, 'Linewidth', 3);
if(1)
    hold on;
    plot(Pout, etaWF, ':r', 'Linewidth', 3); 
end

 

    
    
    
    