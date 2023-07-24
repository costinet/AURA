clear all;
clc;

debug = 0;

HSFET.ron = .05;
HSFET.Cp = 3.4874e-10;
LSFET = HSFET;

inductor.Rl = .01;
inductor.L = 230e-9;

Vg = 5;
Cout = 100e-6;
fs = 4e6;
V = 1.8;
Io = 5;

Imax = 5;
Iorange = linspace(0,Imax,100);

EstEfficiency = Iorange./(Iorange + Imax/100 + Iorange.^2/Imax/8);
maxloc = find(EstEfficiency == max(EstEfficiency));
% plot(Iorange, EstEfficiency)

load('4MHzSweep.mat');
figure(1); hold on;
plot(Iorange, etas);

anseff = max(etas);
ansIo = Iorange(etas == anseff);
plot(ansIo, anseff, 'db', 'linewidth',3)

i=1;

etas = zeros(1,20);
Ios =  zeros(1,20);
Iosteps =  zeros(1,20);

exitflags =  zeros(4,3,20);

while(1)
    [eta, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io, Cout);
    
    dI = max(Io/100, .01);
    [eta_dI, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io+dI, Cout);
    
    dEtadI = (eta_dI-eta)/dI;
    
    if(dEtadI > 0)
        loc = find(abs(diff(EstEfficiency)./diff(Iorange) - dEtadI) == ...
            min(abs(diff(EstEfficiency)./diff(Iorange) - dEtadI)),1);
    else
        loc = find(abs(diff(EstEfficiency)./diff(Iorange) - dEtadI) == ...
            min(abs(diff(EstEfficiency)./diff(Iorange) - dEtadI)),1, 'last');
    end
 
    stepmag = abs(diff(Iorange([loc maxloc])));
    if stepmag == 0
        stepmag = abs(dEtadI*dI/Io*1 + sign(dEtadI)*abs(1-eta)^3*dI*10000);
    end
    
    if(sign(dEtadI) < 0)
        stepmag = min(stepmag, Io/2);
    else
        stepmag = min(stepmag, abs(Imax-Io)/2);
    end
    
    delta_Io = stepmag*sign(dEtadI);
    
%     delta_Io = dEtadI*dI/Io*1 + sign(dEtadI)*abs(1-eta)^3*dI*10000;
    

    %% oscillation detection
    if(i>10 && sum(diff(diff(sign(Iosteps(end-11:end))))) == 0)
        delta_Io = mean(Iosteps(end-10:end)) - Io;
    end
    
    
    etas(i) = eta;
    Ios(i) = Io;
    Iosteps(i) = delta_Io;
    
    exitflags(:,:,i) = exitflag.values;


    
    plot(Io, eta, 'or');
    drawnow
%     ylim([.9 1])
    
    if(abs(dEtadI) < .001/(Imax/100))
        dI = (Imax/100);
        eta_p = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io+dI, Cout);
        eta_m = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io-dI, Cout);
        if(eta_p <= eta && eta_m <= eta)
            break
        end
    end
    
    Io = Io + delta_Io;
    i = i+1
end

plot(Io, eta, 'xg', 'linewidth',2, 'MarkerSize',10);
