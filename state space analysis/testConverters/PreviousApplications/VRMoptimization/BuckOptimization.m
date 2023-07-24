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
Io = 1;

Imax = 10;
Iorange = linspace(0,Imax,100);
EstEfficiencyIo = Iorange./(Iorange + Imax/100 + Iorange.^2/Imax/8);
maxlocIo = find(EstEfficiencyIo == max(EstEfficiencyIo));
% plot(Iorange, EstEfficiencyIo)

fsmax = 4e6;
fsrange = logspace(5, log10(fsmax), 100);
EstEfficiencyfs = 1./(1+1e-8*fsrange + 1./fsrange.^2*max(fsrange)*100);
maxlocF = find(EstEfficiencyfs == max(EstEfficiencyfs));
% plot(fsrange, EstEfficiencyfs);

i=1;
etas = zeros(1,20);
Ios =  zeros(1,20);
Iosteps =  zeros(1,20);
fss = zeros(1,20);
fssteps = zeros(1,20);
exitflags =  zeros(4,3,20);

while(1)
    [eta, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io, Cout);
    
    dI = max(Io/100, .01);
    [eta_dI, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io+dI, Cout);
    dEtadI = (eta_dI-eta)/dI;
    [eta_dI2, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io+2*dI, Cout);
    dEta2d2I = ((eta_dI2-eta_dI)/dI - (eta_dI-eta)/dI)/dI;
    

    
    df = 10e3;
    [eta_df, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs+df, Vg, V, Io, Cout);
    dEtadF = (eta_df-eta)/df;
    [eta_df2, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fs+2*df, Vg, V, Io, Cout);
    dEta2d2F = ((eta_df2-eta_df)/dI - (eta_df-eta)/dI)/dI;
    
    H = [dEta2d2I];
    gradEta = [eta_dI dEtadF];
    
    dir = [dEtadI*Imax, dEtadF*100e6];
    
    
    %% find partial location on estimated curves
    if(dEtadI > 0)
        locIo = find(abs(diff(EstEfficiencyIo)./diff(Iorange) - dEtadI) == ...
            min(abs(diff(EstEfficiencyIo)./diff(Iorange) - dEtadI)),1);
        stepmagIo = abs(diff(Iorange([locIo maxlocIo])))*(Imax-Io)/(Imax-Iorange(locIo));
    else
        locIo = find(abs(diff(EstEfficiencyIo)./diff(Iorange) - dEtadI) == ...
            min(abs(diff(EstEfficiencyIo)./diff(Iorange) - dEtadI)),1, 'last');
        stepmagIo = abs(diff(Iorange([locIo maxlocIo])))*(Io/Iorange(locIo));
    end
    
    if(sign(dEtadI) < 0)
        stepmagIo = min(stepmagIo, Io/2);
    else
        stepmagIo = min(stepmagIo, abs(Imax-Io)/2);
    end
    
    if(dEtadF > 0)
        locF = find(abs(diff(EstEfficiencyfs)./diff(fsrange) - dEtadF) == ...
            min(abs(diff(EstEfficiencyfs)./diff(fsrange) - dEtadF)),1);
        stepmagF = abs(diff(fsrange([locF maxlocF])))*(fsmax-fs)/(fsmax-fsrange(locF));
    else
        locF = find(abs(diff(EstEfficiencyfs)./diff(fsrange) - dEtadF) == ...
            min(abs(diff(EstEfficiencyfs)./diff(fsrange) - dEtadF)),1, 'last');
        stepmagF = abs(diff(fsrange([locF maxlocF])))*fs/fsrange(locF);
    end
    
    delta_Io = stepmagIo*sign(dEtadI);
    delta_f = stepmagF*sign(dEtadF);
    
%     mag = (1-eta)*1000;
    
%     delta_f = mag*dir(2);
%     delta_Io = mag*dir(1);
    
    
    etas(i) = eta;
    Ios(i) = Io;
    fss(i) = fs;
    Iosteps(i) = delta_Io;
    fssteps(i) = delta_f;
    
    exitflags(:,:,i) = exitflag.values;
    
%     plot(Io, fs/1e6, 'or', 'linewidth',3); hold on;
    scatter(Ios, fss/1e6, [],etas);
    colorbar; caxis([.8 1]);
    xlabel('I_o');
    ylabel('f_s [MHz]');
    hold on; plot(Io, fs/1e6, 'or'); hold off;

    drawnow

    
% %     plot(Io, eta, 'or');
% % %     ylim([.9 1])
% %     
% %     if(abs(dEtadI) < .001/(Imax/100))
% %         dI = (Imax/100);
% %         eta_p = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io+dI, Cout);
% %         eta_m = BuckSimulation( inductor, HSFET, LSFET, fs, Vg, V, Io-dI, Cout);
% %         if(eta_p <= eta && eta_m <= eta)
% %             break
% %         end
% %     end
    
    Io = Io + delta_Io;
    fs = fs + delta_f;
    i = i+1
end

F = scatteredInterpolant(Ios', fss'/1e6, etas');
[X, Y] = meshgrid(Iorange, fsrange/1e6);
contourf(X,Y, F(X,Y), 'Levellist', [.97 .95 .93 .91 .89 .87 .85 .83 .81 .79 0]); caxis([.8 1]);
hold on;
scatter(Ios, fss/1e6, [],etas);
colorbar; caxis([.8 1]);
xlabel('I_o');
ylabel('f_s [MHz]');

% plot(Io, eta, 'xg', 'linewidth',2, 'MarkerSize',10);
