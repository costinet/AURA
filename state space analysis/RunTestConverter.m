clear all; 
clc;

plotIterWF = 1;

converter = "HybridDickson";

simulator = SMPSim();
simulator.loadTestConverter(converter);

if converter == "HybridDickson"
    Ioloc = 2;
    Vgloc = 1;
    Iorange = logspace(log10(.5), log10(18), 25);
    Voutloc = 4;
    ylabels = {'V_{C1}', 'V_{C2}', 'V_{C3}', 'V_{out}', 'I_L'};
elseif converter == "Buck"
    Ioloc = 2;
    Vgloc = 1;
    Iorange = logspace(log10(.5), log10(18), 5);
    Voutloc = 4;
    ylabels = {'V_{Chs}', 'V_{Cls}', 'i_L', 'V_{out}'};
end


for j = 1:length(Iorange)
    simulator.u(Ioloc) = Iorange(j);
    Xss = simulator.SS_Soln();
    [xs, t, y] = simulator.SS_WF_Reconstruct( );
    
     if(plotIterWF)
        figure(1)
        ns = size(xs,1);
        for i=1:ns
            subplot(10*ns,1,i*10-9:i*10)
            hold on;
            plot(t,xs(i,:), 'Linewidth', 3); 
            ylabel(ylabels{i})
            box on

            if(i<ns)
                set(gca, 'Xticklabel', []);
            else
                xlabel('t')
            end
        end
        drawnow;
     end
    
     Pout(j) = trapz(t,xs(Voutloc,:).*Iorange(j))/max(t);
     Pin(j) = trapz(t,y*simulator.u(Vgloc))/max(t);
     
     eta(j) = Pout(j)/Pin(j);
end

figure(2);
plot(Iorange, eta);

 

    
    
    
    