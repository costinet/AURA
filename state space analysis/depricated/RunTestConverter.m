clear all; 
clc;

WFanalysis = 1;
plotIterWF = 1;

testDir = [erase(mfilename('fullpath'), mfilename) 'testConverters'];

addpath(testDir);
converter = 'Buck';
    %possible options
        %HybridDickson
        %Buck
        %Buck_nonSingular

        
simulator = SMPSim();
simulator.loadTestConverter(converter);

if strcmp(converter, 'HybridDickson')
    Ioloc = 2;
    Vgloc = 1;
    Iorange = logspace(log10(.5), log10(18), 25);
    Voutloc = 4;
elseif strcmp(converter, 'Buck')
    Ioloc = 2;
    Vgloc = 1;
    Iorange = logspace(log10(.5), log10(18), 25);
    Voutloc = 4;
elseif strcmp(converter, 'Buck_nonSingular')
    Ioloc = 2;
    Vgloc = 1;
    Iorange = logspace(log10(.5), log10(5), 25);
    Voutloc = 3;
end

zeroAlloc = zeros(size(Iorange));
Pout = zeroAlloc;
Pin = zeroAlloc;
eta = zeroAlloc;
PoutWF = zeroAlloc;
PinWF = zeroAlloc;
etaWF = zeroAlloc;

tic
for j = 1:length(Iorange)
    simulator.u(Ioloc) = Iorange(j);
    Xss = simulator.SS_Soln();
    simulator.regulate();
    
     if(WFanalysis)
        [xs, t, y] = simulator.SS_WF_Reconstruct();
        if(plotIterWF)
            figure(1)
            ns = size(xs,1);
            for i=1:ns
                subplot(10*ns,1,i*10-9:i*10)
                hold on;
                plot(t,xs(i,:), 'Linewidth', 3); 
                ylabel(simulator.getstatenames{i})
                box on

                if(i<ns)
                    set(gca, 'Xticklabel', []);
                else
                    xlabel('t')
                end
            end
            drawnow;
        end
        PoutWF(j) = trapz(t,xs(Voutloc,:).*Iorange(j))/max(t);
        PinWF(j) = trapz(t,y*simulator.u(Vgloc))/max(t);
        etaWF(j) = PoutWF(j)/PinWF(j);
     end
     
     [ avgXs, avgYs ] = simulator.ssAvgs(Xss);
     Pout(j) = avgXs(Voutloc)*Iorange(j);
     Pin(j) = avgYs*simulator.u(Vgloc);
     eta(j) = Pout(j)/Pin(j);

end
toc

figure(2);
plot(Pout, eta, 'Linewidth', 3);
if(WFanalysis)
    hold on;
    plot(Pout, etaWF, ':r', 'Linewidth', 3); 
end

 

    
    
    
    