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
        CdsL = 880e-12;
        CdsH = 190e-12;
        
        CdsL = CdsL + CdsH;
        CdsH = 0;
        
%         Lr = 0;
        Lf = 1e-6;
        Co = 1e-6;
        ron = 2.5e-3;

        Ts = .1e-6;
        % ts = Ts/4*ones(1,4);
        dt = 50/100;
        ts = [1/48, dt, 47/48, dt];
        ts = ts/sum(ts)*Ts;
        

        Vf = 1;

        Cbnd = diag([0, 0, 1, 1]);
        Dbnd = [0, 0, Vf, Vf]';
        
        swseq = [0 1; 0 0 ; 1 0];
        
        Iorange = 10;
        Io = Iorange(1);
        u = [Io, Vg]';
    end
end
        
sim = SMPSim();
conv = sim.converter;
top = sim.topology;

% top.Cbnd = Cbnd; top.Dbnd = Dbnd;




top.loadPLECsModel(circuitPath,swseq);
conv.ts = ts;
conv.swseq = [1 2 3 2];

Igloc = find(strcmp(top.outputLabels, 'Ig'));
Vgloc = find(strcmp(top.inputLabels, 'Vg'));

Ioloc = find(strcmp(top.inputLabels, 'Io'));
Voloc = find(strcmp(top.outputLabels, 'Vo'));


Votest = linspace(0, 5, 39);
Iltest = linspace(5, 20, 30);
Vdstest = linspace(-1, 49, 31);

t = linspace(0,Ts,10000);
tdmax = 1/15*Ts;

for i = 1:length(Votest)
    for j = 1:length(Iltest)
        for k = 1:length(Vdstest)
            x0 = [Votest(i), Iltest(j), Vdstest(k)]';
   
            for n = 1:4
                sys = ss(sim.As(:,:,n), sim.Bs(:,:,n), sim.Cs(:,:,n), sim.Ds(:,:,n));
                [Y,T,X] = lsim(sys,repmat(u,1,length(t)),t,x0);
                
                if n == 1
                    ind = find(t>1/48*Ts,1,'first');
                    xs = X(1:ind,:)';
                    ys = Y(1:ind,:)';       
                elseif n == 2
                    ind1 = find(X(:,3) > Vg+Vf,1,'first');
                    ind2 = find(X(:,3) < -Vf,1,'first');
                    ind3 = find(t>tdmax,1,'first');
                    
                    ind = min([ind1,ind2,ind3]);
                    xs = [xs, X(1:ind,:)'];
                    ys = [ys, Y(1:ind,:)'];  
                    
                    tdt1act = t(ind);                   
                elseif n == 3
                    tmax = 47/48*Ts - tdt1act;
                    ind = find(t>tmax,1,'first');
                    xs = [xs, X(1:ind,:)'];
                    ys = [ys, Y(1:ind,:)'];  
                    
                elseif n == 4
                    ind1 = find(X(:,3) > Vg+Vf,1,'first');
                    ind2 = find(X(:,3) < -Vf,1,'first');
                    ind3 = find(t>tdmax,1,'first');
                    
                    ind = min([ind1,ind2,ind3]);
                    xs = [xs, X(1:ind,:)'];
                    ys = [ys, Y(1:ind,:)'];  
                end
                x0 = xs(:,end);
            end
            
            err(i,j,k) = norm(xs(:,1) - xs(:,end),1);
        end
    end
end
                
            
figure(1)
slice(Iltest, Votest, Vdstest, err, [10 15], [.5 1 2], [-1 24 49])
ylabel('Vo')
xlabel('IL')
zlabel('Vds')
colorbar
set(gca, 'CLim', [0 2])

[xi, yi, zi] = find(err == min(err,[],1:3));
hold on;
plot3(Iltest(yi), Votest(xi), Vdstest(zi), 'or')

figure(2)
isosurface(Iltest, Votest, Vdstest, err, .5);
hold on
isosurface(Iltest, Votest, Vdstest, err, 1);
isosurface(Iltest, Votest, Vdstest, err, 2);
isosurface(Iltest, Votest, Vdstest, err, 3);
plot3(Iltest(yi), Votest(xi), Vdstest(zi), 'or')

% 
% tic
% % td = linspace(0, Ts/15,30);
% td = logspace(log10(Ts/10e3), log10(Ts/15),30);
% for i = 1:length(td)
%     for j = 1:length(td)
%         ts = [1/48*Ts, td(i), 47/48*Ts - td(i) - td(j) ,td(j)];
% % if(1); i=1;
% %     if(1); j=1;
% %         ts = [1/48*Ts, Ts/15, 47/48*Ts  - Ts/15 ,0];
% 
% %         for k = 1:length(top.inputLabels)
% %             u(k) = eval(top.inputLabels{k});
% %         end
%         sim.u = u;
%         sim.ts = ts;
%         Xss = sim.SS_Soln();
%     %     sim.regulate();
% 
% 
%         [xs, t, y] = sim.SS_WF_Reconstruct();
%         PoutWF(i,j) = trapz(t,xs(Voloc,:).*Io)/max(t);
%         PinWF(i,j) = trapz(t,y(Igloc,:)*sim.u(Vgloc))/max(t);
%         etaWF(i,j) = PoutWF(i,j)/PinWF(i,j);
% 
% %         sim.plotAllStates(1);
% 
%         for k = 1:length(t)
%             err(:,k) = (Cbnd*xs(:,k) + Dbnd);
%         end
%         
%         errH(i,j) = min(err(4,:));
%         errL(i,j) = min(err(3,:));
% 
% 
% %          [ avgXs, avgYs ] = sim.ssAvgs(Xss);
% %          Pout(i,j) = avgXs(Voloc)*Io;
% %          Pin(i,j) = avgYs(Igloc)*sim.u(Vgloc);
% %          eta(i,j) = Pout(i,j)/Pin(i,j);
%     end
% end
% toc
% 
% figure(1)
% contourf(errH); % No problem because IL always pushes the other way @ high current
% colorbar;
% 
% figure(2)
% contourf(errL);  
% colorbar;
% 
% figure(3)
% surf(td,td,min(errH,errL));
% hold on;
% surf(td,td,zeros(length(td),length(td)))
% 
% % return
% % figure(2);
% % plot(Pout, eta, 'Linewidth', 3);
% % if(1)
% %     hold on;
% %     plot(Pout, etaWF, ':r', 'Linewidth', 3); 
% % end

 

    
    
    
    