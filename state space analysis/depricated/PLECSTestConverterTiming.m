% clear all; 
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
%         Lr = 0;
        Lf = .5e-6;
        Co = 100e-6;
        ron = 2.5e-3;

        Ts = 5e-6;
        % ts = Ts/4*ones(1,4);
        dt = 50/100;
        ts = [1/48, dt, 47/48, dt];
        ts = ts/sum(ts)*Ts;
        

        Vf = 1;

        Cbnd = diag([0, 0, 1, 1]);
        Dbnd = [0, 0, -Vf, -Vf]';
        
        swseq = [0 1; 0 0 ; 1 0];
        
        Iorange = 1;
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
% td = linspace(0, Ts/15,30);
% td = logspace(log10(Ts/500), log10(Ts/20),50);
N = 100;

td = linspace(Ts/1000, Ts/30, N);

ErrD = zeros(N,N);
ErrDfull = zeros(N,N,2);
errH = zeros(N,N);
errL = zeros(N,N);
nI = zeros(N,N);

% convPoint = zeros(N,N,2);

tsmax = [Ts Ts Ts Ts];
for i=1:length(ts)
    maxEig = max(imag(eigs(sim.As(:,:,i))));
    tsmax(i) = min(Ts, 1/2*pi/maxEig);
end
    

for i = 1:length(td)
    for j = 1:length(td)
        ts = [1/48*Ts, td(i), 47/48*Ts - td(i) - td(j) ,td(j)];
% if(1); i=1;
%     if(1); j=1;
%         ts = [1/48*Ts, Ts/15, 47/48*Ts  - Ts/15 ,0];

%         for k = 1:length(top.inputLabels)
%             u(k) = eval(top.inputLabels{k});
%         end
        sim.u = u;
        sim.ts = ts;
        Xss = sim.SS_Soln();
%         sim.regulate();
        
%%
        Err = [Xss(3,3); (Xss(4,1)-(Vg+1))];
        niter = 0;

        display((i-1)*N + j)
%         nI(i,j) = niter;
        
        ErrD(i,j) = max(abs(Err));
        
        if Xss(3,3) < -Vf 
            ErrDfull(i,j,1) = abs(Xss(3,3) + Vf);
        elseif Xss(3,3) > Vg + Vf 
            ErrDfull(i,j,1) = abs(Xss(3,3) - (Vg + Vf));
        else
            ErrDfull(i,j,1) = 0;
        end
        
        if Xss(4,1) < -Vf 
            ErrDfull(i,j,2) = abs(Xss(4,1) + Vf);
        elseif Xss(4,1) > Vg + Vf 
            ErrDfull(i,j,2) = abs(Xss(4,1) - (Vg + Vf));
        else
            ErrDfull(i,j,2) = 0;
        end
        
        

%%
        [xs, t, y] = sim.SS_WF_Reconstruct();
        PoutWF(i,j) = trapz(t,xs(Voloc,:).*Io)/max(t);
        PinWF(i,j) = trapz(t,y(Igloc,:)*sim.u(Vgloc))/max(t);
        etaWF(i,j) = PoutWF(i,j)/PinWF(i,j);

%         sim.plotAllStates(1);

%         err = zeros(4,length(t));
%         for k = 1:length(t)
%             err(:,k) = (Cbnd*xs(:,k) - Dbnd);
%         end
%         
%         errH(i,j) = -min(0,min(err(4,:)));
%         errL(i,j) = -min(0,min(err(3,:))); 

            for k = 1:length(t)
                if xs(3,k) < -Vf 
                    errH(i,j) = max(errH(i,j), abs(xs(3,k) + Vf));
                elseif xs(3,k) > Vg + Vf 
                    errH(i,j) = max(errH(i,j), abs(xs(3,k) - (Vg + Vf)));
                end

                if xs(4,k) < -Vf 
                    errL(i,j) =  max(errL(i,j), abs(xs(4,k) + Vf));
                elseif xs(4,k) > Vg + Vf 
                    errL(i,j) = max(errL(i,j), abs(xs(4,k) - (Vg + Vf)));
                end
            end
%         while(any(abs(Err) > 0.1))
%             [J, J2, XssF, XssB, X0, dt] = sim.discreteJacobian(2);
% 
%             % Need Xss(3,3) > 0, Xss(4,1) < Vg
%             % Modify t2 and t4
% 
%             %Reduce Jacobian to relevant terms
% 
%             %J(state, at time, time changed)
%             Jact = [J(3,3,2), J(4,1,2);
%                     J(3,3,4), J(4,1,4)];
% 
%             Err = [Xss(3,3); (Xss(4,1))];
%             corr = -Jact\Err;
%             corr = [0, corr(1), 0, corr(2)];
%             
%             if any(abs(corr) > Ts)      %% occasionally you hit is just right and the correction is infinity
%                 ind = abs(corr) > Ts;
%                 sgn = sign(corr);
%                 corr(ind) = 0;
%                 corr =  corr + sgn*Ts/100.*(ind);
%             end
% 
%             newts = sim.ts;
% 
%             for ii = 1:length(corr)
%                 while any(newts(ii) + corr(ii) < 0)
%                     corr(ii) = corr(ii)/2;
%                 end
%             end
% 
%             corrShift =  circshift(corr,1);
% 
%             for ii = 1:length(corr)
%                 while any(newts(ii) + corr(ii) - corrShift(ii) < 0)
%                     corrShift(ii) = corrShift(ii)/2;
%                 end
%             end
% 
%             corr = circshift(corrShift,-1);
%             newts = sim.ts + corr - corrShift;
% 
%             sim.ts = newts;
%     %         newts
%             Xss = sim.SS_Soln();
%             
% %             sim.plotAllStates(1);
%             fig = gcf;
% %             hold(fig.Children, 'on')
%             niter = niter+1;
%             
%             if niter > 100
%                 break
%             end
%         end
%         
%         convPoint(i,j,1) = sim.ts(2);
%         convPoint(i,j,2) = sim.ts(4);
        
%         display(niter)


%          [ avgXs, avgYs ] = sim.ssAvgs(Xss);
%          Pout(i,j) = avgXs(Voloc)*Io;
%          Pin(i,j) = avgYs(Igloc)*sim.u(Vgloc);
%          eta(i,j) = Pout(i,j)/Pin(i,j);
    end
end
toc

figure(4)
contourf(log10(abs(errH))); % No problem because IL always pushes the other way @ high current
colorbar;

figure(5)
contourf(log10(abs(errL)));  
colorbar;

figure(6)
surf(td,td,log10(max(abs(errH),abs(errL))));
% hold on;
% surf(td,td,zeros(length(td),length(td)))


figure(7)
surf(td/Ts,td/Ts,log10(ErrD));
hold on;

%% Plots for COMPEL 2020 - DT
AZ = -20.7;
EL = 20.4;

figure(11)
hold off;
s = surf(td/Ts,td/Ts,max(1,(max(ErrDfull,[],3))));
s(1).EdgeColor = [1 1 1]*.2;
s(1).EdgeAlpha = .25;
set(gca, 'Zscale', 'log');
zlim([1 2e3])
xlim([1e-3, 0.0333]);
ylim([1e-3, 0.0333]) 
box on;
xlabel('$t_{d1}/T_s$','Interpreter','latex')
ylabel('$t_{d2}/T_s$','Interpreter','latex')
zlabel('max(Err) (Discrete) [V]')

% Convergence Bound
for i = 1:42
    surfCurve(i,2) = td(i)/Ts;
    surfCurve(i,1) = td(43-i)/Ts;
    surfCurve(i,3) = max(1,(max(ErrDfull(i, 43-i,:),[],3)));
end
hold on;
plot3(surfCurve(:,1),surfCurve(:,2),surfCurve(:,3), 'r', 'LineWidth', 3)

% Contour projection
[X, Y, Z] = meshgrid(td/Ts, td/Ts, 1:100);
V = repmat((max(1,(max(ErrDfull,[],3)))), 1, 1, 100);
cs = contourslice(X,Y,Z,V,[], [], 1, [1 5 10 50 100 200 500 1000]);

pConv = plot3(surfCurve(:,1),surfCurve(:,2),ones(1,length(surfCurve)), 'r', 'LineWidth', 3)
view(AZ,EL);

% eigenvalue boundary
ind = find(td < tsmax(2), 1, 'last');
surfCurve = [];
surfCurve(:,2) = td([1:ind, ind*ones(1,ind)])/Ts;
surfCurve(:,1) = td([ind*ones(1,ind), ind:-1:1])/Ts;
surfCurve(:,3) = [max(1,(max(ErrDfull(1:ind,ind,:),[],3)))', ...
    max(1,(max(ErrDfull(ind, ind:-1:1,:),[],3)))];
plot3(surfCurve(:,1),surfCurve(:,2),surfCurve(:,3), '-.m', 'LineWidth', 3);

pEigs = plot3(surfCurve(:,1),surfCurve(:,2),ones(1,length(surfCurve)), '-.m', 'LineWidth', 3)

legend([pConv pEigs], 'RoC', '$f_N$ Limit','Interpreter','latex')

% set(gcf,'renderer','Painters')
% print('discreteErrBuck','-dpf')

%% Plots for COMPEL 2020 - CT

figure(6)
hold off;
s = surf(td/Ts,td/Ts,max(1,(max(abs(errH),abs(errL)))));
s(1).EdgeColor = [1 1 1]*.2;
s(1).EdgeAlpha = .25;
set(gca, 'Zscale', 'log');
zlim([1 50e3])
xlim([1e-3, 0.0333]);
ylim([1e-3, 0.0333]) 
box on;
xlabel('$t_{d1}/T_s$','Interpreter','latex')
ylabel('$t_{d2}/T_s$','Interpreter','latex')
zlabel('max(Err) (Continuous) [V]')


for i = 1:42
    surfCurve(i,2) = td(i)/Ts;
    surfCurve(i,1) = td(43-i)/Ts;
    surfCurve(i,3) = max(1,(max(errH(i, 43-i),errL(i, 43-i))));
end
hold on;
plot3(surfCurve(:,1),surfCurve(:,2),surfCurve(:,3), 'r', 'LineWidth', 3)

[X, Y, Z] = meshgrid(td/Ts, td/Ts, 1:100);
V = repmat((max(1,(max(errH,errL)))), 1, 1, 100);
cs = contourslice(X,Y,Z,V,[], [], 1, [1 5 10 50 100 200 500 1000]);

pConv = plot3(surfCurve(:,1),surfCurve(:,2),ones(1,length(surfCurve)), 'r', 'LineWidth', 3)
view(AZ,EL);

% eigenvalue boundary
ind = find(td < tsmax(2), 1, 'last');
surfCurve = [];
surfCurve(:,2) = td([1:ind, ind*ones(1,ind)])/Ts;
surfCurve(:,1) = td([ind*ones(1,ind), ind:-1:1])/Ts;
surfCurve(:,3) = [max(1,(max(errH(1:ind,ind),errL(1:ind,ind))))', ...
    max(1,(max(errH(ind, ind:-1:1),errL(ind, ind:-1:1))))];
plot3(surfCurve(:,1),surfCurve(:,2),surfCurve(:,3), '-.m', 'LineWidth', 3);

pEigs = plot3(surfCurve(:,1),surfCurve(:,2),ones(1,length(surfCurve)), '-.m', 'LineWidth', 3)

legend([pConv pEigs], 'RoC', '$f_N$ Limit','Interpreter','latex')

% set(gcf,'renderer','Painters')
% print('discreteErrBuck','-dpf')
% contourf(sqrt(max(ErrDfull,[],3)));
%% ROC (informal -- uses insight from plot)

figure(8)
subplot(2,1,1)
contourf(td/Ts,td/Ts, sqrt(convPoint(:,:,1)/Ts))
subplot(2,1,2)
contourf(td/Ts,td/Ts, sqrt(convPoint(:,:,2)/Ts))

figure(9)
contourf(log10((convPoint(:,:,1)/Ts).^2 + (convPoint(:,:,2)/Ts).^2))

figure(10)
contourf(td/Ts,td/Ts, log10(ErrD))



% [3,43] -- zero error
% [87,43] -- zero error


%% Example Waveforms

% contourf(max(1,log10((max(ErrDfull,[],3)))))
td1 = td(3); td2 = td(5);

ts = [1/48*Ts, td1, 47/48*Ts - td1 - td2 ,td2];
sim.ts = ts;
Xss = sim.SS_Soln();
% sim.plotAllStates(12)
[xs, t, y] = sim.SS_WF_Reconstruct();

ind1 = [find(t>ts(1)-Ts/50 & t<sum(ts(1:2)) + Ts/50)];
ind2 = [find(t>sum(ts(1:3))-Ts/50 & t<sum(ts(1:4))), find(t<Ts/50)];
timeVec2 = t(ind2);
timeVec2(timeVec2 <= Ts/50) = timeVec2(timeVec2 <= Ts/50) +Ts;

figure(15)
subplot(2,2,1)
plot(t(ind1)/Ts, xs(3,ind1), 'LineWidth', 3)
hold on;  box on;
plot(ts(1)/Ts*ones(1,2), [-10 60], '-.k', 'LineWidth', 1)
plot(sum(ts(1:2))/Ts*ones(1,2), [-10 60], '-.k', 'LineWidth', 1)
xls = xlim;
plot(xls, (Vg+Vf)*ones(1,2), '-.r')
plot(xls, (-Vf)*ones(1,2), '-.r')
ylim([-10 60])
plot(sum(ts(1))/Ts, Xss(3,2), 'gd', 'LineWidth',3)
plot(sum(ts(1:2))/Ts, Xss(3,3), 'gd', 'LineWidth',3)

subplot(2,2,2)
plot(timeVec2/Ts, xs(3,ind2), 'LineWidth', 3)
hold on;  box on;
plot(sum(ts(1:3))/Ts*ones(1,2), [-10 60], '-.k', 'LineWidth', 1)
plot(sum(ts(1:4))/Ts*ones(1,2), [-10 60], '-.k', 'LineWidth', 1)
xls = xlim;
plot(xls, (Vg+Vf)*ones(1,2), '-.r')
plot(xls, (-Vf)*ones(1,2), '-.r')
ylim([-10 60])
plot(sum(ts(1:3))/Ts, Xss(3,4), 'gd', 'LineWidth',3)
plot(sum(ts(1:4))/Ts, Xss(3,5), 'gd', 'LineWidth',3)

td1 = td(54); td2 = td(45);
ts = [1/48*Ts, td1, 47/48*Ts - td1 - td2 ,td2];
sim.ts = ts;
Xss = sim.SS_Soln();
[xs, t, y] = sim.SS_WF_Reconstruct();

ind1 = [find(t>ts(1)-Ts/50 & t<sum(ts(1:2)) + Ts/50)];
ind2 = [find(t>sum(ts(1:3))-Ts/50 & t<sum(ts(1:4))), find(t<Ts/50)];
timeVec2 = t(ind2);
timeVec2(timeVec2 <= Ts/50) = timeVec2(timeVec2 <= Ts/50) +Ts;

figure(15)
subplot(2,2,3); hold on; box on;
plot(t(ind1)/Ts, xs(3,ind1), 'LineWidth', 3)
hold on;
plot(ts(1)/Ts*ones(1,2), [-150 60], '-.k', 'LineWidth', 1)
plot(sum(ts(1:2))/Ts*ones(1,2), [-150 60], '-.k', 'LineWidth', 1)
xls = xlim;
plot(xls, (Vg+Vf)*ones(1,2), '-.r')
plot(xls, (-Vf)*ones(1,2), '-.r')
ylim([-150 60])
plot(sum(ts(1))/Ts, Xss(3,2), 'gd', 'LineWidth',3)
plot(sum(ts(1:2))/Ts, Xss(3,3), 'gd', 'LineWidth',3)

subplot(2,2,4); hold on; box on;
plot(timeVec2/Ts, xs(3,ind2), 'LineWidth', 3)
hold on;
plot(sum(ts(1:3))/Ts*ones(1,2), [-150 60], '-.k', 'LineWidth', 1)
plot(sum(ts(1:4))/Ts*ones(1,2), [-150 60], '-.k', 'LineWidth', 1)
xls = xlim;
plot(xls, (Vg+Vf)*ones(1,2), '-.r')
plot(xls, (-Vf)*ones(1,2), '-.r')
ylim([-150 60])
plot(sum(ts(1:3))/Ts, Xss(3,4), 'gd', 'LineWidth',3)
plot(sum(ts(1:4))/Ts, Xss(3,5), 'gd', 'LineWidth',3)

ErrDfull(54, 45, 1:2)
% max(ErrDfull(43,3,:),[],3)




% return
% figure(2);
% plot(Pout, eta, 'Linewidth', 3);
% if(1)
%     hold on;
%     plot(Pout, etaWF, ':r', 'Linewidth', 3); 
% end

 

    
    
    
    