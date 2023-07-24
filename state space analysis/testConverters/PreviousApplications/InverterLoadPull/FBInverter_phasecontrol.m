clear all;
clc;
close all;

debug = 0; % turn on to debug convergence of switch time loop
plotWFs = 0; %Turn on to plot time-domain waveforms (slows speed significantly)

%% Phase shift and dead time to simulate
dt_pct = .2;
phi = pi/2;

%% NV6131 GaNFET
Vbr = 650;
ron = .135;
qg = 2.5e-9;
Rja = 50;
NV6131;
FB_FET = MOSFET(ron, qg, Rja, CossVds, 0, 0, Vbr);

%% Circuit Design
Vdc = 200;
P = 100;
zo = 170;

Cp = FB_FET.CeqQ(find(FB_FET.Vds > Vdc,1));
Conom = 28e-12;
Lonom = 19.68e-6;

% Lr = 2.6e-6;
Lr = 1.3e-6;
Cr = 2.2e-6;

Ron = FB_FET.ron;

Rshunt = 100e3;


%% Timing parameters
fs = 6.78e6;
Ts = 1/fs;

tphi = phi*Ts/2/pi;

dtmax = 2*pi*sqrt(Lr*Cp)/4;

ton = Ts/2*(1-dt_pct);
dt = Ts/2*dt_pct;

% t1= ton/4; t2 = dt/2; t3 = 3*ton/4; t4=dt/2;
t1 = tphi; t2 = dt; t3 = (Ts/2-2*dt-tphi); t4 = dt;
ts = [t1 t2 t3 t4];

%% Get impedances to span smith chart
% generate polar coordinates
r = .1:.1:.9;
theta = [0:pi/16:31*pi/16];%[-pi/4:pi/16:pi/4 pi/4:pi/8:7*pi/4];
[radrange, angrange] = meshgrid(r, theta);
radrange = radrange(:);
angrange = angrange(:);
G = radrange.*exp(1i*angrange);  %Gammas spanning smith chart

%Get back impedances
Rrange = [zo; real(gamma2z(G,zo))];
Zreacrange = [0; imag(gamma2z(G,zo))];

% grab default colors for plotting later
h = figure(1);
colors = get(gca,'ColorOrder');
close(h);

if(~debug)
    h = waitbar(0,'Initializing');
end

errits = [];

for it = 1:1:length(Zreacrange)
    tic
    Zreac = Zreacrange(it);
    
    % Load is series L-R-C, so reactive components are added in series to
    % the ouutput series matching so that states remain independent
    if(Zreac > 0)
        Lo = Lonom + Zreac/2/pi/fs;
        Co = Conom;
    elseif(Zreac < 0)
        Lo = Lonom;
        Co = 1/(1/Conom + 1/abs(1/Zreac/2/pi/fs));
    else
        Lo = Lonom;
        Co = Conom;
    end
    Ro = sqrt(Lonom/Conom);
    
    R = Rrange(it);
    
%     Vpk = Vdc*4/pi;
%     Ipk = Vpk/(R + abs(1i*2*pi*fs*Lo - 1i/(2*pi*fs*Co)));

    %% Subinterval state matrices
    % x = [vp ilr vcr ilo vco]
    % u = [Vdc]

    K = [Cp 0 0 0 0; 0 Lr 0 0 0; 0 0 Cr 0 0; 0 0 0 Lo 0; 0 0 0 0 Co];

    %% Subinterval 1, Q1 and Q4 on
    A1 = [-1/2/Ron      0       0       0       0;
            1           0       -1      0       0;
            0           1       0       0       0;
            1           0       0       -R      -1;
            0           0       0       1       0];

    B1 = [1/2/Ron;       0;       0;       0;       0];

    %% Subinterval 2, All FETs OFF
    A2 = [  -1/Rshunt   -1      0       -1      0;
            1           0       -1      0       0;
            0           1       0       0       0;
            1           0       0       -R      -1;
            0           0       0       1       0];

    B2 = [0;       0;      0;       0;       0];

    
    %% Subinterval 3a, zero state
    A3 = [-1/2/Ron      0       0       0       0;
            1           0       -1      0       0;
            0           1       0       0       0;
            1           0       0       -R      -1;
            0           0       0       1       0];

    B3 = [0;       0;       0;       0;       0];
    
    IHC = -eye(5);

    As = cat(3, A1, A2, A3, A2);
    Bs = cat(3, B1, B2, B3, B2);
    ts = [t1 t2 t3 t4];
    u = Vdc;

    for i = 1:size(As,3)
        As(:,:,i) = K^-1*As(:,:,i);
        Bs(:,:,i) = K^-1*Bs(:,:,i);
    end

    [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
    
    if(~debug)    
        h = waitbar(it/length(Zreacrange),h,['Measuring datapoint ' num2str(it) ' of ' num2str(length(Zreacrange))]);
    end
    
    ts0 = ts;

    % check if diodes conducted or if hard switching occured that was
    % avoidable
    diodcon1 = Xss(1,1) > Vdc+2 || Xss(1,1) < -2;
    diodcon3 = Xss(1,3) > Vdc+2 || Xss(1,3) < -2;
    hardswitching1 = (Xss(1,1) < Vdc-2 && (ts(4) < ts0(4))) && (Xss(2,1)+Xss(4,1)) < 0;
    hardswitching3 = (Xss(1,3) > 2 && (ts(2) < ts0(2)))&& (Xss(2,3)+Xss(4,3)) > 0;


    %% dead time iteration loop
%     try
        while((diodcon1 || diodcon3)  || ((hardswitching1 || hardswitching3) ));
            if diodcon1 || (hardswitching1 && ~diodcon3)
                dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 4, 1e-9, 1, IHC);
                dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 4, -1e-9, 1, IHC);

                dXds = (dXsp - dXsn)/2e-9;
                deltaT = (Xss(1,1)-Vdc)/dXds(1,1);

                ttest = linspace(0,ts(4),100);
                [y, tx, x] = lsim(ss(As(:,:,2), Bs(:,:,2), ones(1,length(As)), 0) , u*ones(size(ttest)), ttest, Xss(:,4));
                if(sum(x(:,1)<-Vdc)>3)
                     deltaT = 10e-9;
                elseif(sum(x(:,1)>0) >length(x(:,1))-3)
                    deltaT = ts(4);
                end
                    


                ts(4) = max(1e-12,min(ts(4) - deltaT, ts0(4)));
                ts(1) = sum(ts0([1 4]))-ts(4);
                [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
            end
            
            if (diodcon3 || (hardswitching3)) && ~(diodcon1)
                dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 2, 1e-9, 3, IHC);
                dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 2, -1e-9, 3, IHC);

                dXds = (dXsp - dXsn)/2e-9;
                deltaT = (Xss(1,3))/dXds(1,3);
%                     deltaT = min(ts0(1)-ts(1), deltaT);  % Cannot make on-time larger

                if(diodcon3 && dXds(1,3) > 0)
                    deltaT = 10e-9;
                end

                ts(2) = max(1e-12,min(ts(2) - deltaT, ts0(2)));
                ts(3) = sum(ts0([2 3]))-ts(2);
                [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
            end

            [ ys, t ] = SS_WF_Reconstruct( Xss(:,2), As(:,:,2), Bs(:,:,2), ts(2), u );

            diodcon1 = Xss(1,1) > Vdc+2 || Xss(1,1) < -2;
            diodcon3 = Xss(1,3) > Vdc+2 || Xss(1,3) < -2;
            hardswitching1 = (Xss(1,1) < Vdc-2 && (ts(4) < ts0(4))) && (Xss(2,1)+Xss(4,1)) < 0;
            hardswitching3 = (Xss(1,3) > 2 && (ts(2) < ts0(2)))&& (Xss(2,3)+Xss(4,3)) > 0;
            
            if(debug)
                hl = findobj('type','line');
                delete(hl)
                [ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u );
                ys = ys(1:size(ys,1),1:length(t));
                ts
                it
                [diodcon1  diodcon3 hardswitching1 hardswitching3]
                plotInvWFs( t, ys, R, colors, Ts);
                breakpoint_x=1;
            end

            timer = toc;
            if(~debug && ((~plotWFs && timer > 1) || (plotWFs && timer > 5)))
                display(['Timed out solving diode behavior at it = ' num2str(it)])
                errits = [errits it];
                break
            end

        end

[ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u );
ys = ys(1:size(ys,1),1:length(t));

if(mean(ys(4,:).^2*R) < 150)
    if(plotWFs)
        plotInvWFs( t, ys, R, colors, Ts);
    end
else
    display(['Output power of Po = ' num2str(mean(ys(4,:).^2*R)) 'at it = ' num2str(it)])
end

%% store values
Pout(it) = mean(ys(4,:).^2*R);
Rout(it) = R;
Zout(it) = Zreac;
Lout(it) = Lo;
Cout(it) = Co;

%     catch ME
%         % on error, continue to other points
%         Pout(it) = -1;
%         Rout(it) = R;
%         Zout(it) = Zreac;
%         Lout(it) = Lo;
%         Cout(it) = Co;
%         display('Error Occured');
%     end

end

close(h)


%% Plot smith chart
addpath('Songnan Code')
figure(3)

Rho=((Rout+1i*Zout)-zo)./((Rout+1i*Zout)+zo);
Greal=real(Rho);
Gimag=imag(Rho);

Greal = Greal(1:length(Pout));
Gimag = Gimag(1:length(Pout));
discard = isnan(Greal);
Greal = Greal(~discard);
Gimag = Gimag(~discard);

[X,Y] = meshgrid(linspace(min(Greal),max(Greal)),linspace(min(Gimag),max(Gimag)));
F = scatteredInterpolant(Greal',Gimag',Pout');
Z=F(X,Y);
h = polar([-pi pi], [0 1]);
delete(h);
hold on;
plotSmithChart;
[C,h] = contour(X,Y,Z,30,'LineWidth',4);
h.LevelList = 10:10:100;
scatter(Greal(Pout<110), Gimag(Pout<110), [], Pout(Pout<110),'filled')
scatter(Greal(errits), Gimag(errits), [],'MarkerFaceColor',[1 0 0])
caxis([0, 100])
colorbar

figure(5)
h = polar([-pi pi], [0 1]);
delete(h);
hold on;
plotSmithChart;
itt = 1:1:length(Zreacrange);
scatter(Greal, Gimag, [], itt,'filled')
scatter(Greal(errits), Gimag(errits), [],'MarkerFaceColor',[1 0 0])
hold on;
text(Greal, Gimag, cellstr(num2str(itt')))
title('Locations of it')


