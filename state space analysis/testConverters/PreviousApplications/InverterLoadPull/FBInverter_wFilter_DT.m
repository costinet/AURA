clear all;
clc;
close all;

debug = 0;
plotWFs = 1;

Cp = 80e-12;
Conom = 28e-12;
Lonom = 19.68e-6;

% Lr = 2.6e-6;
Lr = 1.3e-6;
Cr = 2.2e-6;

Ron = .135;
R = 300;
Rshunt = 100e3;

fs = 6.78e6;
Ts = 1/fs;

dt_pct = .4;

dtmax = 2*pi*sqrt(Lr*Cp)/4;

ton = Ts/2*(1-dt_pct);
dt = Ts/2*dt_pct;

t1= ton/2; t2 = dt/2; t3 = ton/2; t4=dt/2;

Vdc = 200;
P = 100;
zo = 324/4;

%% Farshid's Filter
C5 = 100e-12;
L6 = 2*571.3e-9;
L3 = 739.5e-9;
C1 = 82e-12;
L4 = 2*1.1579e-6;
C2 = 10e-12/2;
L5 = 119.97e-9;
C3 = 180e-12;

s = tf('s');

Zp = 1/s/C5;
Z1 = s*L6;
Z2 = s*L3 + 1/s/C1;
Z3 = 1/(1/s/L4 + s*C2);
Z4 = s*L5 + 1/s/C3;

GVV = Z2/Z1*(R*Z4/(R+Z4))/(R*Z4/(R+Z4) + Z3 + Z1*Z2/(Z1+Z2));

% Zreacrange = linspace(-.5e3,.5e3,20);
% Rrange = linspace(324, 100e3,20);

% Smith chart points
r = .1:.1:.9;
theta = [-pi/4:pi/16:pi/4 pi/4:pi/8:7*pi/4];
[radrange, angrange] = meshgrid(r, theta);
radrange = radrange(:);
angrange = angrange(:);
G = radrange.*exp(1i*angrange);  %Gammas spanning smith chart

%Get back impedances
Rrange = real(gamma2z(G,zo));
Zreacrange = imag(gamma2z(G,zo));

figure(1)
colors = get(gca,'ColorOrder');

% [Zreacrange, Rrange] = meshgrid(linspace(-100,100,5), logspace(log10(324), log10(10e3),5));
% 
% Zreacrange = Zreacrange(:);
% Rrange = Rrange(:);

if(~debug)
    h = waitbar(0,'Initializing');
end

for it = 1:1:length(Zreacrange)
    Zreac = Zreacrange(it);
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
    
    Vpk = Vdc*4/pi;
    Ipk = Vpk/(R + abs(1i*2*pi*fs*Lo - 1i/(2*pi*fs*Co)));

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
    
    
%     %% Subinterval 3, Q2 and Q3 on
%     A3 = [-1/2/Ron      0       0       0       0;
%             1           0       -1      0       0;
%             0           1       0       0       0;
%             1           0       0       -R      -1;
%             0           0       0       1       0];
% 
%     B3 = [-1/2/Ron;       0;       0;       0;       0];

    IHC = -eye(5);

    As = cat(3, A1, A2, A3, A2);%, A3, A2);
    Bs = cat(3, B1, B2, B3, B2);%, B3, B2);
    ts = [t1 t2 t3 t4];% t3 t4];
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

%             diodcon1 = (Xss(1,1) > Vdc);
    diodcon1 = Xss(1,1) > Vdc+2;
    diodcon3 = Xss(1,3) < -2;
    hardswitching1 = (Xss(1,1) < Vdc-2 && (ts(4) < ts0(4)));
    hardswitching3 = (Xss(1,3) > 2 && (ts(2) < ts0(2)));

%     [ ys, t ] = SS_WF_Reconstruct( Xss(:,2), As(:,:,2), Bs(:,:,2), ts(2), u );
%     diodcon2 = (min(ys(1,:)) < -Vdc);

%     hardswitching = Xss(1,3) > -Vdc;

    exitFlag = 0;
    try
        while((diodcon1 || diodcon3)  || ((hardswitching1 || hardswitching3) && exitFlag <20));% || diodcon2);
            if diodcon1 || (hardswitching1 && exitFlag <20)
                dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 4, 1e-9, 1, IHC);
                dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 4, -1e-9, 1, IHC);

                dXds = (dXsp - dXsn)/2e-9;
                deltaT = (Xss(1,1)-Vdc)/dXds(1,1);
%                     deltaT = min(ts0(1)-ts(1), deltaT);  % Cannot make on-time larger

                if(diodcon1 && dXds(1,1) < 0)
                    deltaT = 10e-9;
                end

                ts(4) = max(1e-9,min(ts(4) - deltaT, ts0(4)));
                ts(1) = Ts/4-ts(4);
                [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
            end
            
            if diodcon3 || (hardswitching3 && exitFlag <20)
                dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 2, 1e-9, 3, IHC);
                dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 2, -1e-9, 3, IHC);

                dXds = (dXsp - dXsn)/2e-9;
                deltaT = (Xss(1,3))/dXds(1,3);
%                     deltaT = min(ts0(1)-ts(1), deltaT);  % Cannot make on-time larger

                if(diodcon3 && dXds(1,3) > 0)
                    deltaT = 10e-9;
                end

                ts(2) = max(1e-9,min(ts(2) - deltaT, ts0(2)));
                ts(3) = Ts/4-ts(2);
                [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
            end

%             if diodcon2 && ~diodcon1
% %                     dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 1, -1e-9, 2, IHC);
% %                     dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 1, 1e-9, 2, IHC);
% % 
% %                     dXds = (dXsp - dXsn)/2e-9;
% %                     deltaT = max(min((Xss(1,5)-Vg)/dXds(1,5), ts(4)),0);
% %                     ts(4) = ts(4) - deltaT;
% %                     ts(1) = ts(1) + deltaT;
% %                       [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
%                 [ ys, t ] = SS_WF_Reconstruct( Xss(:,2), As(:,:,2), Bs(:,:,2), ts(2), u );
%                 if(min(ys(1,:)) < -200)
%                      ts(2) = t(find(  ys(1,:) < -200, 1, 'first'));
%                      ts(1) = Ts/2-ts(2);
%                      [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
%                 end
%             end

%                 if hardswitching && (ts(2) < ts0(2))
%                     dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 1, 1e-9, 2, IHC);
%                     dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 1, -1e-9, 2, IHC);
% 
%                     dXds = (dXsp - dXsn)/2e-9;
%                     deltaT = max(min((Xss(1,3)+Vdc)/dXds(1,1), ts(1)),-ts(2));
% %                     deltaT = min(ts0(1)-ts(1), deltaT);  % Cannot make on-time larger
%                     ts(1) = ts(1) - deltaT;
%                     ts(2) = ts(2) + deltaT;
%                     [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
%                 end

%                 hardswitching = Xss(1,3) > -Vdc
            [ ys, t ] = SS_WF_Reconstruct( Xss(:,2), As(:,:,2), Bs(:,:,2), ts(2), u );
%             diodcon1 = Xss(1,1) > Vdc+2;
%             hardswitching = (Xss(1,1) < Vdc-2 && (ts(2) < ts0(2)));
%             diodcon2 = (min(ys(1,:)) < -Vdc-2);
            exitFlag = exitFlag + 1;
            diodcon1 = Xss(1,1) > Vdc+2;
            diodcon3 = Xss(1,3) < -2;
            hardswitching1 = (Xss(1,1) < Vdc-2 && (ts(4) < ts0(4)));
            hardswitching3 = (Xss(1,3) > 2 && (ts(2) < ts0(2)));
            
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



        end


%                 [ Xss] = SS_Soln( As, Bs, ts, u, IHC);

%                 if Voerror
%                     dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 3, 1e-9, 1, IHC);
%                     dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 3, -1e-9, 1, IHC);
% 
%                     dXds = (dXsp - dXsn)/2e-9;
%                     Voavg = sum(Xss(3,1:end-1).*ts)/obj.Ts;
%                     dVdt = mean(dXds(3,:));
%                     deltaT = (Voavg-obj.V)/dVdt;
%                     ts(3) = ts(3) - deltaT;
%                     ts(1) = ts(1) + deltaT;
%                 end
% 
%                 [ Xss] = SS_Soln( As, Bs, ts, u);


%                 diodcon2 = (Xss(1,1) < -Vdc);
%                 Voerror = abs(5- mean(Xss(3,:))) > .01;
%         end

%             [ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u );

[ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u );
ys = ys(1:size(ys,1),1:length(t));

% if(max(ys(1,:)) > 300)
%     Lo
%     Co
%     R
% end

% figure(1);
% colors = get(gca,'ColorOrder');
if(mean(ys(4,:).^2*R) < 120)
    if(plotWFs)
        plotInvWFs( t, ys, R, colors, Ts);
    end
%     subplot(3,1,1);
%     hold on;
%     plot([t t+Ts/2], [ys(1,:), -ys(1,:)], 'Color', colors(1,:));
%     plot([t t+Ts/2], [ys(3,:), -ys(3,:)], 'Color', colors(2,:));
%     plot([t t+Ts/2], [ys(4,:)*R, -ys(4,:)*R], 'Color', colors(4,:));
%     legend('V_p',  'v_{cr}',  'V_{out}');
% 
%     subplot(3,1,2);
%     hold on;
%     plot([t t+Ts/2], [ys(5,:), -ys(5,:)], 'Color', colors(2,:));
%     legend(  'v_{cr}' );
% 
%     subplot(3,1,3);
%     hold on;
%     plot([t t+Ts/2], [ys(2,:), -ys(2,:)], 'Color', colors(1,:));
%     plot([t t+Ts/2], [ys(4,:), -ys(4,:)], 'Color', colors(2,:));
%     legend('i_{lr}',  'i_{lo}');
else
    display(['Output power of Po = ' num2str(mean(ys(4,:).^2*R)) 'at it = ' num2str(it)])
end

Pout(it) = mean(ys(4,:).^2*R);
Rout(it) = R;
% Rout2(it) = R2;
Zout(it) = Zreac;
Lout(it) = Lo;
Cout(it) = Co;

    catch ME
        Pout(it) = -1;
        Rout(it) = R;
        % Rout2(it) = R2;
        Zout(it) = Zreac;
        Lout(it) = Lo;
        Cout(it) = Co;
        display('Error Occured');
    end

end

close(h)



% figure(2)
% [line, hsm] = smithchart(z2gamma(Rout + 1i*Zout + Ro, zo));
% set(line, 'LineWidth', 3);

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
% plot(Greal, Gimag, 'ok')
scatter(Greal(Pout<110), Gimag(Pout<110), [], Pout(Pout<110),'filled')

caxis([0, 100])
colorbar

[FX,FY] = gradient(Z,X,Y);
% quiver(FX,FY);


