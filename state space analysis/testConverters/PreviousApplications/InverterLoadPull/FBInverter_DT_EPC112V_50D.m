clear all;
clc;
close all force;


%% Load in Farshid's Results
addpath('Farshid Results 112V')
parseLTfile

colors = get(gca,'ColorOrder'); % Get Matlab default colors
PlotAll = 0; % Flag to plot all time domain waveforms (slows down run immensely)

% Device parasitic drain-to-source cap
Cp = 120e-12;
Ron = .150;

% Series Resonant L-C
Conom = 23e-12;
Lonom = 23.87e-6;

% ZVS Tanks
Lr = 1e-6;
Cr = 200e-6;


R = 300;
Rshunt = 100e3;

fs = 6.78e6;
Ts = 1/fs;

% gate driver dead time, as percent of period
dt_pct = .4;

% dtmax = 2*pi*sqrt(Lr*Cp)/4;

ton = Ts/2*(1-dt_pct);
dt = Ts/2*dt_pct;

t1= ton; t2 = dt; t3 = ton; t4=dt;

Vdc = 112;
P = 100;
zo = 98;


% Smith chart points
r = .1:.1:.9; 
theta = [-pi/4:pi/16:pi/4 pi/4:pi/8:7*pi/4];
theta = [-3*pi/4:pi/64:3*pi/4 3*pi/4:pi/8:5*pi/4];
[radrange, angrange] = meshgrid(r, theta);
radrange = radrange(:);
angrange = angrange(:);
G = radrange.*exp(1i*angrange);  %Gammas spanning smith chart

%Get back impedances
Rrange = real(gamma2z(G,zo));
Zreacrange = imag(gamma2z(G,zo));

h = waitbar(0,'Initializing');

buff = zeros(1, length(Zreacrange));
Pout = buff;
Rout = buff;
Zout = buff;
Lout = buff;
Cout = buff;

tstart = tic;
for it = 1:1:length(Zreacrange)
    tic
    Zreac = Zreacrange(it);
    %Adjust series resonant L-C values for load reactance
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

    %% Subinterval 3, Q2 and Q3 on
    A3 = [-1/2/Ron      0       0       0       0;
            1           0       -1      0       0;
            0           1       0       0       0;
            1           0       0       -R      -1;
            0           0       0       1       0];

    B3 = [-1/2/Ron;       0;       0;       0;       0];

    IHC = -eye(5);

    As = cat(3, A1, A2);%, A3, A2);
    Bs = cat(3, B1, B2);%, B3, B2);
    ts = [t1 t2];% t3 t4];
    u = Vdc;

    for i = 1:size(As,3)
        As(:,:,i) = K^-1*As(:,:,i);
        Bs(:,:,i) = K^-1*Bs(:,:,i);
    end

    [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
    h = waitbar(it/length(Zreacrange),h,['Measuring datapoint ' num2str(it) ' of ' num2str(length(Zreacrange))]);

    ts0 = ts;

    diodcon1 = Xss(1,1) > Vdc+2;
    hardswitching = (Xss(1,1) < Vdc-2 && (ts(2) < ts0(2)));

    [ ys, t ] = SS_WF_Reconstruct( Xss(:,2), As(:,:,2), Bs(:,:,2), ts(2), u );
    diodcon2 = (min(ys(1,:)) < -Vdc);

    exitFlag = 0;

        while((diodcon1 || diodcon2)  || (hardswitching && exitFlag <20));
            if diodcon1 || (hardswitching && exitFlag <20)
                %% Correct for diode swithing nonlinearity on transition 1
                dXsp  = StateSensitivity( As, Bs, ts, u, 'ts', 1, 1e-9, 2, IHC);
                dXsn  = StateSensitivity( As, Bs, ts, u, 'ts', 1, -1e-9, 2, IHC);

                dXds = (dXsp - dXsn)/2e-9;
                deltaT = max(min((Xss(1,1)-Vdc)/dXds(1,1), ts(1)),-ts(2));

                if(diodcon1 && dXds(1,1) > 0)
                    deltaT = -10e-9;
                end

                ts(2) = min(ts(2) + deltaT, ts0(2));
                ts(1) = Ts/2-ts(2);
                [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
            end

            if diodcon2 && ~diodcon1
                %% Correct for diode switching nonlinearity on transition 2
                [ ys, t ] = SS_WF_Reconstruct( Xss(:,2), As(:,:,2), Bs(:,:,2), ts(2), u );
                if(min(ys(1,:)) < -200)
                     ts(2) = t(find(  ys(1,:) < -200, 1, 'first'));
                     ts(1) = Ts/2-ts(2);
                     [ Xss] = SS_Soln( As, Bs, ts, u, IHC);
                end
            end

            [ ys, t ] = SS_WF_Reconstruct( Xss(:,2), As(:,:,2), Bs(:,:,2), ts(2), u );
            diodcon1 = Xss(1,1) > Vdc+2;
            hardswitching = (Xss(1,1) < Vdc-2 && (ts(2) < ts0(2)));
            diodcon2 = (min(ys(1,:)) < -Vdc-2);
            exitFlag = exitFlag + 1;

            timer = toc;
            if((~PlotAll && timer > 2) || (PlotAll && timer > 5))
                % If code hangs on one point, skip it.
                break
            end
        end

[ ys, t ] = SS_WF_Reconstruct( Xss, As, Bs, ts, u );
ys = ys(1:size(ys,1),1:length(t));


if(mean(ys(4,:).^2*R) < 120 && PlotAll)
    figure(1)
    subplot(3,1,1);
    hold on;
    plot([t t+Ts/2], [ys(1,:), -ys(1,:)], 'Color', colors(1,:));
    plot([t t+Ts/2], [ys(3,:), -ys(3,:)], 'Color', colors(2,:));
    plot([t t+Ts/2], [ys(4,:)*R, -ys(4,:)*R], 'Color', colors(4,:));
    legend('V_p',  'v_{cr}',  'V_{out}');

    subplot(3,1,2);
    hold on;
    plot([t t+Ts/2], [ys(5,:), -ys(5,:)], 'Color', colors(2,:));
    legend(  'v_{cr}' );

    subplot(3,1,3);
    hold on;
    plot([t t+Ts/2], [ys(2,:), -ys(2,:)], 'Color', colors(1,:));
    plot([t t+Ts/2], [ys(4,:), -ys(4,:)], 'Color', colors(2,:));
    legend('i_{lr}',  'i_{lo}');
end

Pout(it) = mean(ys(4,:).^2*R);
Rout(it) = R;
Zout(it) = Zreac;
Lout(it) = Lo;
Cout(it) = Co;

end

toc(tstart)

close(h)


addpath('Songnan Code')
% figure(3)
hold on;

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

Z_LTSpice = LTF(X,Y);

h = polar([-pi pi], [0 1]);
delete(h);
hold on;
plotSmithChart;
% % [C,h] = contour(X,Y,Z,30,'LineWidth',4);
% % h.LevelList = 10:10:100;
% plot(Greal, Gimag, 'ok')
scatter(Greal(Pout<110), Gimag(Pout<110), [], Pout(Pout<110),'filled')

caxis([0, 100])
colorbar

%%
figure(6);
Err = (Z-Z_LTSpice);%./Z_LTSpice*100; 
Err(Z<10) = 0;
Err(Z>100) = 0;
h = polar([-pi pi], [0 1]);
delete(h);
hold on;
plotSmithChart;
% [C,h] = contour(X,Y,Err,5,'LineWidth',4);
% hold on
Xvec = X(Z<100 & Z>5); Yvec = Y(Z<100 & Z>5); Errvec = Err(Z<100& Z>5);
scatter(Xvec,Yvec, [], Errvec(:), 'filled')
% h.LevelList = 0:.25:10;
caxis([0, 10])
colorbar

%%
figure(7)
[FX,FY] = gradient(Z);
FR = sqrt(FX.^2 + FY.^2);
h = polar([-pi pi], [0 1]);
delete(h);
hold on;
plotSmithChart;
[C,h] = contourf(X,Y,FR,30,'LineWidth',4);
hold on;
quiver(X,Y, FX, FY,10)

plot(X(1,:), Y(find(FR(:,1) == min(FR(:,1)),1),1),'LineWidth',4)

