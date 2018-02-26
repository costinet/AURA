close all;


reinit = 1;
reRunSweep = 0;
plotMeasEff = 0;
directFitDPT = 1;
plotAll = 1;
fitLossScale = 0;
fitParamVary = 0;
fitDPTcurves = 0;
plotPDT_extract = 1;


load '../Remote Operation/0-3A/testdata.mat';
MeasData1 = sortrows(LogData,[6 9]);

load '../Remote Operation/3-10A/testdata.mat';
MeasData2 = sortrows(LogData,[6 9]);

load '../Remote Operation/10-15A/testdata.mat';
MeasData3 = sortrows(LogData,[6 9]);

load '../Remote Operation/3-5A/testdata.mat';
MeasData4 = sortrows(LogData,[6 9]);

MeasData = [MeasData1; MeasData2; MeasData3; MeasData4];

Vin = MeasData(:,6);
Iin = MeasData(:,7);
Vout = MeasData(:,8);
Iout = MeasData(:,9);

baddata = Vin>30;

Vin = Vin(baddata);
Iin = Iin(baddata);
Vout = Vout(baddata);
Iout = Iout(baddata);

Pin = Vin.*Iin;
Pout = Vout.*Iout;
eff = Pout./Pin;
Ploss = Pin - Pout;

Vrange = 34:2:56;
Irange = .1:.1:14.6;

[V,I] = meshgrid(34:.5:56, .1:.02:14.6);
F = scatteredInterpolant(Vin, Iout, eff);
E = F(V,I);

F2 = scatteredInterpolant(Vin, Iout, Ploss, 'natural');
PL = F2(V,I);

% plot3(Vin, Iout, eff);
% zlim([.8 .96]);

if(plotMeasEff)
    figure(1)
    surf(V,I,E*100,'EdgeColor','none');
    hold on;   
    contour3(V,I,E*100,[.85 .9 .91 .92 .93]*100, '-k', 'LineWidth', 3);
    hold on;
    plot3(Vin, Iout, eff*100, '.k');
    zlim([.8 .96]*100);
    xlim([34.1 56]);
    ylim([.25 15]);
    xlabel('Input Voltage [V]');
    ylabel('Input Current [A]');
    zlabel('Efficiency [%]');
    caxis([.8 .96]*100);
    hold off;
end

if(reinit || reRunSweep)
    EPC9118_Digitized;
    DLFS2100_Digitized;
    
    
    %% Model
    V = 5;

    L = 1.2e-6;
    Rl = 1e-3;
    Rl_400k = .077;
    Rl_800k = .03;
    Rl_1p2M = .038;

    fs = 400e3;
    Ts = 1/fs;
    dt = 40e-9;

    

    ron_HS = 5.6e-3;%*(.45/12*dsPLoss+1); 
    % ron_HS = 7e-3*(.7/12*dsPLoss+1); 
    ron_LS = 1.8e-3;%*(.5/12*dsPLoss+1); 
    % ron_LS = 2.5e-3*(.7/12*dsPLoss+1);

    % coss_LS = 1500e-12;
    % coss_HS = 0.8e-9;

    qg_LS = 15e-9;
    qg_HS = 8.3e-9;
    Vdr = 5.1;

    rgp_HS = 2.5+.95;
    rgn_HS = 1.5+1.5;

    % % rgp_LS = 2.4+.3;
    % % rgp_LS = 1.1+.3;
    % % 
    % % Vsd_HS = 2;
    % % Rsd_HS = 11e-3;

    %Driver: 25 ns turn on with 3.3nF load
    %Driver: 13/16 ns turn-off with 3.3nF load

    Vf = 0.5;
    Rd = .0714;
    % Cd = 30e-12;

%     Pcore = 0.342;

    %% Waveform reconstruction
    % x = [Vp Il Vo]
    C = 100e-6;
    Rshunt = 100e3;

    A1 = [-1/ron_HS -1 0; 1 0 -1; 0 1 0];
    A2 = [0 -1 0; 1 0 -1; 0 1 -1/Rshunt];
    A3 = [-1/ron_LS -1 0; 1 0 -1; 0 1 0];
    A4 = A2;

    B1 = [1/ron_HS 0; 0 0; 0 -1];
    B2 = [0 0; 0 0; 0 -1];
    B3 = [0 0; 0 0; 0 -1];
    B4 = B2;

    As = cat(3, A1, A2, A3, A4);
    Bs = cat(3,B1, B2, B3, B4);
    
    VfIdT = [DLFS2100_VFIF_T25 25*ones(length(DLFS2100_VFIF_T25),1)];
    VfIdT = [VfIdT; [DLFS2100_VFIF_T75 75*ones(length(DLFS2100_VFIF_T75),1)]];
    VfIdT = [VfIdT; [DLFS2100_VFIF_T125 125*ones(length(DLFS2100_VFIF_T125),1)]];
    
    CrVr = DLFS2100_VC;
    CrVr(:,2) = CrVr(:,2)*1e-12;
    DLFS_Rd = 0;
    DLFS_Vbr = 80;
    
    Rja_HsM = 20;
    Rja_LsM = 54;
    Rja_LsD = 100;
    
    LS_Diode = Diode(DLFS_Vbr, CrVr, VfIdT, DLFS_Rd, Rja_LsD);   
    HS_MOS = MOSFET(ron_HS, qg_HS, Rja_HsM, EPC2021CossVds, EPC2021IsdVsd_T25, EPC2021RonT, 80);
    LS_MOS = MOSFET(ron_LS, qg_LS, Rja_LsM, EPC2001CossVds, 0, EPC2001RonT, 80);
    Buck = BuckConverter(V, C, L, Rl, Rl_400k, Rl_800k, fs, Ts, dt, Vdr, HS_MOS, LS_MOS, As, Bs, LS_Diode);

%     prewarp = [10 .002 .001 1000 1000 1000000 1000 1000 .01 1000 1000 1000 1000 1000 1000];
    prewarp = ones(1,15);
    if(~exist('Eswparams'))
        Eswparams = ones(1, 15);
        prewarp(4:8) = prewarp(4:8)*.00001;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(reRunSweep)
        h = waitbar(0/length(Vin),['Measuring datapoint ' num2str(1) ' of ' num2str(length(Vin)) '. Vin = ' num2str(Vin(1)) ', Iout = ' num2str(Iout(1))]);

        nmeas = length(Vin);
        
        Ploss_model = zeros(1, nmeas);
        Ion_s = zeros(1, nmeas);
        Ioff_s = zeros(1, nmeas);
        Pcond = zeros(1, nmeas);
        Pconst = zeros(1, nmeas);
        Poss = zeros(1, nmeas);
        Esw = zeros(1, nmeas);


        for i = 1:nmeas
                if(Iout(i) > .01 && Vin(i) > 30)

%                     [ Pcond(i), Pg(i), Pbd(i), Poss(i), Pq(i), Pov(i), Pboot(i), Pcore(i), Voavg(i), Temps(:,i)] = EPC9118WhiteBox( Iout(i), Vin(i), Buck);
                        [ Ploss_model(i), Ion_s(i), Ioff_s(i), Pcond(i), Pconst(i), Poss(i), Esw(i)] = ...
                          EPC9118WhiteBox_DPT( Eswparams, [Vin(i) Iout(i)], Buck, prewarp);
                end
            h = waitbar(i/length(Vin),h,['Measuring datapoint ' num2str(i) ' of ' num2str(length(Vin)) '. Vin = ' num2str(Vin(i),3) ', Iout = ' num2str(Iout(i),3)]);
        end
        close(h);
    end
end


options = optimoptions(@lsqnonlin,'TolFun',1e-16, 'DiffMinChange', 1e-9, 'Display', 'iter');

if(directFitDPT)
    DPTLoss = Ploss'-Poss-Pcond - Pconst;
    ind = find(Iout>0);
    ys_fit = [Ioff_s(ind); Ion_s(ind); Vin(ind)'];
    DPTLoss = DPTLoss(ind);
      
    for i = 1:5
        x0 = Eswparams;
%         prewarp = min(1./Eswparams.*prewarp/mean(prewarp),1e-6);
        [x2, resnorm] = lsqnonlin(@(X)fitDPT(Buck, ys_fit, X, prewarp)-DPTLoss, x0, x0*0, 10*ones(size(x0)), options);
        Eswparams = x2;
    end 
    
    [cs, ks] = getCKs(Buck, Eswparams, prewarp);
    [HSEon, HSEoff] = solveEsw(Buck, Ion_s, Ioff_s, Ion_s, cs, ks, Vin');
    Esw = HSEon + HSEoff;
    
    Ploss_model =  Pcond + Pconst + Poss + Esw*Buck.fs*1e-6;
end

if(plotAll)

    if(size(Ploss_model, 2) > 1)
        Ploss_model = Ploss_model';
    end

    [V,I] = meshgrid(34:.5:56, .1:.02:14.6);

    F2_model = scatteredInterpolant(Vin, Iout, Ploss_model, 'natural');
    PL_model = F2_model(V,I);

    figure(2)
%     subplot(3,1,1)
    %% Modeled vs Analytical Losses
%     surf(V,I,PL,'EdgeColor','none');
    hold on;   
%     contour3(V,I,PL,.25:.25:8, '-k', 'LineWidth', 3);
    hold on;
     surf(V,I,PL_model,'EdgeColor','none');
     alpha(.5);
    plot3(Vin, Iout, Ploss, '.k');
   
%     plot3(Vin, Iout, Ploss_model, '*r');
    zlim([0 8.5]);
    xlim([34.1 56]);
    ylim([.25 15]);
    xlabel('Input Voltage [V]');
    ylabel('Input Current [A]');
    zlabel('Power Loss [W]');
%     legend('Measured Power Loss', 'Constant Loss Contours', 'Measured Points', 'Model Predictions', 'Location', 'Best')
    legend('Measured Power Loss', 'Measured Points', 'Model Predictions', 'Location', 'Best')
    view(-45,45)
    hold off;
    box on;
    grid on;
    camlight headlight

%     subplot(3,1,2)
%     %% Modeled Losses
%     surf(V,I,PL_model,'EdgeColor','none');
%     hold on;
% %     contour3(V,I,PL_model,.25:.25:8, '-k', 'LineWidth', 3);
%     zlim([0 8])
%     title('Analytical Power Loss')
%     xlabel('Input Voltage [V]');
%     ylabel('Input Current [A]');
%     zlabel('Power Loss [W]');
% %     camlight headlight
figure(3)
%     subplot(3,1,3)
    %% Model Error 
    resid = (PL - PL_model);%./PL*100;
    resid_filt = filter2(fspecial('gaussian', [10 10], 20), resid);
    surf(V,I,resid_filt,'EdgeColor','none');
    hold on;
    % surf(V,I,resid_filt-resid,'EdgeColor','none');
    % hold on;
    % contour3(V,I,resid,-10:5:30, '-k', 'LineWidth', 3);

%     set(gca, 'CameraPosition', [5 15/2 min(min(resid))])
    view(90,0)
    camlight headlight


    figure(5)
%     subplot(2,1,1)
    %% All Loss Surface Comparison
    V = V(1:10:size(V,1),:);
    I = I(1:10:size(I,1),:);
    Poss_interp = scatteredInterpolant(Vin, Iout, Poss', 'natural');
    Poss_model = Poss_interp(V,I);
    Psw_interp = scatteredInterpolant(Vin, Iout, Esw'*Buck.fs*1e-6, 'natural');
    Psw_model = Psw_interp(V,I);
    Pcond_interp = scatteredInterpolant(Vin, Iout, Pcond', 'natural');
    Pcond_model = Pcond_interp(V,I);

%     subplot(3,1,1)
    surf(V,I,Poss_model,'EdgeColor','none');
    box on; grid on;
    hold on;
    zlim([0 5])
    xlim([30 60])
    ylim([0 15])
    caxis([0 5])
%     subplot(3,1,2)
    surf(V,I,Psw_model,'EdgeColor','none');
%     axis off
%     subplot(3,1,3)
    zlim([0 5])
    xlim([30 60])
    ylim([0 15])
    caxis([0 5])
    surf(V,I,Pcond_model,'EdgeColor','none');
%     axis off
    zlim([0 5])
    xlim([30 60])
    ylim([0 15])
    caxis([0 5])
    alpha(0.5)
    
    Ion_interp = scatteredInterpolant(Vin, Iout, Ion_s', 'natural');
    Ion_model = Ion_interp(V,I);
    Ioff_interp = scatteredInterpolant(Vin, Iout, Ioff_s', 'natural');
    Ioff_model = Ioff_interp(V,I);
    
    
    
%     
%     subplot(2,1,2)
%     %% DPT Loss Surface Comparison
%     surf(V,I,PL-Poss_model-Pcond_model, 'EdgeColor','none');
%     hold on;
%      surf(V,I,Psw_model,'EdgeColor','none');
%      
%      legend('Measured Data - Poss - Pcond', 'Derived DPT Losses');
%      
%     
%     


    %% 36 and 56 V curves
    figure(6)
    for i = 1:2
        subplot(2,1,i)
        if(i==1)
            indexes = find(Vin > 55.5 & Vin <56.5); % 56
        else
            indexes = find(Vin > 35.5 & Vin <36.5); % 36
        end
        plot(Iout(indexes), Ploss(indexes), 'o')
        hold on;
        Ploss_corrected = Ploss_model(indexes)';
        plot(Iout(indexes), Ploss_corrected, 'o')

        plot(Iout(indexes), Ion_s(indexes) > 0, ':')
        plot(Iout(indexes), Esw(indexes)*Buck.fs*10^-6, '-')

        % ylim([.85 .95])
        legend('meas', 'calc')
    end
    
    %% Comparison to datasheet values
    figure(11)
    I = imread('../EPC9118-noaxes.png');
    % hi = imagesc(I);
    xrng = [0 25];
    yrng = [80 96];
    image(xrng,yrng,flipud(I));
    set(gca,'YDir','normal');
    % plot(eta_36V(:,1),eta_36V(:,2),'*', 'linewidth',5);
    % plot(eta_48V(:,1),eta_48V(:,2),'o', 'linewidth',5)
    % plot(eta_56V(:,1),eta_56V(:,2),'d', 'linewidth',5)
    xlim([0 25]);
    xlabel('I_{out} [A]');
    ylabel('Efficiency [%]');
    hold on;
    for i = 1:3
        if(i==1)
            indexes = find(Vin > 55.5 & Vin <56.5); % 56
        elseif(i==2)
            indexes = find(Vin > 47.5 & Vin <48.5); % 48
        else
            indexes = find(Vin > 35.5 & Vin <36.5); % 36
        end
        plot(Iout(indexes), 100*5*Iout(indexes)./(5*Iout(indexes)+Ploss(indexes)), 'o')
    end
%     ylim([80 98])
   
    

end


%% lsqnonlin varying DPT curves
if(fitDPTcurves)
    x0 = Eswparams;
    VI = [Vin Iout];
    VIShort = VI(1:5:length(VI),:);
    PlossShort = Ploss(1:5:length(Ploss));

    [x, resnorm] = lsqnonlin(@(X)EPC9118WhiteBox_DPT(X, VIShort, Buck, prewarp)-PlossShort, x0, [x0(1:end-1)*0 -2], 3*ones(size(x0)), options);
end

%% Extracted DPT Data
if(plotPDT_extract)
    figure(10)
    
    [cs, ks] = getCKs(Buck, Eswparams, prewarp);
    
    Vdcs = 38:6:60;
    Ios = -5:.1:20;
    
    [Vtest, Itest] = meshgrid(Vdcs, Ios);
    Vtest = Vtest(:);
    Itest = Itest(:);

    [HSEon, HSEoff] = solveEsw(obj, Itest, Itest, Itest, cs, ks, Vtest);
    plot(Itest, HSEon, 'LineWidth',2);
    hold on;
    plot(Itest, HSEoff, 'LineWidth',2);
    xlabel('Current [A]');
    ylabel('Switching Energy [{\mu}J]')
    legend('E_{on}', 'E_{off}')
        box on;
    grid on;
end


    

    

