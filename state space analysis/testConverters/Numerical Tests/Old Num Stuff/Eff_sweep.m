function [stick] = Eff_sweep(X)

% Inside GA Loop
% This is going to make a eff plot for the converter Init

if length(X) == 1
    X_adjust_finmoncon = [0.1859    0.1548    0.3997    0.4000    0.0196    0.1075    0.0794  0.3187].*[20 20 20 20 20 20 20 20];
    X = [X_adjust_finmoncon X];
end

try
    penalty  = [];
    adjust = [ones(1,length(X)-1)*20, 1e-6] ;
    X = X./adjust;
    
    W_selection = X(1:end-1);
    fs = X(end);
    
    ron = [];
    Coss = [];
    tot_area = 0;
    
    %% Determining Coss and ron based on w
    if length(W_selection)==8
        FET_selection = [2 2 2 3 3 3 3 3];
    elseif length(W_selection)==9
        FET_selection = [2 2 2 3 3 3 3 3 3];
    else
        fprintf('Something went wrong /n')
    end
    
    for select_FET = 1:length(W_selection)
        
        switch FET_selection(select_FET)
            
            %% 8HVnLDMOS nbl
            case 1
                a =  0.001378;
                b = -1;
                ron(select_FET) = a*W_selection(select_FET)^b;
                
                p1 = 2.925e-9;
                p2 = -2.787e-12;
                Coss(select_FET) = p1*W_selection(select_FET)+p2;
                
                L1=800*10^-9;
                
                tot_area = tot_area+W_selection(select_FET)*L1;
                
                
                %% 8HVnLDMOS iso
            case 2
                a =  0.001453;
                b = -1;
                ron(select_FET) = a*W_selection(select_FET)^b;
                
                p1 = 2.093e-9;
                p2 = -1.058e-12;
                Coss(select_FET) = p1*W_selection(select_FET)+p2;
                
                L2=900*10^-9;
                
                tot_area = tot_area+W_selection(select_FET)*L2;
                
                
                %% 12HVnLDMOS iso hp mac
            case 3
                a =  0.001447;
                b = -1;
                ron(select_FET) = a*W_selection(select_FET)^b;
                
                p1 = 2.487e-9;
                p2 = -6.393e-13;
                Coss(select_FET) = p1*W_selection(select_FET)+p2;
                
                L3=900*10^-9;
                
                tot_area = tot_area+W_selection(select_FET)*L3;
                
                
                %% 12HVnLDMOS iso mac
            case 4
                a =  0.001656;
                b = -1;
                ron(select_FET) = a*W_selection(select_FET)^b;
                
                p1 = 2.281e-9;
                p2 = -3.1e-12;
                Coss(select_FET) = p1*W_selection(select_FET)+p2;
                
                L4=1*10^-6;
                
                
                tot_area = tot_area+W_selection(select_FET)*L4;
                
                %% 12HVnLDMOS nbl hp mac
            case 5
                a =  0.001447;
                b = -1;
                ron(select_FET) = a*W_selection(select_FET)^b;
                
                p1 = 2.487e-9;
                p2 = -6.393e-13;
                Coss(select_FET) = p1*W_selection(select_FET)+p2;
                
                L5=0.9*10^-6;
                
                tot_area = tot_area+W_selection(select_FET)*L5;
                
                
                %% 12HVnLDMOS nbl mac
            case 6
                a =  0.001585;
                b = -1;
                ron(select_FET) = a*W_selection(select_FET)^b;
                
                p1 = 3.176e-9;
                p2 = -7.226e-12;
                Coss(select_FET) = p1*W_selection(select_FET)+p2;
                
                L6=1*10^-6;
                
                tot_area = tot_area+W_selection(select_FET)*L6;
                
                %% 12HVnLDMOS nbl mr mac
            case 7
                a =  0.002556;
                b = -1;
                ron(select_FET) = a*W_selection(select_FET)^b;
                
                p1 = 2.828e-9;
                p2 = -1.703e-12;
                Coss(select_FET) = p1*W_selection(select_FET)+p2;
                
                L7=.9*10^-6;
                
                tot_area = tot_area+W_selection(select_FET)*L7;
                
            otherwise
                
                stick = 100;
                return
                
        end
        
    end
    
    if 2*tot_area > (1.00001e-6)*100
        
        penalty = 2*tot_area*1e6;
        
    else
        penalty = 0;
        
    end
    
    %% Setting up Vales
    
    u = [];
    Vg = 12;
    j=1;
    Vgrange = 5:1:24;
    Vb1 = 4;
    Vb2 = 4;
    
    Rb = 5e-3;
    RL = 5e-3;
    
    %ron = 5e-3;%8.5e-3;
    % Coss = 2e-9;
    %Coss = 0;%.9e-9;
    
    Co = 2e-6; ESRo = 2e-3;
    Cfly1 = 5e-6; ESR1 = 3e-3;  % 5V Cap
    Cfly2 = 2.5e-6; ESR2 = 3e-3; % 10V Cap
    Lc = 100e-9;
    
    % fs = 2e6;
    Ts = 1/fs;
    FETs = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    h3Lb_1x = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    h3Lb_2x = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    h3Lb_bp = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    l3Lb_1x = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
    l3Lb_2x = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
    l3Lb_bp = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
    
    h3Lu_bp = [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
    h3Lu_qs = [0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0];
    % h3Lu_m4 = [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
    h3Lu_0 =  [0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0];
    
    l3Lu_bp = [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
    l3Lu_qs = [0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0];
    % l3Lu_m4 = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
    l3Lu_0 =  [0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0];
    
    h3Lm_1x = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0];
    h3Lm_2x = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0];
    h3Lm_bp = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
    
    l3Lm_1x = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];
    l3Lm_2x = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
    l3Lm_bp = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    
    u0b = h3Lu_0 + h3Lb_1x;
    u0 = h3Lu_0 + h3Lm_1x + h3Lb_1x;
    u4 = h3Lu_bp + h3Lm_1x + h3Lb_1x;
    u8 = h3Lu_bp + h3Lm_bp + h3Lb_2x;
    u12 = h3Lu_bp + h3Lm_2x + h3Lb_1x;
    u16 = h3Lu_bp + h3Lm_2x + h3Lb_2x;
    
    l0b = l3Lu_0 + l3Lb_1x;
    l0 = l3Lu_0 + l3Lb_1x + l3Lm_1x;
    l4 = l3Lu_bp + l3Lm_1x + l3Lb_1x;
    l8 = l3Lu_bp + l3Lm_bp + l3Lb_2x;
    l12 = l3Lu_bp + l3Lm_2x + l3Lb_1x;
    l16 = l3Lu_bp + l3Lm_2x + l3Lb_2x;
    
    
    % swvec = [u12 + l12; u8b + l12; u4 + l12; u12 + l12; u12 + l8b; u12 + l4];
    % ds = [.5, 1, 1, .5, 1, 1];
    
    Dx =.4;
    
    levels = [0,0; 0,4; 4,4; 4,8; 4,12; 4,16; 12,12];
    
    modSchemes = zeros(4,length(u0), size(levels,1));
    
    for ii = 1:size(levels,1)
        if levels(ii,1) == levels(ii,2)
            modSchemes(1,:,ii) = eval(['u',num2str(levels(ii,1))]) + eval(['l',num2str(levels(ii,2))]);
            modSchemes(3,:,ii) = eval(['u',num2str(levels(ii,1))]) + eval(['l',num2str(levels(ii,2))]);
        else
            modSchemes(1,:,ii) = eval(['u',num2str(levels(ii,1))]) + eval(['l',num2str(levels(ii,2))]);
            modSchemes(3,:,ii) = eval(['l',num2str(levels(ii,1))]) + eval(['u',num2str(levels(ii,2))]);
        end
        if ii >1
            modSchemes(2,:,ii-1) = modSchemes(1,:,ii);
            modSchemes(4,:,ii-1) = modSchemes(3,:,ii);
        end
    end
    
      modSchemes(:,:,end) = [];
    
    modSchemes(:,[7:8],:) = [];
    
    swvec = modSchemes(:,:,1);
    % % if Vg == 21
    % %     swvec = [u12 + l12; u12 + l4; l12 + u12; u4 + l12];
    % % elseif Vg <= 24 && Vg> 20
    % %     swvec = [u12 + l12; u16 + l4; l12 + u12; u4 + l16];
    % % elseif Vg == 17 || Vg == 18
    % %     swvec = [u16 + l4; u8 + l4; l16 + u4; u4 + l8];
    % % elseif Vg <= 20 && Vg> 16
    % %     swvec = [u16 + l4; u8 + l8; l16 + u4; u8 + l8];
    % % elseif Vg <= 16 && Vg > 12
    % %     swvec = [u12 + l4; u8 + l4; l12 + u4; u4 + l8];
    % % elseif Vg <= 12 && Vg > 8
    % %     swvec = [u8 + l4; u4 + l4; l8 + u4; u4 + l4];
    % % elseif Vg <= 8 && Vg > 4
    % %     swvec = [u4 + l4; u4 + l0; l4 + u4; u0 + l4];
    % % end
    % %
    % % llim = floor((Vg-.001)/4)*4;
    % % ds = [ (Vg-llim)/4, 1-(Vg-llim)/4+Dx, (Vg-llim)/4, 1-(Vg-llim)/4+Dx];
    % %
    % % swvec(:,[7:8]) = [];
    
    
    % ts = ds/sum(ds)*Ts;
    
    circuitPath = 'Baxter_SCreg2SCharger/QFibonacciESRCoss';
    %circuitPath = 'SCreg2SCharger/2StageESRCoss';
    %circuitPath = 'SCreg2SCharger/3StageESRCoss';
    
    sim = SMPSim();
    conv = sim.converter;
    top = sim.topology;
    
    
    %% Assigning value to the base workspace
    
    % Move DAB parameters to base workspace
    % PLECS will pull from the base workspace, so we need to redefine the
    % variables in that workspace (not this function's) before parsing the
    % state space matrices with PLECS.  An alternative is to use
    % plece('set',...) commands
    
    
    %% FET Pairs
    %{
M1 and M5    (1)
M2 and M4    (2)
M3 and M6    (3)
M9 and M12   (4)
M10 and M11  (5)
M13 and M17  (6)
M14 and M16  (7)
M15 and M18  (8)
    %}
    
    
    
    vars = {'Vg', 'Vb1', 'Vb2', 'Lc', 'Rb', 'RL', 'ESR1','ESR2', 'Cfly1', 'Cfly2','Co', 'ESRo', 'Coss', 'ron'};
    for i = 1:length(vars)
        assignin('base', vars{i},eval(vars{i}));
    end
    
    % This does not work here. I believe the value need to be set in the
    % base workspace
    %{
plecs('set', [circuitPath '/Vg'], 'V', 'Vg')
plecs('set', [circuitPath '/Vb1'], 'V', 'Vb1')
plecs('set', [circuitPath '/Vb2'], 'V', 'Vb2')

plecs('set', [circuitPath '/L1'], 'L', 'Lc')
plecs('set', [circuitPath '/L2'], 'L', 'Lc')

plecs('set', [circuitPath '/R1'], 'R', 'Rb')
plecs('set', [circuitPath '/R2'], 'R', 'Rb')
plecs('set', [circuitPath '/R3'], 'R', 'RL')
plecs('set', [circuitPath '/R4'], 'R', 'RL')
plecs('set', [circuitPath '/R5'], 'R', 'ESR1')
plecs('set', [circuitPath '/R6'], 'R', 'ESR2')
plecs('set', [circuitPath '/R7'], 'R', 'ESR1')
plecs('set', [circuitPath '/R8'], 'R', 'ESR2')
plecs('set', [circuitPath '/R9'], 'R', 'ESRo')
plecs('set', [circuitPath '/R10'], 'R', 'ESRo')

plecs('set', [circuitPath '/C1'], 'C', 'Cfly1')
plecs('set', [circuitPath '/C2'], 'C', 'Cfly1')


plecs('set', [circuitPath '/C5'], 'C', 'Cfly2')
plecs('set', [circuitPath '/C6'], 'C', 'Cfly2')
plecs('set', [circuitPath '/C6'], 'C', 'Cfly2')
plecs('set', [circuitPath '/C6'], 'C', 'Cfly2')
plecs('set', [circuitPath '/C6'], 'C', 'Cfly2')
    %}
    
    top.loadPLECsModel(circuitPath,swvec);
    
    Ib1loc = find(strcmp(top.outputLabels, 'Ib1'));
    Ib2loc = find(strcmp(top.outputLabels, 'Ib2'));
    Vscloc = find(strcmp(top.outputLabels, 'Vsc'));
    
    
    Vgloc = find(strcmp(top.inputLabels, 'Vg'));
    Vb1loc = find(strcmp(top.inputLabels, 'Vb1'));
    Vb2loc = find(strcmp(top.inputLabels, 'Vb2'));
    
    Illoc = find(strcmp(top.stateLabels, 'L1'));
    Vc1loc = find(strcmp(top.stateLabels, 'C1'));
    Vc2loc = find(strcmp(top.stateLabels, 'C2'));
    Vc3loc = find(strcmp(top.stateLabels, 'C3'));
    Vc4loc = find(strcmp(top.stateLabels, 'C4'));
    
    u(Vgloc) = Vg;
    u(Vb1loc) = Vb1;
    u(Vb2loc) = Vb2;
    sim.u = u';
    % conv.ts = ts;
    
    
    
    etaSim = [];
    PlossSim = [];
    PoutSim = [];
    conditions = [];
    PossSim = [];
    drange = linspace(0.05, .95, 8);
    %h = waitbar(0, 'Simulating');
    
   for i = 2:size(modSchemes,3)
  %  for i = 3:6
        swvec = modSchemes(:,:,i);
        %  waitbar((i-1)/size(modSchemes,3), h);
        for d = drange
            ds = [d, 1-d, d, 1-d];
            ts = ds/sum(ds)*Ts;
            conv.ts = ts;
            
            %% Find zero-power Vin
            plecs('set', [circuitPath '/R3'], 'R', '1e6')
            top.loadPLECsModel(circuitPath,swvec);
            Xss = sim.SS_Soln();
            [ avgXs, avgYs ] = sim.ssAvgs(Xss);
            
            OLVin = avgYs(Vscloc);
           % MaxVin = OLVin+1;
            MaxVin = OLVin+2;
            if OLVin<5
                continue
            end
            %{
        if i==2 && MaxVin>10
            continue
        end
        
        if i==3 && MaxVin>15
            continue
        end
        
        if i==4 && MaxVin>16 && OLVin<13
            continue
        end
        
        if i==5 && MaxVin>23 && OLVin<18
            continue
        end
        
        if i==6 && MaxVin>27 && OLVin<19
            continue
        end
            %}
            
            plecs('set', [circuitPath '/R3'], 'R', 'RL')
            top.loadPLECsModel(circuitPath,swvec);
            
            % waitbar((i-1)/size(modSchemes,3) + 1/size(modSchemes,3)*(find(drange==d)-1)/length(drange) , h);
            
            Vinrange = OLVin:.1:MaxVin;
            for Vin = Vinrange
                u(Vgloc) = Vin;
                sim.u = u';
                
                Xss = sim.SS_Soln();
                [ avgXs, avgYs ] = sim.ssAvgs(Xss);
                
                [ xs, t, ys ] = sim.SS_WF_Reconstruct;
                
                RMS_Iin = rms(xs(3,:));
                RMS_Ib1 = rms(ys(1,:));
                RMS_Ib2 = rms(ys(2,:));
                

                
                %             avgYs(1:2)
                %     sim.plotAllStates(1)
                %     sim.plotAllOutputs(2)
                
                %             swActions = sum(abs(diff([swvec; swvec(1,:)],1)))/2;
                %             Vblock = [Vb1 Vb1 Vb1 Vb2 Vb2 Vb2 Vb1 Vb1 Vb2 Vb2 2*Vb1 2*Vb1 2*Vb1 2*Vb2 2*Vb2 2*Vb2];
                %             Poss = sum(1/2*Coss*Vblock.^2.*swActions*fs);
                Poss = 0;
                
                % for when there is no rms calculation
                %{
                eta = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2 - Poss)/(Vin*avgXs(Illoc) - avgYs(Ib1loc)^2*Rb - avgYs(Ib2loc)^2*Rb);
                Ploss = -(avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) + (Vin*avgXs(Illoc)) + Poss;
                Pout = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) - Poss;
                %}
                
                % True input to output eff output of converter
                
                P_calc_Vin = Vin-(RMS_Iin*(RL*2));
                P_calc_Vout1 = Vb1+RMS_Ib1*Rb;
                P_calc_Vout2 = Vb2+RMS_Ib2*Rb;
                
                eta = (P_calc_Vout1*avgYs(Ib1loc) + P_calc_Vout2*avgYs(Ib2loc)) / (P_calc_Vin*avgXs(Illoc));
                Ploss = -(P_calc_Vout1*avgYs(Ib1loc)+P_calc_Vout2*avgYs(Ib2loc))+(P_calc_Vin*avgXs(Illoc));
                Pout = (P_calc_Vout1*avgYs(Ib1loc)+P_calc_Vout2*avgYs(Ib2loc));
                
                if eta>1||eta<.3
                    continue
                end
                
                if Pout<0 || Ploss<0
                    continue
                end
                
                % for specific power
                %{
                if Pout<40 || Pout>60
                    continue
                end
                %}
                
                etaSim = [etaSim; eta];
                PlossSim = [PlossSim; Ploss];
                PoutSim = [PoutSim; Pout];
                %             PossSim = [PossSim; Poss];
                conditions = [conditions; d, i, Vin];
                
                if Pout > 100
                    break;
                end
                
                %             effSim(i,j) = eta;
                %             PlossSim(i,j) = Ploss;
                %             PoutSim(i,j) = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2);
            end
            
            
            % %         %% Modeling
            % %         Ig = avgXs(Illoc);
            % %
            % %         dVc5 = Ig*sum(swvec(:,13-2).*ts')/Cfly2;
            % %         dVc1 = Ig*sum(swvec(:,1).*ts')/Cfly1;
            % %
            % %         Pdir = Vb1*Ig*sum(swvec(:,9-2).*ts')/Ts + Vb2*Ig*sum(swvec(:,12-2).*ts')/Ts;
            % %         Plossdir = 2*Ig^2*( ron*sum(swvec(:,9-2).*ts') + ...
            % %             ron*sum(swvec(:,13-2).*ts') + ...
            % %             2*ron*sum((1-swvec(:,13-2)).*ts') + ...
            % %             ron*sum(swvec(:,1).*ts') + ...
            % %             2*ron*sum((1-swvec(:,1)).*ts') + ...
            % %             Rb*Ts ...
            % %             )/Ts;
            % %
            % %
            % %         Psc = 2*1/2*Cfly2*((8+dVc5)^2 - 8^2)*fs + 2*1/2*Cfly1*((4+dVc1)^2 - 4^2)*fs;
            % %         PlossSc = 2*1/2*Cfly2*(dVc5^2)*fs + 2*1/2*Cfly1*(dVc1^2)*fs;
            % %
            % %         Poutsim = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2);
            % %         Poutmodel = Pdir + Psc;
            % %
            % %     %     Plosssim = Ploss;
            % %     %     Plossmodel = Plossdir + PlossSc;
            % %
            % %         PoutModel(i,j) = Pdir + Psc;
            % %         PlossModel(i,j) =  Plossdir + PlossSc;
            % %         effModel(i,j) = Poutmodel/(Poutmodel + Plossdir + PlossSc);
            
            
        end
    end
    
    %close(h)
    
    % figure(4)
    % subplot(2,1,1)
    % plot(PoutSim, PlossSim, 'ok');
    % hold on;
    % plot(PoutModel, PlossModel)
    %
    % subplot(2,1,2)
    % plot(PoutSim, effSim, 'ok');
    % hold on;
    % plot(PoutModel, effModel)
    %
    
    % % [f,c] = contourf(repmat(Vgrange,length(drange),1),PoutSim,effSim*100);
    % % xlabel('Vg');
    % % ylabel('P_{out}');
    % % ylim([0 100])
    % % colorbar
    % % c.LevelList = [.8 .9 .95 .96 .97 .98 .99]*100;
    %{
mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), etaSim(locs)*100, 'linear','none');

figure(101)

Vgrange = 5:.25:22;
PoutRange = 0:1:80;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
% c.LevelList = [.8 .9 .95 .96 .97 .98 .99]*100;

hold on;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');

figure(102)

Vgrange = 5:.25:22;
PoutRange = 0:1:80;
[VgMesh, PoutMesh] = meshgrid(Vgrange, PoutRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
c.LevelList = [0 1 2 3 4 5 6 7 8 9];

hold on;


Vq = interp2(Vgrange,PoutRange,F(VgMesh,PoutMesh,,Yq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)


%
% figure(2)
% subplot(2,1,1)
% F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');
% [f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
% levels = c.LevelList;
% colorbar
%
% subplot(2,1,2)
% F = scatteredInterpolant(VgSim, PoutSim(locs), PossSim(locs), 'linear','none');
% [f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
% colorbar
%
F2 = scatteredInterpolant(VgSim, PoutSim(locs), conditions(locs,2), 'linear','none');
contourf(Vgrange,PoutRange,F2(VgMesh,PoutMesh),'ShowText','on');
scatter(VgSim, PoutSim(locs),[],conditions(locs,2));
colorbar
[c] = scatter(conditions(:,end),PoutSim,[],etaSim*100);
xlabel('Vg');
ylabel('P_{out}');
ylim([0 100])
colorbar

    %}
    
    mineff = .3;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), PlossSim(locs), 'linear','none');
    

Vgrange = 10:.25:20;
PoutRange = 20:1:60; 
    
sigma1 = 30;
sigma2 = 750;

mu = [16 50];
Sigma = [sigma1, sqrt(sigma1*sigma2-1); sqrt(sigma1*sigma2-1), sigma2];

[X1,X2] = meshgrid(Vgrange',PoutRange');
X = [X1(:) X2(:)];

p = mvncdf(X,mu,Sigma);

Z = reshape(p,length(PoutRange),length(Vgrange));

weighted_vals=F(X1,X2).*(Z+Z([length(PoutRange):-1:1],[length(Vgrange):-1:1]));

    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    Added_ploss = [];
    stick = [];
    stick = mean(weighted_vals,'all');
    stick = stick+penalty;
    %{
    Added_ploss = (1*(median(PlossSim(conditions(:,3)>18 & PoutSim<20)))) + ...
        (4*(median(PlossSim(conditions(:,3)>18 & PoutSim>60)))) + ...
        (3*(median(PlossSim(conditions(:,3)>18 & PoutSim<60 & PoutSim>40)))) + ...
        (2*(median(PlossSim(conditions(:,3)>18 & PoutSim<40 & PoutSim>20)))) + ...
        (4*(median(PlossSim(conditions(:,3)<10 & PoutSim<20)))) + ...
        (0*(median(PlossSim(conditions(:,3)<10 & PoutSim>60)))) + ...
        (1*(median(PlossSim(conditions(:,3)<10 & PoutSim<60 & PoutSim>40)))) + ...
        (2*(median(PlossSim(conditions(:,3)<10 & PoutSim<40 & PoutSim>20)))) + ...
        (3*(median(PlossSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim<20)))) + ...
        (1*(median(PlossSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim>60)))) + ...
        (2*(median(PlossSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim<60 & PoutSim>40)))) + ...
        (4*(median(PlossSim(conditions(:,3)>10 & conditions(:,3)<14 & PoutSim<40 & PoutSim>20)))) + ...
        (2*(median(PlossSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim<20)))) + ...
        (2*(median(PlossSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim>60)))) + ...
        (4*(median(PlossSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim<60 & PoutSim>40)))) + ...
        (3*(median(PlossSim(conditions(:,3)>14 & conditions(:,3)<18 & PoutSim<40 & PoutSim>20))));
    
    stick = Added_ploss/38;
    %}
    
   % stick = median(PlossSim);
    
    
    
    if isnan(stick)||isinf(stick)||stick<0
        stick = 100;
    end
catch ME
    stick = 100;
end
end