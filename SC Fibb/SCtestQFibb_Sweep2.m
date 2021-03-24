Vg = 12;
j=1;
Vgrange = 5:1:24;
Vb1 = 4;
Vb2 = 4;

Rb = 5e-3;
RL = 5e-3;

ron = 5e-3;%8.5e-3;
% Coss = 2e-9;
Coss = 0;%.9e-9;

Co = 2e-6; ESRo = 2e-3;
Cfly1 = 5e-6; ESR1 = 3e-3;  % 5V Cap
Cfly2 = 2.5e-6; ESR2 = 3e-3; % 10V Cap
Lc = 100e-9;

fs = 2e6;
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

circuitPath = 'SCreg2SCharger/Q-FibonacciESRCoss';

sim = SMPSim();
conv = sim.converter;
top = sim.topology;

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
drange = linspace(0, .99, 10);
h = waitbar(0, 'Simulating');
for i = 1:size(modSchemes,3)
    swvec = modSchemes(:,:,i);
    waitbar((i-1)/size(modSchemes,3), h);
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
        MaxVin = OLVin+1;
        
        plecs('set', [circuitPath '/R3'], 'R', 'RL')
        top.loadPLECsModel(circuitPath,swvec);
        
        waitbar((i-1)/size(modSchemes,3) + 1/size(modSchemes,3)*(find(drange==d)-1)/length(drange) , h);

        Vinrange = OLVin:.02:MaxVin;
        for Vin = Vinrange
            u(Vgloc) = Vin;
            sim.u = u';

            Xss = sim.SS_Soln();
            [ avgXs, avgYs ] = sim.ssAvgs(Xss);
%             avgYs(1:2)
    %     sim.plotAllStates(1)
    %     sim.plotAllOutputs(2)
    
%             swActions = sum(abs(diff([swvec; swvec(1,:)],1)))/2;
%             Vblock = [Vb1 Vb1 Vb1 Vb2 Vb2 Vb2 Vb1 Vb1 Vb2 Vb2 2*Vb1 2*Vb1 2*Vb1 2*Vb2 2*Vb2 2*Vb2];
%             Poss = sum(1/2*Coss*Vblock.^2.*swActions*fs);
            Poss = 0;
            

            eta = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2 - Poss)/(Vin*avgXs(Illoc) - avgYs(Ib1loc)^2*Rb - avgYs(Ib2loc)^2*Rb);
            Ploss = -(avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) + (Vin*avgXs(Illoc)) + Poss;
            Pout = (avgYs(Ib1loc)*Vb1 + avgYs(Ib2loc)*Vb2) - Poss;
            
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

close(h)

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

mineff = .5;
etaSim(etaSim<mineff) = mineff;
etaSim(etaSim>1) = mineff;

locs = etaSim > mineff;

VgSim = conditions(locs,end);
F = scatteredInterpolant(VgSim, PoutSim(locs), etaSim(locs)*100, 'linear','none');

figure(1)

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

figure(2)


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

[c] = scatter(conditions(:,end),PoutSim,[],etaSim*100);
xlabel('Vg');
ylabel('P_{out}');
ylim([0 100])
colorbar




