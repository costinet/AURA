
load('EPC_GaN_Table.mat')
load('Sim_exp_3level.mat')
load('Inductor_table.mat')

EPCeGaNFETselectorguide1(1:2,:) =[];

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create scatter3
scatter3(PowerInductorS1.DCRTyp25Cm,PowerInductorS1.LNominalnH,PowerInductorS1.Heightmm,100,'MarkerEdgeColor','none',...
    'MarkerFaceColor',[1 0.509803950786591 0]);

% Create zlabel
zlabel('Height (mm)','FontName','Times New Roman');

% Create ylabel
ylabel('Inductance (nH)','FontName','Times New Roman');

% Create xlabel
xlabel('DC Resistance (m\Omega)','FontName','Times New Roman');

view(axes1,[-224.146874082743 22.8484262048946]);
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',24);



%%


% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

%scatter(FETs(:,1),FETs(:,2),100,'MarkerEdgeColor','none',...
%    'MarkerFaceColor',[1 0.509803950786591 0])
scatter(EPCeGaNFETselectorguide1.MaxRDSonm5VGS,EPCeGaNFETselectorguide1.COSSpF,100,'MarkerEdgeColor','none',...
   'MarkerFaceColor',[1 0.509803950786591 0],'DisplayName','Discrete Devices')

% Create ylabel
ylabel('C_{oss} (pF)','FontName','Times New Roman');

% Create xlabel
xlabel('R_{ds,on} (m\Omega)','FontName','Times New Roman');

set(axes1,'FontName','Times New Roman','FontSize',24);

hold on
plot([FETs(2,1) FETs(2,1)],[FETs(2,2) FETs(2,2) ],'LineWidth',5,'DisplayName','M1')
plot([FETs(3,1) FETs(3,1)],[FETs(3,2) FETs(3,2) ],'LineWidth',5,'DisplayName','M2')
plot([FETs(4,1) FETs(9,1) ],[ FETs(4,2) FETs(9,2)],'LineWidth',5,'DisplayName','M3')
plot([FETs(2,1) FETs(13,1) FETs(12,1) FETs(14,1)],[FETs(2,2) FETs(13,2) FETs(12,2) FETs(14,2)],'LineWidth',5,'DisplayName','M4')
plot([FETs(3,1) FETs(15,1) FETs(4,1)],[FETs(3,2) FETs(15,2) FETs(4,2)],'LineWidth',5,'DisplayName','M5')
plot([FETs(4,1) FETs(9,1) FETs(4,1) ],[FETs(4,2) FETs(9,2) FETs(4,2)],'LineWidth',5,'DisplayName','M6')
% Create legend
legend(axes1,'show');


%% Simulated vs Exp 


% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(Vg,[sim,exp],'MarkerSize',15,'LineWidth',4);
set(plot1(1),'DisplayName','Simulated','Marker','o','Color',[0 0 1]);
set(plot1(2),'DisplayName','Experimental','Marker','+','LineStyle','--',...
    'Color',[1 0 0]);

% Create ylabel
ylabel('Power Loss (W)');

% Create xlabel
xlabel('Input Voltage (V)');

% Create title
title('Experimental vs Simulated');

hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',24);
% Create legend
legend(axes1,'show');











