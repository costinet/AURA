
clear
tic
try
    
    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDB = obj;
    clear obj
    
    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductors.mat')
    indDB = obj;
    clear obj
    
indDBindex = [ 1   0   0   0   1   0   0   0   0   0   0   0   1   1   0   0   1   0   0   0 ]';

transDBindex = [ 0   1   1   1   1   0   0   0   0   1   1   1   0   0   0   0   0   1   1   1 ]';

% transDB = transistorDB;
% transDB.add(transDBig(logical(transDBindex)));
% 
% indDB = inductorDB;
% indDB.add(indDBig(logical(indDBindex)));


    Number_of_FETs = 20;
    FETs = zeros(Number_of_FETs,2);
    for i = 1:Number_of_FETs
        FETs(i,1) = transDB(i).ron.typ*1e-3;
        FETs(i,2) = transDB(i).Coss.typ*1e-12;
        FETs(i,3) = transDB(i).Width.approx*transDB(i).Length.approx;
    end
    
    for i = 1:Number_of_FETs
        FETs(end+1,1) = 0.5*transDB(i).ron.typ*1e-3;
        FETs(end,2) = 2*transDB(i).Coss.typ*1e-12;
        FETs(end,3) = 2*transDB(i).Width.approx*transDB(i).Length.approx;
    end
    
    for i = 1:Number_of_FETs
        FETs(end+1,1) = (1/3)*transDB(i).ron.typ*1e-3;
        FETs(end,2) = 3*transDB(i).Coss.typ*1e-12;
        FETs(end,3) = 3*transDB(i).Width.approx*transDB(i).Length.approx;
    end
    
        for i = 1:Number_of_FETs
        FETs(end+1,1) = (1/4)*transDB(i).ron.typ*1e-3;
        FETs(end,2) = 4*transDB(i).Coss.typ*1e-12;
        FETs(end,3) = 4*transDB(i).Width.approx*transDB(i).Length.approx;
    end
    
    
    for i = 1:length(FETs)
        for j = 1:length(FETs)
            FETs_Ron(i,j) = FETs(i,1) - FETs(j,1);
            FETs_Coss(i,j) = FETs(i,2) - FETs(j,2);
            FETs_Area(i,j) = FETs(i,3) - FETs(j,3);
        end
    end
    
   
    
    Number_of_INDs = 20;
    INDs = zeros(Number_of_INDs,2);
    
    for i = 1:Number_of_INDs
        INDs(i,1) = indDB(i).L.typ*1e-6;
        INDs(i,2) = indDB(i).rdc.typ*1e-3;
        INDs(i,3) = indDB(i).Width.approx*indDB(i).Length.approx;
    end
    
    for i = 1:Number_of_INDs
        INDs(end+1,1) = (1/2)*indDB(i).L.typ*1e-6;
        INDs(end,2) = (1/2)*indDB(i).rdc.typ*1e-3;
        INDs(end,3) = 2*indDB(i).Width.approx*indDB(i).Length.approx;
    end
    
    for i = 1:Number_of_INDs
        INDs(end+1,1) = (1/3)*indDB(i).L.typ*1e-6;
        INDs(end,2) = (1/3)*indDB(i).rdc.typ*1e-3;
        INDs(end,3) = 3*indDB(i).Width.approx*indDB(i).Length.approx;
    end
    
    for i = 1:Number_of_INDs
        INDs(end+1,1) = (1/4)*indDB(i).L.typ*1e-6;
        INDs(end,2) = (1/4)*indDB(i).rdc.typ*1e-3;
        INDs(end,3) = 4*indDB(i).Width.approx*indDB(i).Length.approx;
    end
    
    
    for i = 1:length(INDs)
        for j = 1:length(INDs)
            
            INDs_L(i,j) = INDs(i,1) - INDs(j,1);
            INDs_RL(i,j) = INDs(i,2) - INDs(j,2);
            INDs_Area(i,j) = INDs(i,3) - INDs(j,3);
        end
    end
    
    
    
    
    %{
    [minronCoss] = min(FETs,[],1);
    minron = minronCoss(1);
    minCoss = minronCoss(2);
    FETs(:,1) = FETs(:,1)./minron;
    FETs(:,2) = FETs(:,2)./minCoss;
    
    e_tan_dist = 1;
    e_FET_dist = 1;
    
    for i = 1:length(FETs)
        for j = 1:length(FETs)
            
            distance_between_points(i,j) = sqrt((FETs(j,1)-FETs(i,1))^2+(FETs(j,2)-FETs(i,2))^2);
            
        end
    end
    %}
    saved_fval = 100;
    big_loop = 0;
    max_iteration = 50;
    almost_steady_state = [0 0 0 0 0];
    saved_x = [];
    sf = [ 1 1 1 1e-6];

    Area_adjust = 1e5; % Used as a factor to add area and Ploss for the parato front fuction 
    
    % Best:
    % x = [3     3     3     3     3     3     2];
    % New expanded best : 7.0000    5.0000    6.0000    6.0000   21.0000    8.0000    3.4028
    % Latest: 9/21/22:
    %x =[1 2 3 4 5 4   1  ];
    
     x = [6 6 6 6 6     1     1];
    x = [7.0000    5.0000    6.0000    6.0000   21.0000    8.0000    3.4028];
    saved_x = x;
    FETs_number = 5;
    stick = zeros(FETs_number+1,3);
    deltaRon = -0.1e-3;
    deltaCoss = -50e-12;
    deltaL = -5e-9;
    deltaRL = -0.5e-3;
    
    while(big_loop<max_iteration)
        for i = 1:FETs_number
            
            adjust_FET = [0 0 0 0 0];
            adjust_FET(i) = 1;
            
            Coss_adj = [0 0 0 0 0];
            Ron_adj = [0 0 0 0 0];
            L_adj = [0];
            RL_adj = [0];
            
            [stick(i,1),graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
          
            Coss_adj = [0 0 0 0 0].*adjust_FET;
            Ron_adj = [deltaRon deltaRon deltaRon deltaRon deltaRon ].*adjust_FET;
            
            [stick(i,2),graph2]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            
            Coss_adj = [deltaCoss deltaCoss deltaCoss deltaCoss deltaCoss].*adjust_FET;
            Ron_adj = [0 0 0 0 0].*adjust_FET;
            
            [stick(i,3),graph3]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            % Move in the steepest direction:
            stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];
            
            D_Ron = FETs_Ron./deltaRon;
            D_Coss  = FETs_Coss./deltaCoss;
            
            
            Projected_Ploss = stick(i,1)+D_Ron(x(i),:)*stick_diff(1)+D_Coss(x(i),:)*stick_diff(2);
            
            Projected_PlossArea = Projected_Ploss + FETs(:,3)'.*Area_adjust;
            
          %  Projected_PlossArea = 

            [~,new_x] = min(Projected_PlossArea);
            
            x_test = x;
            x_test(i) = new_x;
            
            Coss_adj_test = [0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0];
            L_adj_test = [0];
            RL_adj_test = [0];
            
            [stick_em,graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);
           

            % If there is not an improvement to the fval of the converter.
           
            PAstick_test = stick_em + FETs(new_x,3).*Area_adjust;
            PAstick = stick(i,1) + FETs(x(i),3).*Area_adjust;
            if PAstick_test >= PAstick
                almost_steady_state(i) = 1;
            else
                x(i) = new_x;
                saved_x(end+1,:) = x;
                almost_steady_state(i) = 0;
                saved_fval(end+1) = stick(i,1);
                i
                toc
            end
           
        end
        
        
          
            i = 6;
            Coss_adj = [0 0 0 0 0];
            Ron_adj = [0 0 0 0 0];
            L_adj = [0];
            RL_adj = [0];
            
            [stick(i,1),graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            L_adj = [deltaL];
            RL_adj = [0];
            
            [stick(i,2),graph2]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            
            L_adj = [0];
            RL_adj = [deltaRL];
            
            [stick(i,3),graph3]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            % Move in the steepest direction:
            stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];
            
            D_Lon = INDs_L./deltaL;
            D_RLoss  = INDs_RL./deltaRL;
            
            
            Projected_Ploss = stick(i,1)+D_Lon(x(i),:)*stick_diff(1)+D_RLoss(x(i),:)*stick_diff(2);
            
            
            Projected_PlossArea = Projected_Ploss + INDs(:,3)'.*1e4;
            

            [~,new_x] = min(Projected_PlossArea);
           
            x_test = x;
            x_test(i) = new_x;
            
            Coss_adj_test = [0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0];
            L_adj_test = [0];
            RL_adj_test = [0];
            
            [stick_em,graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);
            
            % If there is not an improvement to the fval of the converter.
            PAstick_test = stick_em + INDs(new_x,3).*Area_adjust;
            PAstick = stick(i,1) + INDs(x(i),3).*Area_adjust;
            if PAstick_test >= PAstick
                almost_steady_state(i) = 1;
            else
                x(i) = new_x;
                saved_x(end+1,:) = x;
                almost_steady_state(i) = 0;
                saved_fval(end+1) = stick(i,1);
                i
                toc
            end

            fs_pos = 7;
            fs_sweep = linspace(5e6,1e6,10);
            
            if big_loop>0
                fs_sweep = linspace(500e3+x(7)*1e6,-500e3+x(7)*1e6,9);
            end
            for n=1:length(fs_sweep)
                
                Coss_adj_test = [0 0 0 0 0];
                Ron_adj_test = [0 0 0 0 0];
                L_adj_test = [0];
                RL_adj_test = [0];
                
                x(fs_pos) = fs_sweep(n)*1e-6;
                [stick_fs(n),graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead(x,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);
                
                
            end
            [~,pos] = min(stick_fs);
            if fs_sweep(pos)*1e-6 ~= saved_x(end,fs_pos)
            x(fs_pos) = fs_sweep(pos)*1e-6;
            saved_x(end+1,:) = x;
            almost_steady_state(fs_pos) = 0;
            saved_fval(end+1) = stick_fs(pos);
            else
                x(fs_pos) = fs_sweep(pos)*1e-6;
                almost_steady_state(fs_pos) = 1;
            end
            
            
            
            big_loop = big_loop + 1;
            
            
            if sum(almost_steady_state)==length(almost_steady_state) && big_loop>1
                toc
                break
                J = 4564654;
            end
    end
    
    
    
    
catch ME
    J = 5465456465;
    
end

J = 456456456;
figure(9)
plot_ron = FETs(:,1);
plot_Coss = FETs(:,2);
plot_FArea = FETs(:,3);
scatter3(plot_ron,plot_Coss,plot_FArea,'DisplayName','FETs');
hold on
plot3(plot_ron(saved_x(:,1)),plot_Coss(saved_x(:,1)),plot_FArea(saved_x(:,1)),'DisplayName','M1,4,8,9','LineWidth',3,'LineStyle','--');
hold on
plot3(plot_ron(saved_x(:,2)),plot_Coss(saved_x(:,2)),plot_FArea(saved_x(:,2)),'DisplayName','M2,3','LineWidth',3,'LineStyle',':');
hold on
plot3(plot_ron(saved_x(:,3)),plot_Coss(saved_x(:,3)),plot_FArea(saved_x(:,3)),'DisplayName','M5,10,13','LineWidth',3,'LineStyle','-.');
hold on
plot3(plot_ron(saved_x(:,4)),plot_Coss(saved_x(:,4)),plot_FArea(saved_x(:,4)),'DisplayName','M6,11,14','LineWidth',3,'LineStyle','--');
hold on
plot3(plot_ron(saved_x(:,5)),plot_Coss(saved_x(:,5)),plot_FArea(saved_x(:,5)),'DisplayName','M7,12,15','LineWidth',3,'LineStyle',':');
hold on
scatter3(plot_ron(x(:,1:5)),plot_Coss(x(:,1:5)),plot_FArea(x(:,1:5)),50,'r','DisplayName','Ending Point');

% Create ylabel
ylabel('Coss (F)','FontSize',20);

% Create xlabel
xlabel('Ron (\Omega)','FontSize',20);

% Create zlabel
zlabel('Area (m^2)','FontSize',20);

% Create title
title('SC BuckHybridBuck','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend


figure(22)
plot_L = INDs(:,1);
plot_RL = INDs(:,2);
plot_LArea = INDs(:,3);

scatter3(plot_L,plot_RL,plot_LArea,'DisplayName','Inductors');
hold on
plot3(plot_L(saved_x(:,6)),plot_RL(saved_x(:,6)),plot_LArea(saved_x(:,6)),'DisplayName','L1,2,3','LineWidth',3,'LineStyle','--');
hold on
scatter3(plot_L(x(:,6)),plot_RL(x(:,6)),plot_LArea(x(:,6)),50,'r','DisplayName','Ending Point');

% Create ylabel
ylabel('RL (\Omega)','FontSize',20);

% Create xlabel
xlabel('L (H)','FontSize',20);

% Create zlabel
zlabel('Area (m^2)','FontSize',20);

% Create title
title('SC BuckHybridBuck','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend


figure(24)
scatter(saved_x(2:end,7),saved_fval(2:end),'DisplayName','Test Points');
hold on
plot(saved_x(2:end,7),saved_fval(2:end),'DisplayName','fs','LineWidth',3,'LineStyle','--');
hold on
scatter(saved_x(end,7),saved_fval(end),50,'r','DisplayName','Ending Point');


% Create ylabel
ylabel('Pout (W)','FontSize',20);

% Create xlabel
xlabel('fs (MHz)','FontSize',20);

% Create title
title('SC BuckHybridBuck','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend



figure(30)
total_area = zeros(size(saved_x,1),1);
for i = 1:size(saved_x,1)
   for j = 1:size(saved_x,2)
       if j<=5
       total_area(i) = total_area(i) + FETs(saved_x(i,j),3);
       end
       if j==6
          total_area(i) = total_area(i) + INDs(saved_x(i,j),3);
       end
       if j>6
           continue
       end
   end
end

scatter(saved_x(2:end,7),total_area(2:end),'DisplayName','Test Points');
hold on
plot(saved_x(2:end,7),total_area(2:end),'DisplayName','fs','LineWidth',3,'LineStyle','--');
hold on
scatter(saved_x(end,7),total_area(end),50,'r','DisplayName','Ending Point');


% Create ylabel
ylabel('Area (m^2)','FontSize',20);

% Create xlabel
xlabel('fs (MHz)','FontSize',20);

% Create title
title('SC BuckHybridBuck','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend


figure(48)

scatter(saved_fval(2:end),total_area(2:end),'DisplayName','Test Points');
hold on
plot(saved_fval(2:end),total_area(2:end),'DisplayName','fs','LineWidth',3,'LineStyle','--');
hold on
scatter(saved_fval(2:end),total_area(end),50,'r','DisplayName','Ending Point');


% Create ylabel
ylabel('Area (m^2)','FontSize',20);

% Create xlabel
xlabel('Pout (W)','FontSize',20);

% Create title
title('SC BuckHybridBuck','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend

