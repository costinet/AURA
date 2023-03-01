
clear
tic
try
    
    
    
    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDBig = obj;
    clear obj
    
    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductors.mat')
    indDBig = obj;
    clear obj
    
indDBindex = [    
   1
   0
   0
   0
   1
   0
   0
   0
   0
   0
   0
   0
   1
   1
   0
   0
   1
   0
   0
   0];


transDBindex = [    
   0
   0
   0
   0
   0
   0
   0
   0
   0
   1
   1
   1
   0
   0
   0
   0
   0
   1
   1
   1];
    

transDB = transistorDB;
transDB.add(transDBig(logical(transDBindex)));

indDB = inductorDB;
indDB.add(indDBig(logical(indDBindex)));

    Number_of_FETs = 6;
    FETs = zeros(Number_of_FETs,2);
    for i = 1:Number_of_FETs
        FETs(i,1) = transDB(i).ron.typ*1e-3;
        FETs(i,2) = transDB(i).Coss.typ*1e-12;
    end
    
    for i = 1:length(FETs)
        for j = 1:length(FETs)
            FETs_Ron(i,j) = FETs(i,1) - FETs(j,1);
            FETs_Coss(i,j) = FETs(i,2) - FETs(j,2);
        end
    end
    
   

    Number_of_INDs = 5;
    INDs = zeros(Number_of_INDs,2);
    
    for i = 1:Number_of_INDs
        INDs(i,1) = indDB(i).L.typ*1e-6;
        INDs(i,2) = indDB(i).rdc.typ*1e-3;
    end
    
    for i = 1:length(INDs)
        for j = 1:length(INDs)
            
            INDs_L(i,j) = INDs(i,1) - INDs(j,1);
            INDs_RL(i,j) = INDs(i,2) - INDs(j,2);
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
    saved_fval = ones(6,1).*100;
    big_loop = 0;
    max_iteration = 50;
    almost_steady_state = [0 0 0];
    saved_x = [];
    sf = [ 1 1 1 1e-6];
    
    % Best:
     x = [3     3     3     3     3     3     2];
    % Latest: 9/21/22:
    %x =[1 2 3 4 5 4   1  ];
    
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
            
            [stick(i,1),graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            Coss_adj = [0 0 0 0 0].*adjust_FET;
            Ron_adj = [deltaRon deltaRon deltaRon deltaRon deltaRon ].*adjust_FET;
            
            [stick(i,2),graph2]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            
            Coss_adj = [deltaCoss deltaCoss deltaCoss deltaCoss deltaCoss].*adjust_FET;
            Ron_adj = [0 0 0 0 0].*adjust_FET;
            
            [stick(i,3),graph3]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            % Move in the steepest direction:
            stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];
            
            D_Ron = FETs_Ron./deltaRon;
            D_Coss  = FETs_Coss./deltaCoss;
            
            
            Projected_Ploss = stick(i,1)+D_Ron(x(i),:)*stick_diff(1)+D_Coss(x(i),:)*stick_diff(2);
            
            %{
    if  the_sign(1) == 1 && the_sign(1) == 1
        
        errorRon = FETs(:,1)>=FETs(x(i),1);
        errorCoss = FETs(:,2)>=FETs(x(i),2);
    end
    
    if I == 1 && the_sign == -1
        
        errorRon = FETs(:,1)<=FETs(x(i),1);
        errorCoss = FETs(:,2)<=FETs(x(i),2);
    end
    
    if I == 2 && the_sign == 1
        
        errorRon = FETs(:,1)>=FETs(x(i),1);
       
    end
    
    if I == 2 && the_sign == -1
        
        errorRon = FETs(:,1)<=FETs(x(i),1);
        
    end
    
     if I == 3 && the_sign == 1
        

        errorCoss = FETs(:,2)>=FETs(x(i),2);
    end
    
    if I == 3 && the_sign == -1
        

        errorCoss = FETs(:,2)<=FETs(x(i),2);
    end
            %}
            
            [~,new_x] = min(Projected_Ploss);
            
            x_test = x;
            x_test(i) = new_x;
            
            Coss_adj_test = [0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0];
            L_adj_test = [0];
            RL_adj_test = [0];
            
            [stick_em,graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);
            
            % If there is not an improvement to the fval of the converter.
            if stick_em >= stick(i,1)
                almost_steady_state(i) = 1;
            else
                x(i) = new_x;
                saved_x(end+1,:) = x;
                almost_steady_state(i) = 0;
                saved_fval(:,end+1) = stick(:,1);
                i
                toc
            end
           
        end
        
        
          
            i = 6;
            Coss_adj = [0 0 0 0 0];
            Ron_adj = [0 0 0 0 0];
            L_adj = [0];
            RL_adj = [0];
            
            [stick(i,1),graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            L_adj = [deltaL];
            RL_adj = [0];
            
            [stick(i,2),graph2]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            
            L_adj = [0];
            RL_adj = [deltaRL];
            
            [stick(i,3),graph3]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x,Coss_adj,Ron_adj,L_adj,RL_adj);
            
            % Move in the steepest direction:
            stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];
            
            D_Lon = INDs_L./deltaL;
            D_RLoss  = INDs_RL./deltaRL;
            
            
            Projected_Ploss = stick(i,1)+D_Lon(x(i),:)*stick_diff(1)+D_RLoss(x(i),:)*stick_diff(2);
            
            
   
            
            [~,new_x] = min(Projected_Ploss);
           
            x_test = x;
            x_test(i) = new_x;
            
            Coss_adj_test = [0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0];
            L_adj_test = [0];
            RL_adj_test = [0];
            
            [stick_em,graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);
            
            % If there is not an improvement to the fval of the converter.
            if stick_em >= stick(i,1)
                almost_steady_state(i) = 1;
            else
                x(i) = new_x;
                saved_x(end+1,:) = x;
                almost_steady_state(i) = 0;
                saved_fval(:,end+1) = stick(:,1);
                i
                toc
            end

            
            fs_sweep = linspace(2e6,500e3,10);
            for n=1:length(fs_sweep)
                
                Coss_adj_test = [0 0 0 0 0];
                Ron_adj_test = [0 0 0 0 0];
                L_adj_test = [0];
                RL_adj_test = [0];
                
                x(7) = fs_sweep(n)*1e-6;
                [stick_fs(n),graph1]=COMPEL_2023_SCBuckHybridBuck_D_Sweep(x,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);
                
                
            end
            [~,pos] = min(stick_fs);
            x(7) = fs_sweep(pos)*1e-6;
            big_loop = big_loop + 1;
            
            
            if sum(almost_steady_state)==length(almost_steady_state)
                toc
                break
                J = 4564654;
            end
    end
    
    
    
    
catch ME
    J = 5465456465;
    
end

J = 456456456;
figure
plot_ron = FETs(:,1);
plot_Coss = FETs(:,2);
scatter(plot_ron,plot_Coss,'DisplayName','FETs');
hold on
plot(plot_ron(saved_x(:,1)),plot_Coss(saved_x(:,1)),'DisplayName','M1,4,8,9','LineWidth',3,'LineStyle','--');
hold on
plot(plot_ron(saved_x(:,2)),plot_Coss(saved_x(:,2)),'DisplayName','M2,3','LineWidth',3,'LineStyle',':');
hold on
plot(plot_ron(saved_x(:,3)),plot_Coss(saved_x(:,3)),'DisplayName','M5,10,13','LineWidth',3,'LineStyle','-.');
hold on
plot(plot_ron(saved_x(:,4)),plot_Coss(saved_x(:,4)),'DisplayName','M6,11,14','LineWidth',3,'LineStyle','--');
hold on
plot(plot_ron(saved_x(:,5)),plot_Coss(saved_x(:,5)),'DisplayName','M7,12,15','LineWidth',3,'LineStyle',':');
hold on
scatter(plot_ron(x(:,1:5)),plot_Coss(x(:,1:5)),50,'r','DisplayName','Ending Point');

% Create ylabel
ylabel('Coss (F)','FontSize',20);

% Create xlabel
xlabel('Ron (\Omega)','FontSize',20);

% Create title
title('SC BuckHybridBuck','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend
