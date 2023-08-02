% This code completes an optimization on a buck converter transistor,
% inductor and frequency. Note: the load files will need to be changed to
% the directory the database files are located

clear
tic
try

    % Load in FETs and Inductors
    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDB = obj;
    clear obj

    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductors.mat')
    indDB = obj;
    clear obj

    % Calculate Ron and Coss Differences
    Number_of_FETs = 20;
    FETs = zeros(Number_of_FETs,2);

    Number_of_FETs = 16;
    FETs = zeros(Number_of_FETs,2);
    for i = 1:Number_of_FETs
        FETs(i,1) = transDB(i).ron.typ*1e-3;
        FETs(i,2) = transDB(i).Coss.typ*1e-12;
        FETs(i,3) = transDB(i).Width.approx*transDB(i).Length.approx;
    end

    for i = 1:length(FETs)
        for j = 1:length(FETs)
            FETs_Ron(i,j) = FETs(i,1) - FETs(j,1);
            FETs_Coss(i,j) = FETs(i,2) - FETs(j,2);
            FETs_Area(i,j) = FETs(i,3) - FETs(j,3);
        end
    end


  % Calculate RL and L Differences
    Number_of_INDs = 20;
    INDs = zeros(Number_of_INDs,2);

    for i = 1:Number_of_INDs
        INDs(i,1) = indDB(i).L.typ*1e-6;
        INDs(i,2) = indDB(i).rdc.typ*1e-3;
        INDs(i,3) = indDB(i).Width.approx*indDB(i).Length.approx;
    end

    for i = 1:length(INDs)
        for j = 1:length(INDs)

            INDs_L(i,j) = INDs(i,1) - INDs(j,1);
            INDs_RL(i,j) = INDs(i,2) - INDs(j,2);
            INDs_Area(i,j) = INDs(i,3) - INDs(j,3);
        end
    end

    
    saved_fval = 100;
    big_loop = 0;
    max_iteration = 50;
    almost_steady_state = [0 0 0 0];
    saved_x = [];
    sf = [ 1 1 1 1e-6];

    % Set inital conditions
    x =[3.0000   3.0000  3.0000   1 ];

    FETs_number = 2;
    stick = zeros(FETs_number+1,3);
    % Set Pertubations
    deltaRon = -0.1e-3;
    deltaCoss = -50e-12;
    deltaL = -5e-9;
    deltaRL = -0.5e-3;
    saved_x(end+1,:) = x;

    while(big_loop<max_iteration)
        for i = 1:FETs_number
            
            % Set adjustment values for pertubations
            adjust_FET = [0 0];
            adjust_FET(i) = 1;

            Coss_adj = [0 0];
            Ron_adj = [0 0];
            L_adj = [0];
            RL_adj = [0];

            % Normal Value
            [stick(i,1)]=COMPEL_2023_Buck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);

            Coss_adj = [0 0 ].*adjust_FET;
            Ron_adj = [deltaRon deltaRon ].*adjust_FET;
            
            % Perturb Ron 
            [stick(i,2)]=COMPEL_2023_Buck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);


            Coss_adj = [deltaCoss deltaCoss].*adjust_FET;
            Ron_adj = [0 0 ].*adjust_FET;

            % Perturb Coss 
            [stick(i,3)]=COMPEL_2023_Buck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);

            % Move in the steepest direction:
            stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];


            D_Ron = FETs_Ron./deltaRon;
            D_Coss  = FETs_Coss./deltaCoss;

            % Project power loss to other transistors
            Projected_Ploss = stick(i,1)+D_Ron(x(i),:)*stick_diff(1)+D_Coss(x(i),:)*stick_diff(2);

            % Find minimum power loss
            [~,new_x] = min(Projected_Ploss);

            x_test = x;
            x_test(i) = new_x;

            Coss_adj_test = [0 0];
            Ron_adj_test = [0 0];
            L_adj_test = [0];
            RL_adj_test = [0];

            % Test new transistor
            [stick_em]=COMPEL_2023_Buck_D_Sweep_dead(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);

            % If there is not an improvement to the fval of the converter then assume optimum point is reached else update and move on.
            if stick_em >= stick(i,1)
                almost_steady_state(i) = 1;
            else
                x(i) = new_x;
                saved_x(end+1,:) = x;
                almost_steady_state(i) = 0;
                saved_fval(:,end+1) = stick(i,1);
                i
                toc
            end

        end



        i = 3;
        Coss_adj = [0 0];
        Ron_adj = [0 0];
        L_adj = [0];
        RL_adj = [0];

        [stick(i,1)]=COMPEL_2023_Buck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);

        L_adj = [deltaL];
        RL_adj = [0];

        [stick(i,2)]=COMPEL_2023_Buck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);


        L_adj = [0];
        RL_adj = [deltaRL];

        [stick(i,3)]=COMPEL_2023_Buck_D_Sweep_dead(x,Coss_adj,Ron_adj,L_adj,RL_adj);

        % Move in the steepest direction:
        stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];

        D_Lon = INDs_L./deltaL;
        D_RLoss  = INDs_RL./deltaRL;

         
        Projected_Ploss = stick(i,1)+D_Lon(x(i),:)*stick_diff(1)+D_RLoss(x(i),:)*stick_diff(2);


        [~,new_x] = min(Projected_Ploss);

        x_test = x;
        x_test(i) = new_x;

        Coss_adj_test = [0 0];
        Ron_adj_test = [0 0];
        L_adj_test = [0];
        RL_adj_test = [0];

        [stick_em]=COMPEL_2023_Buck_D_Sweep_dead(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);

        % If there is not an improvement to the fval of the converter.
        if stick_em >= stick(i,1)
            almost_steady_state(i) = 1;
        else
            x(i) = new_x;
            saved_x(end+1,:) = x;
            almost_steady_state(i) = 0;
            saved_fval(:,end+1) = stick(i,1);
            i
            toc
        end

        fs_pos = 4;
        fs_sweep = linspace(5e6,1e6,10);

        if big_loop>0
            fs_sweep = linspace(500e3+x(fs_pos)*1e6,-500e3+x(fs_pos)*1e6,9);
        end

        if big_loop>0 && x(fs_pos) < 0.550

            fs_sweep = linspace((x(fs_pos)*1e6-100e3)+x(fs_pos)*1e6,100e3,9);
        end

        for n=1:length(fs_sweep)

            Coss_adj_test = [0 0 ];
            Ron_adj_test = [0 0 ];
            L_adj_test = [0];
            RL_adj_test = [0];

            x(fs_pos) = fs_sweep(n)*1e-6;
            [stick_fs(n)]=COMPEL_2023_Buck_D_Sweep_dead(x,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);

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
plot(plot_ron(saved_x(:,1)),plot_Coss(saved_x(:,1)),'DisplayName','M1','LineWidth',3,'LineStyle','--');
hold on
plot(plot_ron(saved_x(:,2)),plot_Coss(saved_x(:,2)),'DisplayName','M2','LineWidth',3,'LineStyle',':');
hold on
scatter(plot_ron(x(:,1:2)),plot_Coss(x(:,1:2)),50,'r','DisplayName','Ending Point');

% Create ylabel
ylabel('Coss (F)','FontSize',20);

% Create xlabel
xlabel('Ron (\Omega)','FontSize',20);

% Create title
title('Buck Converter Optimization','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend

figure
plot_L = INDs(:,1);
plot_RL = INDs(:,2);
scatter(plot_L,plot_RL,'DisplayName','Inductors');
hold on
plot(plot_L(saved_x(:,3)),plot_RL(saved_x(:,3)),'DisplayName','L1','LineWidth',3,'LineStyle','--');
hold on
scatter(plot_L(x(:,3)),plot_RL(x(:,3)),50,'r','DisplayName','Ending Point');

% Create ylabel
ylabel('RL (\Omega)','FontSize',20);

% Create xlabel
xlabel('L (H)','FontSize',20);

% Create title
title('Buck Converter Optimization','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend


figure
scatter(saved_x(2:end,4),saved_fval(2:end),'DisplayName','Test Points');
hold on
plot(saved_x(2:end,4),saved_fval(2:end),'DisplayName','fs','LineWidth',3,'LineStyle','--');
hold on
scatter(saved_x(end,4),saved_fval(end),50,'r','DisplayName','Ending Point');
% Create ylabel
ylabel('Pout (W)','FontSize',20);

% Create xlabel
xlabel('fs (MHz)','FontSize',20);

% Create title
title('Buck Converter Optimization','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend