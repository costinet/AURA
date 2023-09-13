function [stick] = COMPEL_2023_Optim_SCBuckHybridBuck_Buck_D_Final_Multiple_Loops(X)
% This function outputs to the GA. It completes the optimization of FETs,
% Inductors, and frequency.
% The input variable X is the input from from the GA.
% It should be in the format

% X = [# of FETs in Parallel, # of Inductors in Parallel, # of Capacitors in Parallel, Iout]

% stick is the output to the GA it should have the form.
% stick = [Area in mm/Iout,  1-efficiency]
% stick will be used as the parato front variable

tic
try

    debug = 1;


    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDB = obj;
    clear obj

    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductorsFull.mat')
    indDB = obj;
    clear obj

    load('D:\GitHub\AURA\AURAdb\@topologyDB\topology.mat')
    topDB = obj;
    clear obj

    load('D:\GitHub\AURA\AURAdb\databases\@capacitorDB\capacitors.mat')
    capDB = obj;
    clear obj



    % indDBindex = [ 1   0   0   0   1   0   0   0   0   0   0   0   1   1   0   0   1   0   0   0 ]';

    % transDBindex = [ 0   1   1   1   1   0   0   0   0   1   1   1   0   0   0   0   0   1   1   1 ]';

    % transDB = transistorDB;
    % transDB.add(transDBig(logical(transDBindex)));
    %
    % indDB = inductorDB;
    % indDB.add(indDBig(logical(indDBindex)));


    %% Extract parameters from databases

    FETs_number = 5;
    FETs_in_parallel = X(1:FETs_number);
    Number_of_FETs = 20;
    FETs = zeros(Number_of_FETs,2,FETs_number);

    for j = 1:FETs_number
        for i = 1:Number_of_FETs
            FETs(i,1,j) = (1/FETs_in_parallel(j))*transDB(i).ron.typ*1e-3;
            FETs(i,2,j) = FETs_in_parallel(j)*transDB(i).Coss.typ*1e-12;
            %FETs(i,3,j) =
            %FETs_in_parallel(j)*transDB(i).Width.approx*transDB(i).Length.approx;
            % Used if I am looking a changing area here as well
        end
    end
    %{
    for k = 1:FETs_number
        for i = 1:length(FETs)
            for j = 1:length(FETs)
                FETs_Ron(i,j,k) = FETs(i,1,k) - FETs(j,1,k);
                FETs_Coss(i,j,k) = FETs(i,2,k) - FETs(j,2,k);
                % FETs_Area(i,j) = FETs(i,3) - FETs(j,3); %Used if I am looking a changing area here as well
            end
        end
    end
    %}

    INDs_number = 1;
    INDs_in_parallel = X(FETs_number+1:FETs_number+INDs_number);
    Number_of_INDs = 36;
    INDs = zeros(Number_of_INDs,2);

    for j = 1:INDs_number
        for i = 1:Number_of_INDs
            INDs(i,1,j) = (1/INDs_in_parallel(j))*indDB(i).L.typ*1e-6;
            INDs(i,2,j) = (1/INDs_in_parallel(j))*indDB(i).rdc.typ*1e-3;
            %INDs(i,3) = 2*indDB(i).Width.approx*indDB(i).Length.approx;
            % Used if I am looking a changing area here as well
        end
    end
    %{
    for k = 1:INDs_number
        for i = 1:length(INDs)
            for j = 1:length(INDs)
                INDs_L(i,j,k) = INDs(i,1,k) - INDs(j,1,k);
                INDs_RL(i,j,k) = INDs(i,2,k) - INDs(j,2,k);
                %INDs_Area(i,j) = INDs(i,3) - INDs(j,3);
                % Used if I am looking a changing area here as well
            end
        end
    end
    %}


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

    %% Set up inital condidtions 
    saved_fval = 100;
    big_loop = 0;
    max_iteration = 50;
    almost_steady_state = [0 0 0 0 0];
    saved_x = [];
    sf = [ones(1,5), 1, 1e-6] ;


    %                  |  |
    %x = [8 10 1 1 1 11 0.5 X];

    x = [9.0000    9.0000   10.0000    1.0000   4  6.0000   0.5 X];

    saved_x = x;
    FETs_number = 5;
    stick = zeros(FETs_number+1,3);
    deltaRon = -0.1e-3;
    deltaCoss = -50e-12;
    deltaL = -5e-9;
    deltaRL = -0.5e-3;
    frequencyIndex = 7;
    deltaF = -1e3*1e-6;
    scaleFETs = [1e2 1e9 1e6];
    scaleINDs = [1e6 1e3 1e6];

    %% Debug map showing direction of greatest descent on each iteration 
    if debug
        for i = 1:FETs_number
            figure(100+i)
            scatter(FETs(:,1,i).*scaleFETs(1),FETs(:,2,i).*scaleFETs(2),'k','filled')
            hold on
            h = voronoi(FETs(:,1,i).*scaleFETs(1),FETs(:,2,i).*scaleFETs(2));
            for j = 1:length(h)
                h(j).LineWidth = 2;
            end
           % labelpoints(FETs(:,1,i).*scaleFETs(1),FETs(:,2,i).*scaleFETs(2),[1:Number_of_FETs],'E',0.1)
            xlabel('Ron (scaled)')
            ylabel('Coss (scaled)')
            findfigs
        end
        for i = 1:INDs_number
            figure(100+FETs_number+i)
            scatter(INDs(:,1).*scaleINDs(1),INDs(:,2,i).*scaleINDs(2),'k','filled')
            hold on

            h = voronoi(INDs(:,1).*scaleINDs(1),INDs(:,2,i).*scaleINDs(2));
            for j = 1:length(h)
                h(j).LineWidth = 2;
            end
            labelpoints(INDs(:,1).*scaleINDs(1),INDs(:,2,i).*scaleINDs(2),[1:Number_of_INDs],'E',0.1)
            xlabel('L (scaled)')
            ylabel('RL (scaled)')
            findfigs
        end

    end
    
    %% Set inital values
    Coss_adj = [0 0 0 0 0];
    Ron_adj = [0 0 0 0 0];
    L_adj = [0];
    RL_adj = [0];

    [Inital_Eff]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);
    saved_fval = (1-Inital_Eff(2));

%% Optimization Loop
    while(big_loop<max_iteration)

        %% FETs Loop
        for i = 1:FETs_number
            
            NoFETImprovFLAG = 0;
            adjust_FET = [0 0 0 0 0];
            adjust_FET(i) = 1;

            Coss_adj = [0 0 0 0 0];
            Ron_adj = [0 0 0 0 0];
            L_adj = [0];
            RL_adj = [0];

            [Control_Eff]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);

            tempx = x;
            x(frequencyIndex) = x(frequencyIndex)-deltaF;
            [PerturbF]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);

            x = tempx;


            Coss_adj = [0 0 0 0 0].*adjust_FET;
            Ron_adj = [deltaRon deltaRon deltaRon deltaRon deltaRon ].*adjust_FET;

            [PerturbRon]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);


            Coss_adj = [deltaCoss deltaCoss deltaCoss deltaCoss deltaCoss].*adjust_FET;
            Ron_adj = [0 0 0 0 0].*adjust_FET;

            [PerturbCoss]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);


            ic = x(i);
            oldpt = ic;
            ts = 1/(x(frequencyIndex)/sf(frequencyIndex))*scaleFETs(3);
            data = FETs(:,:,i);
            data = data .* scaleFETs(1:2);
            ptdist = data-data(ic,:);

            dETAdRON  = ((1-PerturbRon(2)) - (1-Control_Eff(2))) / (deltaRon*scaleFETs(1));
            dETAdCOSS  = ((1-PerturbCoss(2)) - (1-Control_Eff(2))) / (deltaCoss*scaleFETs(2));
            dETAdF  = ((1-PerturbF(2)) - (1-Control_Eff(2))) / (deltaF*scaleFETs(3));
            minvec = [dETAdRON dETAdCOSS dETAdF];

            mindist = sum(ptdist.^2,2) ./ sum(ptdist.*minvec(1:2),2);
            mindist(mindist < 0) = inf;
            [deltaDist,nextpt] = mink(mindist,2);

            if nextpt(1) == oldpt
                nextpt = nextpt(2);
                deltaDist = deltaDist(2);
            else
                nextpt = nextpt(1);
                deltaDist = deltaDist(1);
            end
            %%{
            if isinf(mindist(nextpt))
                dT = (1+.05*minvec(3)*10);
                dT = min(max(dT,.5),1.5);
                ts = ts*dT;
                NoFETImprovFLAG = 1;
            end
            dT = (1+deltaDist*minvec(3)*10);
            dT = min(max(dT,.5),1.5);
            ts = ts*dT;
            %%}


            if debug
                figure(100+i)
                hold on
                quiver(data(ic,1),data(ic,2),5*minvec(1),5*minvec(2),'r','LineWidth',3)
            end



            x_test = x;
            x_test(i) = nextpt;
            x_test(frequencyIndex) = 1/(ts/scaleFETs(3))*sf(frequencyIndex);
            Coss_adj_test = [0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0];
            L_adj_test = [0];
            RL_adj_test = [0];

            [Tested_Eff]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);

            if (Tested_Eff(2)>Control_Eff(2))||NoFETImprovFLAG
                % If there is not an improvement to the fval of the converter.
                almost_steady_state(i) = 1;
            else
                x = x_test;
                saved_x(end+1,:) = x;
                almost_steady_state(i) = 0;
                saved_fval(end+1) = Tested_Eff(2);
                i
                toc
            end
        end


        %% INDs Loop

        for i = FETs_number+1:FETs_number+INDs_number
            NoINDImprovFLAG = 0;
            Coss_adj = [0 0 0 0 0];
            Ron_adj = [0 0 0 0 0];
            L_adj = [0];
            RL_adj = [0];

            [Control_Eff]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);


            tempx = x;
            x(frequencyIndex) = x(frequencyIndex)-deltaF;
            [PerturbF]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);

            x = tempx;

            L_adj = [deltaL];
            RL_adj = [0];

            [PerturbL]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);


            L_adj = [0];
            RL_adj = [deltaRL];

            [PerturbRL]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);



            ic = x(i);
            oldpt = ic;
            ts = 1/(x(frequencyIndex)/sf(frequencyIndex))*scaleINDs(3);
            data = INDs;
            %data(:,1) = 1./(INDs(:,1).^2); % Uses 1/L2 as the factor for inductance
            data = data .* scaleINDs(1:2);
            ptdist = data-data(ic,:);

            dETAdL  = ((1-PerturbL(2)) - (1-Control_Eff(2))) / (deltaL*scaleINDs(1));
            dETAdRL  = ((1-PerturbRL(2)) - (1-Control_Eff(2))) / (deltaRL*scaleINDs(2));
            dETAdF  = ((1-PerturbF(2)) - (1-Control_Eff(2))) / (deltaF*scaleINDs(3));
            minvec = [dETAdL dETAdRL dETAdF];

            mindist = sum(ptdist.^2,2) ./ sum(ptdist.*minvec(1:2),2);
            mindist(mindist < 0) = inf;
            [deltaDist,nextpt] = mink(mindist,2);

            if nextpt(1) == oldpt
                nextpt = nextpt(2);
                deltaDist = deltaDist(2);
            else
                nextpt = nextpt(1);
                deltaDist = deltaDist(1);
            end
            %%{
            if isinf(mindist(nextpt))
                dT = (1+.05*minvec(3)*10);
                dT = min(max(dT,.5),1.5);
                ts = ts*dT;
                NoFETImprovFLAG = 1;
            end
            dT = (1+deltaDist*minvec(3)*10);
            dT = min(max(dT,.5),1.5);
            ts = ts*dT;
            %%}


            if debug
                figure(100+i)
                hold on
                quiver(data(ic,1),data(ic,2),5*minvec(1),5*minvec(2),'r','LineWidth',3)
            end



            x_test = x;
            x_test(i) = nextpt;
            x_test(frequencyIndex) = 1/(ts/scaleINDs(3))*sf(frequencyIndex);
            Coss_adj_test = [0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0];
            L_adj_test = [0];
            RL_adj_test = [0];

            [Tested_Eff]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x_test,Coss_adj_test,Ron_adj_test,L_adj_test,RL_adj_test);

            if (Tested_Eff(2)>Control_Eff(2))||NoFETImprovFLAG
                % If there is not an improvement to the fval of the converter.
                almost_steady_state(i) = 1;
            else
                x = x_test;
                saved_x(end+1,:) = x;
                almost_steady_state(i) = 0;
                saved_fval(end+1) = Tested_Eff(2);
                i
                toc
            end
        end

        if sum(almost_steady_state)==length(almost_steady_state)
            toc
            break
            J = 4564654;
        end

    end


    %% Set final values
    Coss_adj = [0 0 0 0 0];
    Ron_adj = [0 0 0 0 0];
    L_adj = [0];
    RL_adj = [0];

    [stick]=COMPEL_2023_SCBuckHybridBuck_D_Sweep_dead_Final(x,Coss_adj,Ron_adj,L_adj,RL_adj);

    return
catch ME
    J = 5465456465;
end

%{
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

%}