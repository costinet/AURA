
clear
tic
try
    Number_of_FETs = 100;
    FETs = zeros(Number_of_FETs,2);
    for i = 1:Number_of_FETs
        [FETs(i,1),FETs(i,2)]=Select_IC_FET(i);
    end
    
    for i = 1:length(FETs)
        for j = 1:length(FETs)
            FETs_Ron(i,j) = FETs(i,1) - FETs(j,1);
            FETs_Coss(i,j) = FETs(i,2) - FETs(j,2);
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
    almost_steady_state = [0 0 0 0 0 0 0 0];
    
    sf = [1 1 1 1 1 1e-6 100 100];
    
    %  0.1389    0.0793    0.2542    0.3055    0.0271    0.0763    0.0424    0.1863
    
    x = [ 34    19    62    76    05   19    11    47   1.2951   0.009 0.009];
    
    %  x = [0.1389    0.0793    0.2542    0.3055    0.0271     1.2951   0.009 0.009];
    
    %  x = [500   500   900    900    100     1.2951   0.009 0.009];
    
    % x = [535+13 438-13 1000 1000 1 1.2951   0.009 0.009];
    
    saved_x = x;
    FETs_number = 8;
    stick = zeros(FETs_number,FETs_number);
    deltaRon = -1e-3;
    deltaCoss = -0.01e-9;
    deltaW = 1;
    % coss_adjust_values = [0 deltaCoss/minCoss];
    %  ron_adjust_values = [deltaRon/minron 0];
    Coss_adj = [0 0 0 0 0 0 0 0];
    Ron_adj = [0 0 0 0 0 0 0 0];
    
    while(big_loop<max_iteration)
        for i = 1:FETs_number
            if x(i) == 100
                continue
            end
            for j = 1:FETs_number
                if x(j) == 1
                    continue
                end
                if i == j
                    [stick(i,j),graph1]=Texas_SC_Fib_4_IC_Sweep(x,Coss_adj,Ron_adj);
                    
                else
                    x_test = x;
                    x_test(i) = x(i)+1;
                    x_test(j) = x(j)-1;
                    [stick(i,j),graph1]=Texas_SC_Fib_4_IC_Sweep(x_test,Coss_adj,Ron_adj);
                end
                
            end
            Invalid_index = 1;
            while Invalid_index
                [V,I]=min(stick(i,:));
                if x(I) == 1
                    stick(i,I) = 100;
                else
                    Invalid_index = 0;
                end
            end
            if I == i
                %  [V,I]=max(stick(i,:));
                %  x_test = x;
                %  x_test(i) = x(i)-1;
                %  x_test(j) = x(I)+1;
                %  [stick_test,graph1]=Texas_SC_Fib_4_IC_Sweep(x_test,Coss_adj,Ron_adj);
                %  if stick_test < stick(i,i)
                %     x = x_test;
                %      saved_x(end+1,:) = x;
                %     almost_steady_state(I) = 0;
                %      saved_fval(end+1) = stick_test;
                % else
                almost_steady_state(i) = 1;
                % end
                
            else
                x(i) = x(i)+1;
                x(I) = x(I)-1;
                saved_x(end+1,:) = x;
                almost_steady_state(I) = 0;
                almost_steady_state(i) = 0;
                saved_fval(end+1) = stick(i,I);
            end
            
            
            %{
            adjust_FET(i) = 1;
            
            
            
            Coss_adj = [0 0 0 0 0 0 0 0];
            Ron_adj = [0 0 0 0 0 0 0 0];
            
            [stick(i,1),graph1]=Texas_SC_Fib_4_IC_Sweep(x,Coss_adj,Ron_adj);
            
            
            
            [stick(i,2),graph2]=Texas_SC_Fib_4_IC_Sweep(x,Coss_adj,Ron_adj);
            
            
            Coss_adj = [deltaCoss deltaCoss deltaCoss deltaCoss deltaCoss].*adjust_FET;
            Ron_adj = [0 0 0 0 0 ].*adjust_FET;
            
            [stick(i,3),graph3]=Texas_SC_Fib_4_IC_Sweep(x,Coss_adj,Ron_adj);
            
            % Move in the steepest direction:
            stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];
            
            D_Ron = FETs_Ron./deltaRon;
            D_Coss  = FETs_Coss./deltaCoss;
            
            
            Projected_Ploss = stick(i,1)+D_Ron(x(i),:)*stick_diff(1)+D_Coss(x(i),:)*stick_diff(2);
            
            %{
            stick_delta_ploss_over_delta_x = stick_diff./([deltaRon deltaCoss]);
            
            slope = stick_delta_ploss_over_delta_x(1)/stick_delta_ploss_over_delta_x(2);
            
            
            %[V,I] = max(abs(stick_diff));
            
            the_sign = sign(stick_diff);
            
            XY0 = FETs(x(i),:);
            m = slope;
            %XY1 = XY0 + [ron_adjust_values(I)  coss_adjust_values(I)];
            
            lineA = -m;
            lineB = 1;
            lineC = m*XY0(1)-XY0(2);
            
            for j = 1:length(FETs)
                FET_distance(j) = abs(lineA*FETs(j,1)+lineB*FETs(j,2)+lineC)/sqrt(lineA^2+lineB^2);
            end
            
            combine_distance = e_tan_dist.*FET_distance+e_FET_dist.*distance_between_points(x(i),:);
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
            
            
            % Angle calculation
            perp_m = -1/m; % slope
            
            perp_b =  -perp_m*XY0(1)+XY0(2); % y interscept for perpendicular line
           
                perp_y = perp_m.*FETs(:,1)+ perp_b;
                
               GraterThan90 = perp_y<=FETs(:,2); % 1 if angle should have 90 degrees added to it

               
            DegreeFromGreatest = rad2deg(asin(FET_distance./distance_between_points(x(i),:)));
            
            Angle_Error = (DegreeFromGreatest'+GraterThan90.*90)/18;
            
            error_adjusted =combine_distance' +Angle_Error; %+ errorRon*100 + errorCoss*100;
            error_adjusted(x(i)) =  error_adjusted(x(i)) + 100;
            [~,new_xs]=sort(error_adjusted);
            
            
            new_Ploss(1)=sum((XY0-FETs(new_xs(1),:)).*stick_delta_ploss_over_delta_x);
            new_Ploss(2)=sum((XY0-FETs(new_xs(2),:)).*stick_delta_ploss_over_delta_x);
            new_Ploss(3)=sum((XY0-FETs(new_xs(3),:)).*stick_delta_ploss_over_delta_x);
            new_Ploss(4)=sum((XY0-FETs(new_xs(4),:)).*stick_delta_ploss_over_delta_x);
            new_Ploss(5)=sum((XY0-FETs(new_xs(5),:)).*stick_delta_ploss_over_delta_x);
            
            
            [~,new_x_ind] = min(new_Ploss);
            
            new_x=new_xs(new_x_ind);
            x_test = x;
            
            x_test(i) = new_x;
            %}
            [~,new_x] = min(Projected_Ploss);
            
            x_test = x;
            x_test(i) = new_x;
            
            Coss_adj_test = [0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0];
            
            [stick_em,graph1]=Texas_SC_Fib_4_IC_Sweep(x_test,Coss_adj_test,Ron_adj_test);
            
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
            
            
            
            %}
            
            
        end
        big_loop = big_loop + 1;
        
        if sum(almost_steady_state)==length(almost_steady_state)
            toc
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
plot(plot_ron(saved_x(:,1)),plot_Coss(saved_x(:,1)),'DisplayName','M1,M5,M11,M15','LineWidth',3,'LineStyle','--');
hold on
plot(plot_ron(saved_x(:,2)),plot_Coss(saved_x(:,2)),'DisplayName','M2,M4,M12,M14','LineWidth',3,'LineStyle',':');
hold on
plot(plot_ron(saved_x(:,3)),plot_Coss(saved_x(:,3)),'DisplayName','M3,M6,M13,M16','LineWidth',3,'LineStyle','-.');
hold on
plot(plot_ron(saved_x(:,4)),plot_Coss(saved_x(:,4)),'DisplayName','M7&M10','LineWidth',3,'LineStyle','--');
hold on
plot(plot_ron(saved_x(:,5)),plot_Coss(saved_x(:,5)),'DisplayName','M8&M9','LineWidth',3,'LineStyle',':');
hold on
scatter(plot_ron(x(:,1:5)),plot_Coss(x(:,1:5)),50,'r','DisplayName','Ending Point');

% Create ylabel
ylabel('Coss (F)','FontSize',20);

% Create xlabel
xlabel('Ron (\Omega)','FontSize',20);

% Create title
title('SC Fib 4V Converter Optimization','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
legend

