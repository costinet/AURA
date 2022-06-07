
clear
tic
try
    FETs = [
        1.45e-3   1530e-12
        16e-3      150e-12
        4e-3       710e-12
        3.6e-3     408e-12
        2.4e-3    1120e-12
        1.5e-3    1620e-12
        2.6e-3     980e-12
        2.2e-3    1020e-12
        3.6e-3     534e-12
        3.2e-3     820e-12
        2.2e-3    1100e-12
        9e-3       267e-12
        13.5e-3     195e-12 
        6e-3        304e-12
        3.2e-3      562e-12
        ];
    
    FETs(:,1) = FETs(:,1)./1.45e-3;
    FETs(:,2) = FETs(:,2)./150e-12;
    
    e_tan_dist = 1;
    e_FET_dist = 1;
    
    for i = 1:length(FETs)
        for j = 1:length(FETs)
            
            distance_between_points(i,j) = sqrt((FETs(j,1)-FETs(i,1))^2+(FETs(j,2)-FETs(i,2))^2);
            
        end
    end
    
    saved_fval = ones(6,1).*100;
    big_loop = 0;
    max_iteration = 3;
    almost_steady_state = [0 0 0 0 0 0];
    
    
    x = [2 3 4 2 3 4    2.0000    2.4372    0.2737    0.7965];
    
    saved_x = x;
    FETs_number = 6;
    stick = zeros(6,3);
    coss_adjust_values = [0 -0.1e-9/150e-12];
    ron_adjust_values = [-1e-3/1.45e-3 0];
    deltaRon = -1e-3;
    deltaCoss = -0.1e-9;
    
    
    while(big_loop<max_iteration)
        for i = 1:FETs_number
            
            adjust_FET = [0 0 0 0 0 0];
            adjust_FET(i) = 1;
            
            Coss_adj = [0 0 0 0 0 0];
            Ron_adj = [0 0 0 0 0 0];
            
            [stick(i,1),graph1]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);
            
            Coss_adj = [0 0 0 0 0 0];
            Ron_adj = [-1e-3 -1e-3 -1e-3 -1e-3 -1e-3 -1e-3];
            
            [stick(i,2),graph2]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);
            
            
            Coss_adj = [-0.1e-9 -0.1e-9 -0.1e-9 -0.1e-9 -0.1e-9 -0.1e-9];
            Ron_adj = [0 0 0 0 0 0];
            
            [stick(i,3),graph3]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);
            
            % Move in the steepest direction:
            stick_diff = [stick(i,1)-stick(i,2) stick(i,1)-stick(i,3)];
            stick_delta_ploss_over_delta_x = stick_diff./(coss_adjust_values+ron_adjust_values);
            
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
            
            perp_b =  -perp_m*XY0(1)+XY0(2); % y interscept
           
                perp_y = perp_m.*FETs(:,1)+ perp_b;
                
               GraterThan90 = perp_y<=FETs(:,2); % 1 if angle shoul have 90 degrees added to it 

               
            DegreeFromGreatest = rad2deg(asin(FET_distance./distance_between_points(x(i),:)));
            
            Angle_Error = (DegreeFromGreatest'+GraterThan90.*90)/18;
            
            error_adjusted =combine_distance' +Angle_Error; %+ errorRon*100 + errorCoss*100;
            error_adjusted(x(i)) =  error_adjusted(x(i)) + 100;
            [~,new_x]=min(error_adjusted);
            
            x_test = x;
            
            x_test(i) = new_x;
            
            Coss_adj_test = [0 0 0 0 0 0];
            Ron_adj_test = [0 0 0 0 0 0];
            
            [stick_em,graph1]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x_test,Coss_adj_test,Ron_adj_test);
            
            % If there is not an improvement to the fval of the converter.
            if stick_em > stick(i,1)
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
        big_loop = big_loop + 1;
        
        if sum(almost_steady_state)==length(almost_steady_state)
            toc
            J = 4564654;
        end
    end
    
    
    
    
catch ME
    J = 5465456465;
    
end





