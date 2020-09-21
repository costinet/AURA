% This is the code to control the buck converter used for COMPEL numerical
% comparison

% Jared Baxter
% June 7, 2019

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\


clear
try
    Ls = 6e-6;
    Ron = 0.0015;
    Vg = 12;
    Vout = 5;
    % R_L should be around 0.16 for experimetnal stuff
    I_out = 3;
    Cout = (4*100+680)*10^-6;
    R_out = 1;
    R_L = 0.16;
    Coss = 1620e-12;
    Vfwd = 1;
    
    
    D = 5/12;
    fs_adj = 200e3;
    Ts = 1/fs_adj;
    dead = 0.002;
    M2_D = (1-D)-dead;
    M1_D = D-dead;
    M2_delay = Ts*(D);
    M1_delay = 0;
    
    
    rtol = 1e-6;
    
    %{
% takes 0.02 seconds to reach steady state
tic
Sim_out=sim('Buck_COMPEL_2019',0.003);
toc

tic
plsteadystate('Buck_COMPEL_2019/Steady-State Analysis');
toc
    %}
    
    
    M1_V = 0;
    M2_V = 0;
    L1_I = 0;
    C1_V = 0;
    
    M1_V_data = 0;
    M2_V_data = 0;
    L1_I_data = 0;
    C1_V_data = 0;
    
    
    counter = 0;
    
    type = 2;
    
    tic
    Time_stamp = toc;
    switch type
        %% Strait Lsim
        
        case 1
            
            sim('Buck_COMPEL_2019'); % This is like pressing play in Simulink
            Time_stamp(end+1) = toc;
            C1_V_data(end+1) = C1sim.data(end);
            L1_I_data(end+1) = L1sim.data(end);
            M1_V_data(end+1) = M1sim.data(end);
            M2_V_data(end+1) = M2sim.data(end);
            
            C1_V = C1sim.data(end);
            L1_I = L1sim.data(end);
            M1_V = M1sim.data(end);
            M2_V = M2sim.data(end);
            Total_time = L1sim.time;
            
            
            
            C1_V_data_all = C1sim.data;
            L1_I_data_all = L1sim.data;
            M1_V_data_all = M1sim.data;
            M2_V_data_all = M2sim.data;
            Total_time = L1sim.time;
            
            
            while abs(C1_V_data(end)-C1_V_data(end-1))>1e-6 && abs(L1_I_data(end)-L1_I_data(end-1))>1e-6
                
                sim('Buck_COMPEL_2019'); % This is like pressing play in Simulink
                Time_stamp(end+1) = toc;
                C1_V_data(end+1) = C1sim.data(end);
                L1_I_data(end+1) = L1sim.data(end);
                M1_V_data(end+1) = M1sim.data(end);
                M2_V_data(end+1) = M2sim.data(end);
                
                
                C1_V_data_all(end+1:end+length(C1sim.data)-1) = C1sim.data(2:end);
                L1_I_data_all(end+1:end+length(C1sim.data)-1) = L1sim.data(2:end);
                M1_V_data_all(end+1:end+length(C1sim.data)-1) = M1sim.data(2:end);
                M2_V_data_all(end+1:end+length(C1sim.data)-1) = M2sim.data(2:end);
                Total_time(end+1:end+length(C1sim.data)-1) = (counter*Ts)+L1sim.time(2:end);
                
                counter = counter+1;
                
                
                C1_V = C1sim.data(end);
                L1_I = L1sim.data(end);
                M1_V = M1sim.data(end);
                M2_V = M2sim.data(end);
                
                
                
                
            end
            
            
            C1_error = (abs(C1_V_data(:)-C1_V_data(end))./abs(C1_V_data(end)));
            
            CM1_error =(abs(M1_V_data(:)-M1_V_data(end))./abs(M1_V_data(end)));
            
            CM2_error =(abs(M2_V_data(:)-M2_V_data(end))./abs(M2_V_data(end)));
            
            L_error= (abs(L1_I_data(:)-L1_I_data(end))./abs(L1_I_data(end)));
            
            Errorcode = [C1_error CM1_error CM2_error L_error];
            
            Timecode = [Time_stamp' Time_stamp' Time_stamp' Time_stamp'];
            
            figure(10)
            subplot(3,1,[1:2])
            % Create axes
            
            % Create multiple lines using matrix input to plot
            plot1 = plot(Timecode,Errorcode,'MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerSize',10,...
                'Marker','o',...
                'LineWidth',3);
            set(plot1(1),'DisplayName','C_{out}');
            set(plot1(2),'DisplayName','C_{M1}');
            set(plot1(3),'DisplayName','C_{M2}');
            set(plot1(4),'DisplayName','L');
            
            % Create ylabel
            ylabel('$E_{x_{ss}}$','Interpreter','latex');
            
            set(gca, 'Xticklabel', []);
            set(gca,'FontName','Times New Roman','FontSize',24);
            % Create title
            title({'Error per Time to Reach Steady-State'});
            legend
            
            subplot(3,1,3)
            plot2 = plot(Timecode(:,1:2),zeros(length(Errorcode),2),'MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerSize',10,...
                'Marker','o',...
                'LineWidth',3);
            
            set(plot2(1),'DisplayName','Deadtime1');
            set(plot2(2),'DisplayName','Deadtime2');
            
            
            ylim([-0.49 1.49]);
            
            % Create ylabel
            
            ylabel('$E_{V_{F}}$','Interpreter','latex')
            
            
            % Create xlabel
            xlabel('Real Time [s]');
            
            set(gca,'FontName','Times New Roman','FontSize',24,'YTick',[0 1]);
            legend({},'Orientation','horizontal')
            
            
            figure(2)
            plot(Total_time,C1_V_data_all,Total_time,L1_I_data_all)
            
            
        case 2
            %% Bryodens Method
            
            iterations = 0;
            
            C1_V_data_all = 0;
            L1_I_data_all= 0;
            M1_V_data_all= 0;
            M2_V_data_all= 0;
            Total_time= 0;
            
            while true
                
                while true
                    
                    sim('Buck_COMPEL_2019');
                    Time_stamp(end+1) = toc;
                    C1_V_data(end+1) = C1sim.data(end);
                    L1_I_data(end+1) = L1sim.data(end);
                    M1_V_data(end+1) = M1sim.data(end);
                    M2_V_data(end+1) = M2sim.data(end);
                    
                    C1_V_data_all(end+1:end+length(C1sim.data)-1) = C1sim.data(2:end);
                    L1_I_data_all(end+1:end+length(C1sim.data)-1) = L1sim.data(2:end);
                    M1_V_data_all(end+1:end+length(C1sim.data)-1) = M1sim.data(2:end);
                    M2_V_data_all(end+1:end+length(C1sim.data)-1) = M2sim.data(2:end);
                    Total_time(end+1:end+length(C1sim.data)-1) = (counter*Ts)+L1sim.time(2:end);
                    counter = counter+1;
                    
                    iterations =iterations+1;
                    
                    C1_V = C1sim.data(end);
                    L1_I = L1sim.data(end);
                    M1_V = M1sim.data(end);
                    M2_V = M2sim.data(end);
                    
                    if sum(Diode_condsim.Data(1,:)==Diode_condsim.Data(end,:))==length(Diode_condsim.Data(1,:))
                        break
                    end
                end
                
                fprintf('Circular Topology found\n')
                x = [M1sim.data C1sim.data L1sim.data M2sim.data]';
                max_x = max(abs(x),[],2);
                
                % min_x = [12 4 -3];
                
                min_x = [1 1 1 1];
                
                % Estimate in relative error
                % eta = [abs(12-x(1,end)) abs(4-x(2,end)) abs(-3-x(3,end)) ]';
                eta = [1e-6 1e-6 1e-6 1e-6]';
                
                delta_x(:) = (eta).^0.5.*max(max_x(:),min_x(:));
                
                % From what the default is in PLECS
                delta_x(:) = 1e-4;
                
                x_big = x;
                F_x = x(:,end);
                I = eye(length(F_x));
                J = zeros(length(F_x));
                x = [M1sim.data(1);C1sim.data(1);L1sim.data(1);M2sim.data(1)];
                for i = 1:1:length(x)
                    todeal =  (x+I(:,i)*delta_x(i))';
                    [M1_V, C1_V, L1_I, M2_V]=deal(todeal(1),todeal(2),todeal(3),todeal(4));
                    
                    sim('Buck_COMPEL_2019');
                    
                    % Added to emphasize the delay in find the
                    % jacobians
                    % Time_stamp(end+1) = toc;
                    % C1_V_data(end+1) = C1_V_data(end);
                    % L1_I_data(end+1) = L1_I_data(end);
                    % M1_V_data(end+1) = M1_V_data(end);
                    % M2_V_data(end+1) = M2_V_data(end);
                    
                    
                    
                    new_Fx = [M1sim.data(end) C1sim.data(end) L1sim.data(end) M2sim.data(end)]';
                    % norm(new_x-F_x)/norm(F_x);
                    % J(:,i) = I(:,i)-((new_x-x)./delta_x(i));
                    J(:,i) = I(:,i)-((new_Fx-F_x)./delta_x(i));
                    
                    
                end
                x_k = x;
                
                Alli = 0;
                
                error = true;
                
                while error
                    Alli=Alli+1;
                    x_k = x;
                    
                    x_plus_1 = x_k-(J^-1)*(x_k-F_x);
                    todeal=x_plus_1;
                    [M1_V, C1_V, L1_I, M2_V]=deal(todeal(1),todeal(2),todeal(3),todeal(4));
                    
                    sim('Buck_COMPEL_2019');
                    Time_stamp(end+1) = toc;
                    C1_V_data(end+1) = C1sim.data(end);
                    L1_I_data(end+1) = L1sim.data(end);
                    M1_V_data(end+1) = M1sim.data(end);
                    M2_V_data(end+1) = M2sim.data(end);
                    
                    C1_V_data_all(end+1:end+length(C1sim.data)-1) = C1sim.data(2:end);
                    L1_I_data_all(end+1:end+length(C1sim.data)-1) = L1sim.data(2:end);
                    M1_V_data_all(end+1:end+length(C1sim.data)-1) = M1sim.data(2:end);
                    M2_V_data_all(end+1:end+length(C1sim.data)-1) = M2sim.data(2:end);
                    Total_time(end+1:end+length(C1sim.data)-1) = (counter*Ts)+L1sim.time(2:end);
                    
                    counter = counter+1;
                    
                    
                    
                    F_x_tot = [M1sim.data C1sim.data L1sim.data M2sim.data]';
                    F_x = F_x_tot(:,end);
                    % Broyden's Update
                    
                    J = J+(((x_plus_1-F_x)*(x_plus_1-x_k)')/(norm(x_plus_1-x_k))^2);
                    
                    
                    
                    
                    
                    if sum(Diode_condsim.Data(1,:)==Diode_condsim.Data(end,:))~=length(Diode_condsim.Data(1,:))
                        C1_V = C1_V_data(end);
                        L1_I = L1_I_data(end);
                        M1_V = M1_V_data(end);
                        M2_V = M2_V_data(end);
                        fprintf('Non-Circular Topology found\n')
                        disp(counter);
                        break
                    end
                    
                    
                    error  = any(abs(x_plus_1-F_x)./max(F_x_tot,[],2)>rtol & abs((x_plus_1-x_k)./x_k)>rtol);
                    
                    x = x_plus_1;
                    
                end
                
                if sum(Diode_condsim.Data(1,:)==Diode_condsim.Data(end,:))==length(Diode_condsim.Data(1,:))
                    break
                end
                
            end
            
            trace = [M1_V_data' M2_V_data' C1_V_data' L1_I_data'];
            
            
            
            C1_error = (abs(C1_V_data(:)-C1_V_data(end))./abs(C1_V_data(end)));
            
            CM1_error =(abs(M1_V_data(:)-M1_V_data(end))./abs(M1_V_data(end)));
            
            CM2_error =(abs(M2_V_data(:)-M2_V_data(end))./abs(M2_V_data(end)));
            
            L_error= (abs(L1_I_data(:)-L1_I_data(end))./abs(L1_I_data(end)));
            
            Errorcode = [C1_error CM1_error CM2_error L_error];
            
            
            Timecode = [Time_stamp' Time_stamp' Time_stamp' Time_stamp'];
            
            
            figure(11)
            
            subplot(3,1,[1:2])
            % Create axes
            
            % Create multiple lines using matrix input to plot
            plot1 = plot(Timecode,Errorcode,'MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerSize',10,...
                'Marker','o',...
                'LineWidth',3);
            set(plot1(1),'DisplayName','C_{out}');
            set(plot1(2),'DisplayName','C_{M1}');
            set(plot1(3),'DisplayName','C_{M2}');
            set(plot1(4),'DisplayName','L');
            
            % Create ylabel
            ylabel('$E_{x_{ss}}$','Interpreter','latex');
            
            set(gca, 'Xticklabel', []);
            set(gca,'FontName','Times New Roman','FontSize',24);
            % Create title
            title({'Error per Time to Reach Steady-State'});
            legend
            
            subplot(3,1,3)
            plot2 = plot(Timecode(:,1:2),zeros(length(Errorcode),2),'MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerSize',10,...
                'Marker','o',...
                'LineWidth',3);
            
            set(plot2(1),'DisplayName','Deadtime1');
            set(plot2(2),'DisplayName','Deadtime2');
            
            
            ylim([-0.49 1.49]);
            
            % Create ylabel
            
            ylabel('$E_{V_{F}}$','Interpreter','latex')
            
            
            % Create xlabel
            xlabel('Real Time [s]');
            
            set(gca,'FontName','Times New Roman','FontSize',24,'YTick',[0 1]);
            legend({},'Orientation','horizontal')
            
            
            
            
            
            figure(4)
            plot(Total_time,C1_V_data_all,Total_time,L1_I_data_all,Total_time,M2_V_data_all)
            
            
            
    end
    
    %% Mesh grid of the Buck Converter
    
    % Span of the Variables that are to be investigated
    %{
    
    M1_span = linspace(11.9,12.1,3);
    M1_span = 11.9958;
    C1_span = linspace(3.5,5,250);
    L1_span = linspace(-5,0,250);
    M2_span = Vg-M1_span;
    
    normall = zeros(length(M1_span),length(C1_span),length(L1_span));
    for M1 = 1:1:length(M1_span)
        for C1 = 1:1:length(C1_span)
            for L1 = 1:1:length(L1_span)
                
                
                M1_V = M1_span(M1);
                M2_V = M2_span(M1);
                C1_V = C1_span(C1);
                L1_V = L1_span(L1);
                sim('Buck_COMPEL_2019');
                
                error_M1 = M1_V-M1sim.data(end);
                error_C1 = C1_V-C1sim.data(end);
                error_L1 = L1_V-L1sim.data(end);
                error_all = [error_M1 error_C1 error_L1];
                erros=norm(error_all);
                
                normall(M1,C1,L1)=erros;
                
            end
        end
    end
    
    new_var(:,:)=(normall(1,:,:));
    figure
    s = surf(L1_span,C1_span,new_var);
    
    
    
    %     light               % add a light
    %     lighting gouraud    % preferred lighting for a curved surface
    %     view(40,30)         % set viewpoint
    %     camzoom(1.5)        % zoom into scen
    
    
    %     for i = 1:1:5
    %         new_var(:,:)=(normall(i,:,:));
    %         s.ZData = new_var;    % replace surface z values
    %         pause(1)
    %     end
    
    %}
    Alex = 4564821;
    
    
catch ME
    rethrow(ME)
end

return


%{
figure(1)

ts_hist1 = ts_hist(:,1);
ts_hist2 = ts_hist(:,2);
ts_hist3 = ts_hist(:,3);
ts_hist4 = ts_hist(:,4);
ts_hist5 = ts_hist(:,5);
ts_hist6 = ts_hist(:,6);


plot(abs(ts_hist1(:)-ts_hist1(end))./abs(ts_hist1(end)))
hold on
plot(abs(ts_hist2(:)-ts_hist2(end))./abs(ts_hist2(end)))
hold on
plot(abs(ts_hist3(:)-ts_hist3(end))./abs(ts_hist3(end)))
hold on
plot(abs(ts_hist4(:)-ts_hist4(end))./abs(ts_hist4(end)))
hold on
plot(abs(ts_hist5(:)-ts_hist5(end))./abs(ts_hist5(end)))
hold on
plot(abs(ts_hist6(:)-ts_hist6(end))./abs(ts_hist6(end)))
%}

load('COMPEL_2020_Digest_AURA_Data.mat')
play = Xs_hist;


%play  = play(:,:,2:end);
play(:,:,1)  = play(:,:,2);

error1 = abs(squeeze(play(3,3,:)+1));
error2 = abs(squeeze(play(3,6,:)+1));


Errorcode = [error1 error2];

Time_stamp = [0 Time_stamp];

Timecode = [Time_stamp' Time_stamp' Time_stamp' Time_stamp'];

figure(13)




subplot(3,1,[1:2])
% Create axes

% Create multiple lines using matrix input to plot

plot1 = plot(Timecode(:,1:2),Errorcode,'MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerSize',10,...
                'Marker','o',...
                'LineWidth',3);

set(plot1(1),'DisplayName','Deadtime1');
set(plot1(2),'DisplayName','Deadtime2');

xlim([0 1]);
ylabel('$E_{V_{F}}$','Interpreter','latex')




% Create ylabel

set(gca, 'Xticklabel', []);
set(gca,'FontName','Times New Roman','FontSize',24);
% Create title
title({'Error per Time to Reach Steady-State'});
legend

subplot(3,1,3)
plot2 = plot(Timecode,zeros(length(Errorcode),4),'MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerSize',10,...
                'Marker','o',...
                'LineWidth',3);
set(plot2(1),'DisplayName','C_{out}');
set(plot2(2),'DisplayName','C_{M1}');
set(plot2(3),'DisplayName','C_{M2}');
set(plot2(4),'DisplayName','L');


xlim([0 1]);
ylim([-0.49 1.49]);

ylabel('$E_{x_{ss}}$','Interpreter','latex');

% Create xlabel
xlabel('Real Time [s]');

set(gca,'FontName','Times New Roman','FontSize',24,'YTick',[0 1]);
legend({},'Orientation','horizontal')



