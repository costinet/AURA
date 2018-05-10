% This script tests the circuit parsing code for the boost converter
% "Boost.net"

%{
[filename,path]=uigetfile % opens gui to select file
%}
clear
filename = 'Buck.net';

[A,B,C,D,NLnets,StateNameAB,StateNameCD]=ABCD(filename);



for i=1:1:size(A,3)
    L1 = 16e-6;
    C1 = 40e-6;
    R1 = 10;
    M1_C = 1e-9;
    D1_C = 1e-9;
    M1_R = 0.01;
    D1_R = 0.01;
    R2 = 0.001;
    
    % Switch cases for Boost converter
    
    switch i
        case 1 % M and D off
            %  Inital condidtions
            % V_M1  V_C1  I_L1  V_D1
            
            X = [0 25 9 -50];
                        
        case 2 % M on D off
            
            X = [0.09 50 9 -49.91];
            
        case 3 % M off D on
            
            X = [50 50 9 0];
            
        case 4 % M and D on
            
             X = [25 50 0 -25];
            
    end
    
    
    sys  = ss(eval(A(:,:,i)),eval(B(:,:,i)),eval(C(:,:,i)),eval(D(:,:,i))); % Create state space
    %sys  = ss(eval(A),eval(B),C,D);
%     sys.StateName = StateName;
     sys.OutputName = StateNameCD(:,i);
    
    % Create input and time:
    t = linspace(0,10e-6,100000);
    u = 50*ones(size(t));
    
    Y = lsim(sys,u,t,X);
    
    figure
    h = lsimplot(sys,u,t,X);
    p = getoptions(h);
    for j = 1:1:size(Y,2)
        p.YLim(i) = {[min(Y(100:end,i)) max(Y(100:end,i))]};
    end
    switch i
        case 1 % M and D off
            p.Title.String='M1-OFF and D1-OFF';
        case 2 % M on D off
            p.Title.String='M1-ON and D1-OFF';
        case 3 % M off D on
            p.Title.String='M1-OFF and D1-ON';
        case 4 % M and D on
            p.Title.String='M1-ON and D1-ON';
    end
    setoptions(h,p);
    
end
%{
%     
%     
%     figure
%     p = plot(time,VN001SW,t,Y(:,3));
%     title('Inductor Voltage (L1)','FontSize',14);
%     ylabel('Voltage (V)','FontSize',12);
%     xlabel('Time (s)','FontSize',12);
%     xlim([0 1e-5])
%     legend({'LTSpice','MATLAB'},'FontSize',12);
%     p(1).Color = 'r';
%     p(2).Color = 'b';
%     p(1).LineStyle = '--';
%     p(2).LineStyle = ':';
%     p(1).LineWidth = 2;
%     p(2).LineWidth = 2;
%     
%     figure
%     p = plot(time,IC2,t,Y(:,1));
%     title('FET Capacitor Current (C2)','FontSize',14);
%     ylabel('Current (A)','FontSize',12);
%     xlabel('Time (s)','FontSize',12);
%     xlim([0 1e-5])
%     legend({'LTSpice','MATLAB'},'FontSize',12);
%     p(1).Color = 'r';
%     p(2).Color = 'b';
%     p(1).LineStyle = '--';
%     p(2).LineStyle = ':';
%     p(1).LineWidth = 2;
%     p(2).LineWidth = 2;
%     
%     figure
%     p = plot(time,IC1,t,Y(:,2));
%     title('Output Capacitor Current (C1)','FontSize',14);
%     ylabel('Current (A)','FontSize',12);
%     xlabel('Time (s)','FontSize',12);
%     xlim([0 1e-5])
%     legend({'LTSpice','MATLAB'},'FontSize',12);
%     p(1).Color = 'r';
%     p(2).Color = 'b';
%     p(1).LineStyle = '--';
%     p(2).LineStyle = ':';
%     p(1).LineWidth = 2;
%     p(2).LineWidth = 2;
% 
%     figure
%     p = plot(time,-IC3,t,Y(:,4));
%     title('Diode Capacitor Current (C3)','FontSize',14);
%     ylabel('Current (A)','FontSize',12);
%     xlabel('Time (s)','FontSize',12);
%     xlim([0 1e-5])
%     legend({'LTSpice','MATLAB'},'FontSize',12);
%     p(1).Color = 'r';
%     p(2).Color = 'b';
%     p(1).LineStyle = '--';
%     p(2).LineStyle = ':';
%     p(1).LineWidth = 2;
%     p(2).LineWidth = 2;
%     
%     clear time VN001SW IC1 IC2 IC3
%     
  
  %}