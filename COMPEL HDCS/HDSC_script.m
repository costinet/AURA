clear
tic
Vg = 24;
V = 5;
Io = 10;

L = 2*72e-9;
RL = (0.007-2.65e-3 + 2*.35e-3);
Da = 1;
M = 5/6;

fs = 500e3;
Ts = 1/fs;

dt = 1/fs/Ts;

Cout = 100e-6;
C1 = 13.38e-6;
% C1 = 1e-6;
C2 = 12.49e-6;
C3 = 10.27e-6;

ESRx = 3e-3;
ESR1 = 1.1e-3+ESRx;
ESR2 = 2.3e-3+ESRx;
ESR3 = 2.3e-3+ESRx;
ron = 2.5e-3;

Lp = 0e-9;
Lp1 = Lp;
Lp2 = 2*Lp;
Lp3 = 3*Lp;


kc = 1;
Coss = kc*2.5e-9;
Coss6 = kc*2.26e-9;
Coss12 = kc*1.99e-9;
Qg = 18e-9;

Phase1a = [1 0 1 1 0 0 0 1];
Phase1b = [1 0 0 1 0 0 0 1];
Phase2a = 1 - Phase1a;
Phase2b = Phase2a; Phase2b(6) = 0;
Phase3 = [0 0 0 1 1 0 1 1];
swseq = [Phase1a; Phase1b; Phase3;  Phase2a; Phase2b;  Phase3];

history = [];
stateElem = [L Cout C1 C2 C3 Coss12 Coss12 Coss6 Coss6 Coss6 Coss6 Coss6 Coss6 ];
top = SMPStopology;
% top.loadPLECsModel('COMPEL_2019_HDSC/HDSC',swseq)

%{
u = [Vg Io]';
conv.u = u;
sim.converter.ts = [Da*Ts/2*M, (1-Da)*Ts/2*M, Ts/2*(1-M), Da*Ts/2*M, (1-Da)*Ts/2*M, Ts/2*(1-M)];
%}

rounded = 20;

FET4Duty = round((Da*Ts/2*M)/Ts,rounded);
FET4Delay = 0;

FET3Duty = round((Da*Ts/2*M+(1-Da)*Ts/2*M)/Ts,rounded);
FET3Delay = round(Da*Ts/2*M+(1-Da)*Ts/2*M+Ts/2*(1-M),rounded);

FET2Duty = round((Da*Ts/2*M+(1-Da)*Ts/2*M)/Ts,rounded);
FET2Delay = 0;

FET1Duty = round((Da*Ts/2*M)/Ts,rounded);
FET1Delay = round(Da*Ts/2*M+(1-Da)*Ts/2*M+Ts/2*(1-M),rounded);

FET58Duty = round((Da*Ts/2*M+(1-Da)*Ts/2*M)/Ts,rounded);
FET58Delay = round(Da*Ts/2*M+(1-Da)*Ts/2*M+Ts/2*(1-M),rounded);

FET67Duty = round((Ts/2*(1-M)+Da*Ts/2*M+(1-Da)*Ts/2*M+Ts/2*(1-M))/Ts,rounded);
FET67Delay = round(Da*Ts/2*M+(1-Da)*Ts/2*M,rounded);


% Sim_out=sim('COMPEL_2019_HDSC');


type = 1;
switch type
    %% Strait Lsim
    
    case 1
        
        
        C1_V_data = 0;
        L1_I_data = 0;
        Coss1_V_data = 0;
        Coss2_V_data = 0;
        Coss3_V_data = 0;
        Coss4_V_data = 0;
        Coss5_V_data = 0;
        Coss6_V_data = 0;
        Coss7_V_data = 0;
        Coss8_V_data = 0;
        VC1_V_data = 0;
        VC2_V_data = 0;
        VC3_V_data = 0;
        
        
        C1_V = 0;
        L1_I = 0;
        Coss1_V = 0;
        Coss2_V = 0;
        Coss3_V = 0;
        Coss4_V = 0;
        Coss5_V = 0;
        Coss6_V = 0;
        Coss7_V = 0;
        Coss8_V = 0;
        VC1_V = 0;
        VC2_V = 0;
        VC3_V = 0;
        
        
        
        
        
        
        sim('COMPEL_2019_HDSC'); % This is like pressing play in Simulink
        toc
        C1_V_data(end+1) = C1sim.data(end);
        L1_I_data(end+1) = L1sim.data(end);
        Coss1_V_data(end+1) = Coss1sim.data(end);
        Coss2_V_data(end+1) = Coss2sim.data(end);
        Coss3_V_data(end+1) = Coss3sim.data(end);
        Coss4_V_data(end+1) = Coss4sim.data(end);
        Coss5_V_data(end+1) = Coss5sim.data(end);
        Coss6_V_data(end+1) = Coss6sim.data(end);
        Coss7_V_data(end+1) = Coss7sim.data(end);
        Coss8_V_data(end+1) = Coss8sim.data(end);
        VC1_V_data(end+1) = VC1sim.data(end);
        VC2_V_data(end+1) = VC2sim.data(end);
        VC3_V_data(end+1) = VC3sim.data(end);
        
        
        C1_V = C1sim.data(end);
        L1_I = L1sim.data(end);
        Coss1_V = Coss1sim.data(end);
        Coss2_V = Coss2sim.data(end);
        Coss3_V = Coss3sim.data(end);
        Coss4_V = Coss4sim.data(end);
        Coss5_V = Coss5sim.data(end);
        Coss6_V = Coss6sim.data(end);
        Coss7_V = Coss7sim.data(end);
        Coss8_V = Coss8sim.data(end);
        VC1_V = VC1sim.data(end);
        VC2_V = VC2sim.data(end);
        VC3_V = VC3sim.data(end);
        
        
        while abs(C1_V_data(end)-C1_V_data(end-1))>1e-6 && abs(L1_I_data(end)-L1_I_data(end-1))>1e-6
            
            sim('COMPEL_2019_HDSC'); % This is like pressing play in Simulink
            
            C1_V_data(end+1) = C1sim.data(end);
            L1_I_data(end+1) = L1sim.data(end);
            Coss1_V_data(end+1) = Coss1sim.data(end);
            Coss2_V_data(end+1) = Coss2sim.data(end);
            Coss3_V_data(end+1) = Coss3sim.data(end);
            Coss4_V_data(end+1) = Coss4sim.data(end);
            Coss5_V_data(end+1) = Coss5sim.data(end);
            Coss6_V_data(end+1) = Coss6sim.data(end);
            Coss7_V_data(end+1) = Coss7sim.data(end);
            Coss8_V_data(end+1) = Coss8sim.data(end);
            VC1_V_data(end+1) = VC1sim.data(end);
            VC2_V_data(end+1) = VC2sim.data(end);
            VC3_V_data(end+1) = VC3sim.data(end);
            
            
            C1_V = C1sim.data(end);
            L1_I = L1sim.data(end);
            Coss1_V = Coss1sim.data(end);
            Coss2_V = Coss2sim.data(end);
            Coss3_V = Coss3sim.data(end);
            Coss4_V = Coss4sim.data(end);
            Coss5_V = Coss5sim.data(end);
            Coss6_V = Coss6sim.data(end);
            Coss7_V = Coss7sim.data(end);
            Coss8_V = Coss8sim.data(end);
            VC1_V = VC1sim.data(end);
            VC2_V = VC2sim.data(end);
            VC3_V = VC3sim.data(end);
            
        end
        toc
        figure(1)
        
        plot(abs(C1_V_data(:)-C1_V_data(end))./abs(C1_V_data(end)))
        hold on
        plot(abs(L1_I_data(:)-L1_I_data(end))./abs(L1_I_data(end)))
        hold on
        plot(abs(Coss1_V_data(:)-Coss1_V_data(end))./abs(Coss1_V_data(end)))
        hold on
        plot(abs(Coss2_V_data(:)-Coss2_V_data(end))./abs(Coss2_V_data(end)))
        hold on
        plot(abs(Coss3_V_data(:)-Coss3_V_data(end))./abs(Coss3_V_data(end)))
        hold on
        plot(abs(Coss4_V_data(:)-Coss4_V_data(end))./abs(Coss4_V_data(end)))
        hold on
        plot(abs(Coss5_V_data(:)-Coss5_V_data(end))./abs(Coss5_V_data(end)))
        hold on
        plot(abs(Coss6_V_data(:)-Coss6_V_data(end))./abs(Coss6_V_data(end)))
        hold on
        plot(abs(Coss7_V_data(:)-Coss7_V_data(end))./abs(Coss7_V_data(end)))
        hold on
        plot(abs(Coss8_V_data(:)-Coss8_V_data(end))./abs(Coss8_V_data(end)))
        hold on
        plot(abs(VC1_V_data(:)-VC1_V_data(end))./abs(VC1_V_data(end)))
        hold on
        plot(abs(VC2_V_data(:)-VC2_V_data(end))./abs(VC2_V_data(end)))
        hold on
        plot(abs(VC3_V_data(:)-VC3_V_data(end))./abs(VC3_V_data(end)))
        
    case 2
        %% Bryodens Method
        
        C1_V_data = 0;
        L1_I_data = 0;
        Coss1_V_data = 0;
        Coss2_V_data = 0;
        Coss3_V_data = 0;
        Coss4_V_data = 0;
        Coss5_V_data = 0;
        Coss6_V_data = 0;
        Coss7_V_data = 0;
        Coss8_V_data = 0;
        VC1_V_data = 0;
        VC2_V_data = 0;
        VC3_V_data = 0;
        
        
        C1_V = 0;
        L1_I = 0;
        Coss1_V = 0;
        Coss2_V = 0;
        Coss3_V = 0;
        Coss4_V = 0;
        Coss5_V = 0;
        Coss6_V = 0;
        Coss7_V = 0;
        Coss8_V = 0;
        VC1_V = 0;
        VC2_V = 0;
        VC3_V = 0;
        
        
        sim('COMPEL_2019_HDSC'); % This is like pressing play in Simulink
        
        C1_V_data(end+1) = C1sim.data(end);
        L1_I_data(end+1) = L1sim.data(end);
        Coss1_V_data(end+1) = Coss1sim.data(end);
        Coss2_V_data(end+1) = Coss2sim.data(end);
        Coss3_V_data(end+1) = Coss3sim.data(end);
        Coss4_V_data(end+1) = Coss4sim.data(end);
        Coss5_V_data(end+1) = Coss5sim.data(end);
        Coss6_V_data(end+1) = Coss6sim.data(end);
        Coss7_V_data(end+1) = Coss7sim.data(end);
        Coss8_V_data(end+1) = Coss8sim.data(end);
        VC1_V_data(end+1) = VC1sim.data(end);
        VC2_V_data(end+1) = VC2sim.data(end);
        VC3_V_data(end+1) = VC3sim.data(end);
        
        
        x = [C1sim.data L1sim.data Coss1sim.data Coss2sim.data Coss3sim.data Coss4sim.data Coss5sim.data Coss6sim.data Coss7sim.data Coss8sim.data VC1sim.data VC2sim.data VC3sim.data]';
        max_x = max(abs(x),[],2);
        
        
        
        min_x = [5 5 5 5 5 5 5 5 5 5 5 5 5];
        
        % Estimate in relative error
        
        eta = [1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3]';
        
        delta_x(:) = (eta).^0.5.*max(max_x(:),min_x(:));
        
        
        x_big = x;
        F_x = x(:,end);
        I = eye(length(F_x));
        J = zeros(length(F_x));
        x = [0;0;0;0;0;0;0;0;0;0;0;0;0];
        for i = 1:1:length(x)
            todeal =  (x+I(:,i)*delta_x(i))';
            [C1_V L1_I Coss1_V Coss2_V Coss3_V Coss4_V Coss5_V Coss6_V Coss7_V Coss8_V VC1_V VC2_V VC3_V]=deal(todeal(1),todeal(2),todeal(3),todeal(4),todeal(5),todeal(6),todeal(7),todeal(8),todeal(9),todeal(10),todeal(11),todeal(12),todeal(13));
            sim('COMPEL_2019_HDSC');
            new_x = [ C1sim.data(end)  L1sim.data(end)  Coss1sim.data(end)  Coss2sim.data(end)  Coss3sim.data(end)  Coss4sim.data(end)  Coss5sim.data(end)  Coss6sim.data(end)  Coss7sim.data(end)  Coss8sim.data(end)  VC1sim.data(end) VC2sim.data(end) VC3sim.data(end)]';
            norm(new_x-F_x)/norm(F_x);
            J(:,i) = I(:,i)-((new_x-x)./delta_x(i));
            
        end
        x_k = zeros(length(new_x),1);
        
        Alli = 0;
        while abs(C1_V_data(end)-C1_V_data(end-1))>1e-6 && abs(L1_I_data(end)-L1_I_data(end-1))>1e-6
            Alli=Alli+1;
            x_plus_1 = x_k-(J^-1)*(x_k-F_x);
            todeal=x_plus_1;
            [C1_V L1_I Coss1_V Coss2_V Coss3_V Coss4_V Coss5_V Coss6_V Coss7_V Coss8_V VC1_V VC2_V VC3_V]=deal(todeal(1),todeal(2),todeal(3),todeal(4),todeal(5),todeal(6),todeal(7),todeal(8),todeal(9),todeal(10),todeal(11),todeal(12),todeal(13));
            sim('COMPEL_2019_HDSC');
            
            C1_V_data(end+1) = C1sim.data(end);
            L1_I_data(end+1) = L1sim.data(end);
            Coss1_V_data(end+1) = Coss1sim.data(end);
            Coss2_V_data(end+1) = Coss2sim.data(end);
            Coss3_V_data(end+1) = Coss3sim.data(end);
            Coss4_V_data(end+1) = Coss4sim.data(end);
            Coss5_V_data(end+1) = Coss5sim.data(end);
            Coss6_V_data(end+1) = Coss6sim.data(end);
            Coss7_V_data(end+1) = Coss7sim.data(end);
            Coss8_V_data(end+1) = Coss8sim.data(end);
            VC1_V_data(end+1) = VC1sim.data(end);
            VC2_V_data(end+1) = VC2sim.data(end);
            VC3_V_data(end+1) = VC3sim.data(end);
            
            
            F_x = [C1sim.data L1sim.data Coss1sim.data Coss2sim.data Coss3sim.data Coss4sim.data Coss5sim.data Coss6sim.data Coss7sim.data Coss8sim.data VC1sim.data VC2sim.data VC3sim.data]';
            F_x = F_x(:,end);
            % Broyden's Update
            
            J = J+(((x_plus_1-F_x)*(x_plus_1-x_k)')/(norm(x_plus_1-x_k))^2);
            
            x_k = x_plus_1;
            
            
        end
        
        
        % trace = [M1_V_data' M2_V_data' C1_V_data' L1_I_data'];
        
        
        figure(2)
        
        plot(abs(C1_V_data(:)-C1_V_data(end))./abs(C1_V_data(end)))
        hold on
        plot(abs(L1_I_data(:)-L1_I_data(end))./abs(L1_I_data(end)))
        hold on
        plot(abs(Coss1_V_data(:)-Coss1_V_data(end))./abs(Coss1_V_data(end)))
        hold on
        plot(abs(Coss2_V_data(:)-Coss2_V_data(end))./abs(Coss2_V_data(end)))
        hold on
        plot(abs(Coss3_V_data(:)-Coss3_V_data(end))./abs(Coss3_V_data(end)))
        hold on
        plot(abs(Coss4_V_data(:)-Coss4_V_data(end))./abs(Coss4_V_data(end)))
        hold on
        plot(abs(Coss5_V_data(:)-Coss5_V_data(end))./abs(Coss5_V_data(end)))
        hold on
        plot(abs(Coss6_V_data(:)-Coss6_V_data(end))./abs(Coss6_V_data(end)))
        hold on
        plot(abs(Coss7_V_data(:)-Coss7_V_data(end))./abs(Coss7_V_data(end)))
        hold on
        plot(abs(Coss8_V_data(:)-Coss8_V_data(end))./abs(Coss8_V_data(end)))
        hold on
        plot(abs(VC1_V_data(:)-VC1_V_data(end))./abs(VC1_V_data(end)))
        hold on
        plot(abs(VC2_V_data(:)-VC2_V_data(end))./abs(VC2_V_data(end)))
        hold on
        plot(abs(VC3_V_data(:)-VC3_V_data(end))./abs(VC3_V_data(end)))
        
        
end


