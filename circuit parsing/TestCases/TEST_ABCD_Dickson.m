function [] = TEST_ABCD_Dickson(parse)

A = parse.Asym;
B = parse.Bsym;
C = parse.Csym;
D = parse.Dsym;

StateNameAB = parse.StateNames;
StateNameCD = parse.OutputNames;

SortedTree = parse.SortedTree;
SortedCoTree = parse.SortedCoTree;

L1 = 16e-6;
C1 = 40e-6;
R1 = 10;
M1_C = 1e-9;
M2_C = 1e-9;
M3_C = 1e-9;
M4_C = 1e-9;
M5_C = 1e-9;
M6_C = 1e-9;
M7_C = 1e-9;
M8_C = 1e-9;
M9_C = 1e-9;
M10_C = 1e-9;
M11_C = 1e-9;
M12_C = 1e-9;
M1_R = 0.01;
M2_R = 0.01;
M3_R = 0.01;
M4_R = 0.01;
M5_R = 0.01;
M6_R = 0.01;
M7_R = 0.01;
M8_R = 0.01;
M9_R = 0.01;
M10_R = 0.01;
M11_R = 0.01;
M12_R = 0.01;
C2 = 40e-6;
C3 = 40e-6;
C4 = 40e-6;
C5 = 40e-6;
C6 = 40e-6;
C7 = 40e-6;
C8 = 40e-6;


if isempty(A)
    for k = 1:1:size(parse.HtempAB,3)
        HtempAB(:,:,k) = eval(parse.HtempAB(:,:,k));
        HtempCD(:,:,k) = eval(parse.HtempCD(:,:,k));
        dependsAB(:,:,k) = eval(parse.dependsAB(:,:,k));
        savedCD(:,:,k) = eval(parse.savedCD(:,:,k));
        for j = 1:1:size(parse.DependentNames(:,k),1)
            DependentNames(j,k) = eval(parse.DependentNames{j,k});
        end
        for j = 1:1:size(parse.OutputNames(:,k),1)
            OutputNames(j,k) = eval(parse.OutputNames{j,k});
        end
    end
    for k = 1:1:size(parse.HtempAB,3)
        [A,B,C,D] = parse.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
        [C,D] = parse.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),SortedTree(:,:,k),SortedCoTree(:,:,k));
        parse.Anum(:,:,k)=A;
        parse.Bnum(:,:,k)=B;
        parse.Cnum(:,:,k)=C;
        parse.Dnum(:,:,k)=D;
    end
    
else
    for k = 1:1:size(A,3)
        parse.Anum(:,:,k) = eval(A(:,:,k));
        parse.Bnum(:,:,k) = eval(B(:,:,k));
        parse.Cnum(:,:,k) = eval(C(:,:,k));
        parse.Dnum(:,:,k) = eval(D(:,:,k));
    end
end


for i=1:1:size(parse.Anum,3)
    
    
    sys  = ss(parse.Anum(:,:,i),parse.Bnum(:,:,i),parse.Cnum(:,:,i),parse.Dnum(:,:,i)); % Create state space
    
    
    % Switch cases for Boost converter
    
    switch i
        case 1 % M and D off
            %  Inital condidtions
            % V_M1  V_C1  I_L1  V_D1
            
            X = [50 50 9 0];
            
        case 2 % M on D off
            
            X = [0.09 50 9 -49.91];
            
        case 3 % M off D on
            
            X = [50 50 9 0];
            
        case 4 % M and D on
            
            X = [25 50 0 -25];
            
    end
    
    
    
    %sys  = ss(eval(A),eval(B),C,D);
    %     sys.StateName = StateName;
    sys.OutputName = parse.OutputNamesCD(:,i);
    
    % Create input and time:
    t = linspace(0,10e-6,100000);
    u = 25*ones(size(t));
    
    [Y] = lsim(sys,u,t,X);
    
    
    
    figure
    h = lsimplot(sys,u,t,X);
    p = getoptions(h);
    for j = 1:1:size(Y,2)
        if Y(1,j) ~= Y(end,j)
            p.YLim(j) = {[min(Y(100:end,j)) max(Y(100:end,j))]};
        end
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




