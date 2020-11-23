clear all; 
clc;

addpath('E:\Dropbox\UTK\Research\Hybrid Dickson\old AURA')
addpath('C:\Users\dcostine\Dropbox\UTK\Research\Hybrid Dickson\old AURA')

plotIterWF = 1;

Vg = 24;
V = 5;
Io = 4.9;
Iorange = logspace(log10(.5), log10(18), 25);

eta = [];
Psw = [];
Ploss = [];
 
% R = V/Io;
it = 1;
% for ron = [1e-3]
% for Rl = [0.49e-3 0.51e-3]
% for C1 = [10e-6*3e6/fs*.9 10e-6*3e6/fs*1.1]
fs = .5e6;
Ts = 1/fs;

MaeveData = [2.579472	0.839671875;
20.76567	0.926377141;
24.77316	0.923269231;
30.681054	0.915743016;
34.88422	0.910491836;
48.143252	0.89233489;
59.17632	0.873425434;
66.39405	0.857804264;
74.51991	0.839188176];

 
L = 143e-9;
Cout = 100e-6;
C1 =12e-6;
C2 = 14e-6;
C3 = 16e-6;
% ESR1 = 35e-3;
% ESR2 = 35e-3;
% ESR3 = 35e-3;
ESR1 = 5e-3;
ESR2 = 5e-3;
ESR3 = 5e-3;
ron = 1.2e-3;
Rl = 6e-3;
Coss = 2.5e-9;
Qg = 18e-9;

origParams = [Rl Qg Coss ron L C1 fs ESR1 ESR2 ESR3];

partials = [Rl*.99 Rl*1.01;
            Qg*.99 Qg*1.01;
            Coss*.99 Coss*1.01;
            ron*.99 ron*1.01;
            L*.99 L*1.01;
            C1*.99 C1*1.01;
            fs*.99 fs*1.01
            ESR1*.99 ESR1*1.01
            ESR2*.99 ESR2*1.01
            ESR3*.99 ESR3*1.01]';
        
simulator = SMPSim();

        
for it = 1:length(partials(:))

    t = num2cell(origParams);
    [Rl, Qg, Coss, ron, L, C1, fs, ESR1, ESR2, ESR3] = deal(t{:});

    if it <= 2
        Rl = partials(it);
    elseif it <= 4
        Qg = partials(it);
    elseif it <= 6
        Coss = partials(it);
    elseif it<=8
        ron = partials(it);
    elseif it<=10
        L = partials(it);
    elseif it<=12
        C1 = partials(it);
        C2 = C1;
        C3 = C1;
    elseif it<=14
        fs = partials(it);
        Ts = 1/fs;
    elseif it<=16
        ESR1 = partials(it);
    elseif it<=18
        ESR2 = partials(it);
    elseif it<=20
        ESR3 = partials(it);
    end
    
    it2=1; 

    % ts = [Ts*.5-dt dt Ts*.5-dt dt];


    %% State Space Construction
    % x = [V1 V2 V3 Vout IL]
    u = [Vg Io]';

    ic3 = [1, -1, -1, 0, 3*ron + ESR1 + ESR2]/(5*ron + ESR1 + ESR2 + ESR3);
    ic1 = -ic3 + [0 0 0 0 1]; %iL-ic3
    ic2 = -ic1;

    A1a = cat(1, ...
    ic1, ...
    ic2, ...
    ic3, ...
    [0, 0, 0, 0, 1], ...
    -ic3*(2*ron + ESR3) - [0 0 1 1 Rl] );

    B1a = cat(1, ...
    [-1/(5*ron + ESR1 + ESR2 + ESR3), 0], ...
    [1/(5*ron + ESR1 + ESR2 + ESR3), 0], ...
    [1/(5*ron + ESR1 + ESR2 + ESR3), 0], ...
    [0, -1], ...
    [-1/(5*ron + ESR1 + ESR2 + ESR3)*(2*ron + ESR3) + 1, 0]);


    %             X0 = [6.089, 12.183, 17.432, 4.894, 1.6655];
    %             vl = A1a(end,:)*X0'+B1a(end,:)*u
    % 
    %             stateElem = [C1 C2 C3 Cout L];
    %             K = diag(stateElem);
    %             A1a = K^-1*A1a;
    %             B1a = K^-1*B1a;
    % 
    %             phase1a = ss(A1a, B1a, zeros(size(A1a)), zeros(size(B1a)));
    %             t = linspace(0,550e-9,5000);
    %             us = repmat(u,1,5000);
    %             
    %             [y,t,x] = lsim(phase1a, us, t, X0);
    %             plot(t, x)
    %             legend('vc1','vc2','vc3','vout','il')


    Cig1a = ic3;
    Dig1a = [1/(5*ron + ESR1 + ESR2 + ESR3), 0];

    ic3 = [1, 1, -1, 0, -2*ron + ESR1]/(5*ron + ESR1 + ESR2 + ESR3);
    ic1 = -ic3 + [0 0 0 0 -1]; 
    ic2 = -ic3;

    A2a = cat(1, ...
    ic1, ...
    ic2, ...
    ic3, ...
    [0, 0, 0, 0, 1], ...
    ic1*(2*ron + ESR1) + [1 0 0 -1 -Rl] );

    B2a = cat(1, ...
    [0, 0], ...
    [0, 0], ...
    [0, 0], ...
    [0, -1], ...
    [0, 0]);


    Cig2a = [0 0 0 0 0];
    Dig2a = [0 0];

    %             X0 = [6.533, 11.729, 17.8777, 4.8924, 1.378];
    %             vl = A2a(end,:)*X0'+B2a(end,:)*u
    % 
    %             stateElem = [C1 C2 C3 Cout L];
    %             K = diag(stateElem);
    %             A2a = K^-1*A2a;
    %             B2a = K^-1*B2a;
    % 
    %             phase1a = ss(A2a, B2a, zeros(size(A1a)), zeros(size(B1a)));
    %             t = linspace(0,550e-9,5000);
    %             us = repmat(u,1,5000);
    %             
    %             [y,t,x] = lsim(phase1a, us, t, X0);
    %             plot(t, x)
    %             legend('vc1','vc2','vc3','vout','il')



    A1b = [0, 0, 0, 0, 1;
      0, 0, 0, 0, -1;
      0, 0, 0, 0, 0;
      0, 0, 0, 0, 1;
      -1, 1, 0, -1, (-3*ron-Rl-ESR1-ESR2)];

    B1b = [0, 0;
    0, 0;
    0, 0;
    0, -1;
    0, 0];

    Cig1b = [0 0 0 0 0];
    Dig1b = [0 0];

    A2b = [0, 0, 0, 0, 0;
      0, 0, 0, 0, 1;
      0, 0, 0, 0, -1;
      0, 0, 0, 0, 1;
      0, -1, 1, -1, (-3*ron-Rl-ESR3-ESR2)];

    B2b = [0, 0;
    0, 0;
    0, 0;
    0, -1;
    0, 0];

    Cig2b = [0 0 0 0 0];
    Dig2b = [0 0];

    A3 = [0, 0, 0, 0, 0;
      0, 0, 0, 0, 0;
      0, 0, 0, 0, 0;
      0, 0, 0, 0, 1;
      0, 0, 0, -1, -(ron)-Rl];

    B3 = [0, 0;
    0, 0;
    0, 0;
    0, -1;
    0, 0];

    Cig3 = [0, 0,0, 0,0];
    Dig3 = [0, 0];

    stateElem = [C1 C2 C3 Cout L];
    K = diag(stateElem);


    M = 5/6;
    phi = 0;%-M*0*Ts;
    th = 0;%-M*0.05*Ts;
    
    Da = 2/3;

    ts = [Da*Ts/2*M-phi+th, (1-Da)*Ts/2*M+phi+th, Ts/2*(1-M), Da*Ts/2*M-phi-th, (1-Da)*Ts/2*M+phi-th, Ts/2*(1-M)];


    As = cat(3, A1a, A1b, A3, A2a, A2b, A3);
    Bs = cat(3, B1a, B1b, B3, B2a, B2b, B3);
    Cs = cat(3, Cig1a, Cig1b, Cig3, Cig2a, Cig2b, Cig3);
    Ds = cat(3, Dig1a, Dig1b, Dig3, Dig2a, Dig2b, Dig3);


    for i = 1:size(As,3)
        As(:,:,i) = K^-1*As(:,:,i);
        Bs(:,:,i) = K^-1*Bs(:,:,i);
    end
    
    simulator.settopology(As, Bs, Cs, Ds);
    simulator.setmodulation(ts);


    
    tic
    for it2 = 1:length(Iorange)
        Io = Iorange(it2);
        u = [Vg Io]';
        Xi = [Vg/4, Vg/2, 3*Vg/4, 6*M Io];
%         Xi = [2 6 9 6*M/2 Io];
        Bi = [.9 .9 .9 .9 .9];
        
        simulator.setinputs(u);
        
        Xs = simulator.SS_Soln( Xi, Bi);
        
%         [ Xs] = SS_Soln( As, Bs, ts, u, Xi, Bi);

        [ xs, t, y ] = simulator.SS_WF_Reconstruct( );

%         [ xs, t, y ] = SS_WF_Reconstruct( Xs, As, Bs, ts, u, Cig, Dig );

        ylabels = {'V_{C1}', 'V_{C2}', 'V_{C3}', 'V_{out}', 'I_L'};
        if(plotIterWF)
            figure(1)
            for i=1:5
                subplot(50,1,i*10-9:i*10)
                hold on;

                plot(t/Ts,xs(i,:), 'Linewidth', 3); 
                ylabel(ylabels{i})
                box on

                if(i<5)
                    set(gca, 'Xticklabel', []);
                else
                    xlabel('t/T_s')
                end
            end
            drawnow;
        end


        Vout(it, it2) = trapz(t,xs(4,:))/Ts;
        Pout(it, it2) = Vout(it, it2)*Io;
        Pin(it, it2) = trapz(t,y*Vg)/Ts + (.5*stateElem*xs(:,1).^2 - .5*stateElem*xs(:,end).^2)/Ts;
        
            
        Pcoss = 2*(1/2)*(2*V)^2*Coss*fs + 5*(1/2)*(V)^2*Coss*fs;
        Pg = Qg*5*fs*8;

        eta(it, it2) = Pout(it, it2)/(Pin(it, it2)+Pcoss+Pg);
        Psw(it, it2) = Pcoss + Pg;
        Ploss(it, it2) = Pin(it, it2)+Psw(it, it2) - Pout(it, it2);
    end
    toc


    figure(2); hold on;
    plot( Pout(it,:),eta(it,:), 'b','LineWidth', 2)
    ylim([.8 1])
    ylabel('\eta [%]')
    xlabel('P_o [W]')
    box on
    hline = findobj(gcf, 'type', 'line');
    
    plot(MaeveData(:,1), MaeveData(:,2),'or','LineWidth', 2)
    
    params_store(it,:) = [Rl Qg Coss ron L C1 fs ESR1 ESR2 ESR3];
    % drawnow
end


%%
figure(7);
    for i = 1:2:size(Ploss,1)-1
        if i < length(partials(:))/2
            sty = '-';
        else
            sty = '-.';
        end
%         xi = length(partials(:))/2;
%         DpDx = diff(Ploss(i+xi:xi+i+1, :))/diff(partials(:,(i+1)/2));
%         DpDx = diff(Ploss(i+xi:xi+i+1, :))/.02;
        DpDx = (Ploss(i+1,:) - Ploss(i,:))./Pout(i,:)/.02;
        plot(Iorange, (DpDx), sty, 'LineWidth',3);
        hold on;
    end

% [Rl Qg Coss ron L C1 fs ESR1 ESR2 ESR3];
legend('R_L', 'Q_g', 'C_{oss}', 'r_{on}', 'L', 'C_{fly}', 'f_s', 'ESR1', 'ESR2', 'ESR3')
% ylim([-10 10])
xlabel('I_{out}');
ylabel('\partial (P_{loss}/P_{out}) / \partial (x/x_{nom})')


% plot( Iorange, 1-Psw./Pin)