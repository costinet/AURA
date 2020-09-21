% time vector
% voltage vector
% currnt vector
% frequency
% gOpt... 1 = graph, 0 = no graph
% figNum = integer: specify the number of the figure to graph

function [phi,Vfund,Cfund,THDV,THDC] = wavesPhi(time,V,C,f,gOpt,figNum)

    % A + B parts of Fourier Series for waveform "V"
    anv = (2/(time(end)-time(1)))*trapz(time, V.*cos(2*pi*f*time)); 
    bnv = (2/(time(end)-time(1)))*trapz(time, V.*sin(2*pi*f*time)); 
    Vfund = sqrt(anv^2 + bnv^2); 
    phiV = atan2(anv, bnv)*180/pi; 

    % A + B parts of Fourier Series for waveform "C"
    anc = (2/(time(end)-time(1)))*trapz(time, C.*cos(2*pi*f*time)); 
    bnc = (2/(time(end)-time(1)))*trapz(time, C.*sin(2*pi*f*time));  
    Cfund = sqrt(anc^2 + bnc^2); 
    phiC = atan2(anc, bnc)*180/pi; 

    % Phase between "V" and "C" in [degrees]
    phi = phiV - phiC;
    if(phi > 180) phi = phi - 360; 
    elseif(phi < -180) phi = phi + 360; 
    end 

    % THD in percentage
    RMSV = sqrt(mean(V.*V)); 
    THDV = 100*sqrt(RMSV^2-(Vfund/sqrt(2))^2)/(Vfund/sqrt(2)); 
    
    % THD in percentage
    RMSC = sqrt(mean(C.*C)); 
    THDC = 100*sqrt(RMSC^2-(Cfund/sqrt(2))^2)/(Cfund/sqrt(2)); 

    if(gOpt == 1)
        figure(figNum); 
        plot(time, V);
        grid on;
        yyaxis right; 
        plot(time, C);
        set(gca, 'XTick', []);  % Supress X axis tick marks
        ylabel('Current [A]', 'fontSize', 20); yyaxis left; 
        ylabel('Voltage [V]', 'fontSize', 20);
    end; 
    
end