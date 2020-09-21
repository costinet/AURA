% void function that plots FFT of waveform

% wave =    array of waveform values
% f =       fundamental frequency
% sam =     sampling frequency
% figNum =  specify the figure number of plot
% xlims =   min and max vector for x-axis (kHz) Ex: [10, 1000]
% ylims =   min and max vector for y-axis (dB) Ex: [-40, 30]


function [hh] = plotFFT(wave,f,fsam,figNum,xlims,ylims)

if(nargin < 6)  % check number of args
    lims = 0;   % decide on graph limits
else 
    lims = 1; 
end; 

    % FIND FFT Spectrum
    L = length(wave);               % Data length.
    MYfft = fft(wave);              % Take FFT.
    P2 = abs(MYfft/L);              % 2-sided spectrum. 
    P1 = P2(1:L/2+1);               % 1-sided spectrum. 
    P1(2:end-1) = 2*P1(2:end-1); 
    
    freq = fsam*(0:(L/2))/L;        % Frequency. 
    
    % PLOT FFT Spectrumfont = 'Times';
    font = 'Times';
    ff = figure(figNum); hold off; kHz = 1e-3; 
    set(ff, 'DefaultTextFontName', font, 'DefaultAxesFontName', font);
    set(ff, 'position', [650.0 500 550 300]);
    semilogx(freq*kHz, abs(P1), 'linewidth', 2); 
%     semilogx(freq*kHz, 10*log10(abs(P1)), 'linewidth', 2); 
    grid on; 
    if(lims)            
        xlim(xlims); ylim(ylims); end; 
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 14);
    xlabel('Frequency [kHz]', 'fontSize', 16);
    ylabel('Magnitude', 'fontSize', 16); 

end