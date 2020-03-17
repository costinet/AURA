% This code is very well written and plots some sine waveforms for the
% COMPEL paper in the year MMXX

close all
Omega = linspace(0, 6*pi(), 10000);

t0 = linspace(2,-2,size(Omega,2));

t=sin(Omega);

t1=0.1*cos(Omega*20);
t2=cos(Omega);


figure(1)




%subplot(3,1,1)
hold on
plot(t0,':r','LineWidth',3)

%subplot(3,1,2)
hold on
plot(t,'--b','LineWidth',3)

%subplot(3,1,3)
hold on
plot(t1+t2,'-k','LineWidth',3)
ylim([-2 2])


% hold on
% too = -0.5*ones(size(Omega,2));
% plot(too,'-m','LineWidth',3)


legend({'Linear','Single Frequency','Mutli-Frequncy'},'FontSize',12)

