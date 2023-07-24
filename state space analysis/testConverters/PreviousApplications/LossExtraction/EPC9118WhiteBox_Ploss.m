function [ Ploss] = EPC9118WhiteBox_Ploss( param, u, Buck )

plotWFs = 0;


Buck.dt = 30e-9*param(1);
Buck.Rl = 1e-3*param(2)*5;

EPC9118_Digitized;
Rja_HsM = 20;
qg_HS = 8.3e-9*param(3);
ron_HS = 5.6e-3;
newCoss = EPC2021CossVds;
newCoss(:,2) = newCoss(:,2)*param(4);
HS_MOS = MOSFET(ron_HS, qg_HS, Rja_HsM, newCoss, EPC2021IsdVsd_T25, EPC2021RonT, 80);
Buck.HS_FET = HS_MOS;


display(param);

for i = 1:size(u,1);
    Vg = u(i,1);
    Io = u(i,2);
    [Xss, ys, t, ts] = simulate(Buck, Io, Vg);
%     try
        [ Pcond, Pg, Pbd, Poss, Pq, Pov, Pboot, Pcore, Voavg, Temps] = CalculateLosses(Buck, ys, t, ts, Xss, Vg, Io);
%     catch
%         x=1;
%     end
    Plosses =  [Pcond; Pg; Pbd; Poss; Pq; Pov; Pboot; Pcore];
    Ploss(i) =  sum(Plosses,1);
end

Ploss = Ploss';


if(plotWFs)
    figure(11);
    for i = 1:size(As,1)
        subplot(size(As,1),1,i)
        plot(t, ys(i,:));
        hold on;
        ylimits = ylim;
        for j = 1:length(ts);
            plot([sum(ts(1:j)) sum(ts(1:j))], [ylimits(1) ylimits(2)], ':r');
        end
%         hold off;
    end
end
end

