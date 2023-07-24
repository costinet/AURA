function [ Pcond, Pg, Pbd, Poss, Pq, Pov, Pboot, Pcore, Voavg, Temps] = EPC9118WhiteBox( Io, Vg, Buck )

plotWFs = 0;

if(nargin == 0)
    Io = 12.9987;%1.6993;
    Vg = 34.0558;%38.0912;
end

[Xss, ys, t, ts] = simulate(Buck, Io, Vg);
[ Pcond, Pg, Pbd, Poss, Pq, Pov, Pboot, Pcore, Voavg, Temps] = CalculateLosses(Buck, ys, t, ts, Xss, Vg, Io);

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

