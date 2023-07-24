function [ Ploss, Ion_s, Ioff_s, Pcond, Pconst, Poss, Esw] = EPC9118WhiteBox_DPT( param, u, Buck, prewarp)

plotWFs = 0;

% display(param);

[cs, ks] = getCKs(Buck, param, prewarp);

Ploss = zeros(1,size(u,1));
Ion_s = zeros(1,size(u,1));
Ioff_s = zeros(1,size(u,1));
Pcond = zeros(1,size(u,1));
Pconst = zeros(1,size(u,1));
Poss = zeros(1,size(u,1));
Esw = zeros(1,size(u,1));

for i = 1:size(u,1);
    Vg = u(i,1);
    Io = u(i,2);
    [Xss, ys, t, ts] = simulate(Buck, Io, Vg);
    [ Pcond, Pconst, Poss, Esw, Ion, Ioff] = FindEsw(Buck, ys, t, ts, Xss, Vg, Io, cs, ks);

    Ploss(i) =  Pcond + Pconst + Poss + Esw*Buck.fs*1e-6;
    Ion_s(i) = Ion;
    Ioff_s(i) = Ioff;
end

Ploss = Ploss';
% display(max(Ploss));
% display(min(Ploss));


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

