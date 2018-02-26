function plotAllStates(obj, fn)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    [ xs, t, ~] = obj.SS_WF_Reconstruct;
    figure(fn);
    ns = size(xs,1);
    for i=1:ns
        subplot(10*ns,1,i*10-9:i*10)
        plot(t,xs(i,:), 'Linewidth', 3);
        hold on;
        ylims = ylim;
        ylim(ylims)
        for j = 1:length(obj.ts)
            plot(sum(obj.ts(1:j))*ones(1,2), ylims, ':k');
        end
        hold off;
        ylabel(obj.getstatenames{i});
        box on
        if(i<ns)
            set(gca, 'Xticklabel', []);
        else
            xlabel('t')
        end
    end
end