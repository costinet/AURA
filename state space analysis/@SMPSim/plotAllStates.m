function plotAllStates(obj, fn, subplots)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    [ xs, t] = obj.SS_WF_Reconstruct;
    figure(fn);
    ns = size(xs,1);
    
    if(nargin == 2)
        subplots = 1;
    end
    
    if(subplots)
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
    else
        lines = findobj(gca,'Type','Line');
        firstrun = (numel(lines) == 0);
        for i = 1:numel(lines)
          lines(i).LineWidth = 1;
          lines(i).LineStyle = ':';
        end

        plot(t,xs, 'linewidth',2);
        if(firstrun)
            legend(obj.converter.topology.stateLabels);
        end
        hold on;
%         plot(t, obj.converter.topology.constraints.regtarget*ones(size(t)), '-.k', 'linewidth',3);
        ylims = [min(min(xs)) max(max(xs))];
        ylim(ylims)
        for i = 1:length(obj.ts)
            plot(sum(obj.ts(1:i))*ones(1,2), ylims, ':r');
        end
    end

end