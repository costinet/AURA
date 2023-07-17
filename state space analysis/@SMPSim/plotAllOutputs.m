function plotAllOutputs(obj, fn, oSelect, subplots)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    [ ~, t, ys] = obj.SS_WF_Reconstruct;
    fig = figure(fn);
    ns = size(ys,1);
    
    if(nargin <= 2)
        subplots = 1;
        oSelect = 1:ns; 
    elseif(nargin <=3)
        subplots = 1;
    end
    
    ns = length(oSelect);
    
    if(subplots)
        for i=1:ns
            subplot(10*ns,1,i*10-9:i*10)
            plot(t,ys(oSelect(i),:), 'Linewidth', 3);
            hold on;
            ylims = ylim;
            ylim(ylims)
            for j = 1:length(obj.ts)
                plot(sum(obj.ts(1:j))*ones(1,2), ylims, ':k');
            end
            hold off;
            try
                ylabel(obj.converter.topology.outputLabels{oSelect(i)});
            catch
                warning('Output Labels not set in topology subclass');
            end
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

        plot(t,ys(oSelect,:), 'linewidth',2);
        if(firstrun)
            legend(obj.converter.topology.outputLabels);
        end
        hold on;
%         plot(t, obj.converter.topology.constraints.regtarget*ones(size(t)), '-.k', 'linewidth',3);
        ylims = [min(min(ys(oSelect,:))) max(max(ys(oSelect,:)))];
        ylim(ylims)
        for i = 1:length(obj.ts)
            plot(sum(obj.ts(1:i))*ones(1,2), ylims, ':r');
        end
    end

    %% Link all time axes
    children = get(fig,'Children');
    ax = [];
    for i = 1:length(children)
        if strcmp(class(children(i)),'matlab.graphics.axis.Axes')
            ax = [ax, children(i)];
        end
    end
    linkaxes(ax,'x');

end