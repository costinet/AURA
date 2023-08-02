function plotAllOutputs(obj, fn, subplots)
%PLOTALLOUTPUTS is a SMPSim function that plots all of the output functions described in output equation y=Cx+Du
%   plotAllOutputs(plot#) plots all of the output functions as subplots under the plot
%   number specified
%
%   See also PLOTALLSTATES


[ ~, t, ys] = obj.SS_WF_Reconstruct;
figure(fn);
ns = size(ys,1);

if(nargin == 2)
    subplots = 1;
end
%ns = 16;

if(subplots)
    for i=1:ns
        subplot(10*ns,1,i*10-9:i*10)
        plot(t,ys(i,:), 'Linewidth', 3);
        hold on;
        ylims = ylim;
        ylim(ylims)
        for j = 1:length(obj.ts)
            plot(sum(obj.ts(1:j))*ones(1,2), ylims, ':k');
        end
        hold off;
        try
            ylabel(obj.converter.topology.outputLabels{i});
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

    plot(t,ys, 'linewidth',2);
    if(firstrun)
        legend(obj.converter.topology.outputLabels);
    end
    hold on;
    %         plot(t, obj.converter.topology.constraints.regtarget*ones(size(t)), '-.k', 'linewidth',3);
    ylims = [min(min(ys)) max(max(ys))];
    ylim(ylims)
    for i = 1:length(obj.ts)
        plot(sum(obj.ts(1:i))*ones(1,2), ylims, ':r');
    end
end

end