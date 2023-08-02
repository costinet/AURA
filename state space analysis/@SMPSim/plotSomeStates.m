function plotSomeStates(obj, fn, states ,subplots)
%PLOTSOMESTATES is a SMPSim function that plots some of the state variables across the entire period
%   plotSomeStates(plot#,[state#]) plots all of the state variables as subplots under the plot
%   number specified
%
%   See also PLOTSOMEOUTPUTS

[ xs, t] = obj.SS_WF_Reconstruct;
figure(fn);
ns = length(states);

if(nargin == 3)
    subplots = 1;
end

counter = 1;

if(subplots)
    nrMax = 10;
    nGrid = 10;
    nc = ceil(ns/nrMax);
    nr = min(ns,nrMax);
    for j = 1:nc
        for i=1:nr
            nsi = states(counter);

            inds = [((i-1)*nGrid+1):(i*nGrid)]*nc - (nc-1) + (j-1);
            subplot(nGrid*nr,nc,inds)
            plot(t,xs(nsi,:), 'Linewidth', 3);
            hold on;
            ylims = ylim;
            ylim(ylims);
            xlim([min(t) max(t)]);
            for k = 1:length(obj.ts)
                plot(sum(obj.ts(1:k))*ones(1,2), ylims, ':k');
            end
            hold off;
            try
                ylabel(obj.stateNames{nsi});
            catch
                warning('State Names not set in topology subclass');
            end
            box on
            if(i<nrMax) && states(counter) ~= states(end)
                set(gca, 'Xticklabel', []);
            else
                xlabel('t [\mus]')
            end
            
            if states(counter) == states(end)
                break
            end
            counter = counter+1;
        end
        if states(counter) == states(end)
            break
        end
    end
    %         for i=1:ns
    %             subplot(10*ns,1,i*10-9:i*10)
    %             plot(t,xs(i,:), 'Linewidth', 3);
    %             hold on;
    %             ylims = ylim;
    %             ylim(ylims)
    %             for j = 1:length(obj.ts)
    %                 plot(sum(obj.ts(1:j))*ones(1,2), ylims, ':k');
    %             end
    %             hold off;
    %             ylabel(obj.getstatenames{i});
    %             box on
    %             if(i<ns)
    %                 set(gca, 'Xticklabel', []);
    %             else
    %                 xlabel('t')
    %             end
    %         end
else
    lines = findobj(gca,'Type','Line');
    firstrun = (numel(lines) == 0);
    for i = 1:numel(lines)
        lines(i).LineWidth = 1;
        lines(i).LineStyle = ':';
    end

    plot(t,xs, 'linewidth',2);
    if(firstrun)
        legend(obj.stateNames);
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