function plotWaveforms(obj, type, fn, oSelect, subplots)
%plot all states/outputs over one period from the current steady-state solution
%
%   plotWaveforms(obj, type, fn, oSelect, subplots)
%   for SMPSim object, plots all states (type==1) or outputs (type==2) to
%   figure number fn.  oSelect selects only certain signals from the
%   state/output vector.  subplots is a boolean variable to determine 
%   whether to plot each signal on its own subplot or plot them all on a 
%   single plot.
%
%   In most cases, plotAllStates or plotAllOutputs should be used instead
%   of directly using plotWaveforms due to the lack of input validation
%   here.
%
%   See Also SMPSim.SS_WF_Reconstruct,  SMPSim.plotAllStates, SMPSim.plotAllOutputs  

    if type == 1
        [ sigs, t] = obj.SS_WF_Reconstruct;
        names = obj.stateNames;

        if ~all(obj.IHC == eye(size(obj.IHC)),'all')
            sigs = [sigs, sigs.*repmat(diag(obj.IHC),[1,size(sigs,2)])];
            t = [t, t+t(end)];
        end
    elseif type == 2
        [ ~, t, sigs] = obj.SS_WF_Reconstruct;
        names = obj.outputNames;

%         if ~all(obj.IHC == eye(size(obj.IHC)),'all')
%             sigs = [sigs, sigs.*repmat(diag(obj.IHC),[1,size(sigs,2)])];
%             t = [t, t+t(end)];
%         end
    end 
    fig = figure(fn);

    ns = size(sigs,1);
    if isempty(oSelect)
        oSelect = 1:ns;
    else
        ns = length(oSelect);
    end
    
    if(subplots)
        nrMax = 10;
        nGrid = 10;
        nc = ceil(ns/nrMax);
        nr = min(ns,nrMax);
        for j = 1:nc
            for i=1:nr
                nsi = i+nr*(j-1);
                if nsi > ns
                    break
                end
                inds = [((i-1)*nGrid+1):(i*nGrid)]*nc - (nc-1) + (j-1);
                ax = subplot(nGrid*nr,nc,inds, 'Parent', fig);
                plot(ax, t,sigs(oSelect(nsi),:), 'Linewidth', 3);
                hold on;
                ylims = ylim;
                ylim(ylims);
                xlim([min(t) max(t)]);
                for k = 1:length(obj.ts)
                    plot(ax,sum(obj.ts(1:k))*ones(1,2), ylims, ':k');
                end
                hold off;
                try
                    sigName = names{oSelect(nsi)};
                    underInd = strfind(sigName,'_');
                    if ~isempty(underInd)
                        sigName = strrep(sigName, '_', '_{');
                        sigName = [sigName repmat('}',1,numel(underInd))];
                    end
                    ylabel(sigName);
                catch
                    warning('State Names not set in topology subclass');
                end
                box on
                if(i<nr)
                    set(ax, 'Xticklabel', []);
                else
                    xlabel('t')
                end

                %% remove yticks near ends
                closeYTicks = min(abs(ax.YTick - ax.YLim')/diff(ax.YLim),[],1) < 0.1;
                if sum(~closeYTicks) >= 3
                    ax.YTick(closeYTicks) = [];
                else
                    ax.YTick = ax.YLim(1) + diff(ax.YLim)*[.25 .5 .75];
                end


            end
            if nsi > ns
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

        plot(t,sigs(oSelect,:), 'linewidth',2);
        if(firstrun)
            legend(names);
        end
        hold on;
%         plot(t, obj.converter.topology.constraints.regtarget*ones(size(t)), '-.k', 'linewidth',3);
        ylims = [min(min(sigs(oSelect,:))) max(max(sigs(oSelect,:)))];
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
    if numel(ax) >= 2
        linkaxes(ax,'x');
    end
end