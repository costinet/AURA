function [tps] = validateTimePerturbations2(obj, dts)
        assert(length(obj.ts) == length(dts), 'vectors conv.ts and dts must have the same number of elements');
        odts = dts;

        % Scale all values to make sure we're below 
        % a reasonable approximation to deltaT
        [~,dtLims] = getDeltaT(obj);
        dtLims = min(dtLims,obj.ts);

        tr = dts./dtLims;
        tr = max(tr(tr>1));
        
        if ~isempty(tr)
            dts = dts/tr;
        end
        
        
        [ots, ints, subInts] = getIntervalts(obj);

        fullDTs = zeros(size(obj.fullts));
        activeTimeInds = sub2ind(size(fullDTs),subInts,ints);
        fullDTs(activeTimeInds) = dts;

        if all(sum(fullDTs,1)<eps)
            % already formatted appropriately
            % just get it to go to zero exactly
            lastSubInt = max(1,sum((fullDTs ~=0),1));
            inds = sub2ind(size(fullDTs),lastSubInt,1:length(lastSubInt));
            fullDTs(inds) = fullDTs(inds) - sum(fullDTs,1);
            tps = fullDTs(activeTimeInds);
            if size(tps,1) > 1
                tps = tps';
            end
        else
            %pairwise equal perturbationbs
             % deal with controlled switching timing
            tis = 1:length(obj.ts);
            numSubInts = zeros(length(tis),1);
            for i = 1:length(tis)
                numSubInts(i) = sum(ints == ints(tis(i)));
                if tis(i) == length(ints) || ints(tis(i)+1) ~= ints(tis(i))
                    tis(i) = tis(i)-1;
                    dts(i) = -dts(i);
                end
            end

            dts(numSubInts ==1) = 0;
            tis(dts == 0) = []; dts(dts == 0) = [];

            tps = zeros(1,length(ots));
            tps(tis)= dts;
            tps(tis+1) = tps(tis+1)-dts;
        end
        
        % Remove any zero intervals
%             tis(dts == 0) = []; dts(dts == 0) = [];
%      
%             dts(dts > 0) = min([dts(dts > 0), ots(tis(dts > 0)+1)]);
%             dts(dts < 0) = max([dts(dts < 0), -ots(tis(dts < 0))]);
        
        
        
            
        newTs = ots + tps;
        errLocs = newTs< 0;
        if(any(errLocs))
            if all(abs(newTs(errLocs)) < abs(tps(errLocs)))
                overRatio = abs(ots(errLocs)./tps(errLocs));
                tps = tps.*(min(overRatio) - eps);
            else
                error('a time is negative or zero, already');
            end
        end
      
        try
            assert(all(ots + tps >= 0), 'something went wrong, a time is negative');
            assert(abs(sum(tps)) <= eps, 'something went wrong, time perturbations sum to non-zero value');
        catch
            x=1;
        end
    end