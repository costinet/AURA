classdef SMPSconverter < handle
    %Converter object for use with AURA
    %   SMPSconverter is contained within the SMPSim class.  SMPSconverter
    %   contains methods and handles information related to the modulation
    %   pattern of the converter.  In most cases, SMPSconverter is
    %   automatically constructed when SMPSim() is constructed.
    %
    %   See also SMPSim, SMPStopology, AURAdb
    
    
    properties
        topology
        u
        
%         tcomp % compensation indexes. When ts(i) decreases by dt, ts(tcomp(i)) should increase by dt
%         tsmax % max permissible values for ts 

    end
    
    properties (Dependent = true)
        As
        Bs
        Cs
        Ds
        Is
        swseq % as integer indexes
        swvec % as binary switch states
        
        ts
        swind %indices into the 3rd dimension of A/B/C/D etc. matrices stored
        % by the topology that form the switching pattern. e.g. swseq = [1 3] 
        % will form a period using As(:,:,1) and As(:,:,3)
    end
    
    properties (SetAccess = private)
        expAts
        cachets
        cacheAs
    end
    
    properties (SetAccess = protected, Hidden)
        fullts
        fullswind
        controlledts
        previousts
    end
    
    properties (Hidden)
        caching = 0
        timingThreshold = 1e-15;

        maximumTimeSteps = 150;
        
    end
    
    methods
        %% Externally Defined
        [altered, newSwInd] = addUncontrolledSwitching(obj, interval, beforeAfter, initialTime, switches, newStates, force)
        [tps] = validateTimePerturbations2(obj, tis, dts, lim)
        
%         function dts = adjustUncontrolledTimingVector(obj, ti, dt)
%             
%             %%%%%%%%%%%%%%%%%%%%
%             
% % %             newAttempt
%                     dts = zeros(size(obj.ts));
%                     tis = 1:length(dts);
%                     
%                     tcomps = dts;
%                     dts(ti) = dt;
% 
%                     [~, ints, subInts] = getIntervalts(obj);
%                     numSubInts = histcounts(ints, max(ints));
%                     
%                     dts(numSubInts == 1) = 0;
%                     nextint = circshift(ints,-1);
%                     swapsies = ((ti' == length(ints)) | (nextint(ti) ~= ints(ti))) ...
%                         & dts' ~= 0;
%                     normalInts = ((ti' ~= length(ints)) & (nextint(ti) == ints(ti))) ...
%                         & dts' ~= 0;
%                     tcomps(find(swapsies)-1) = -dts(swapsies);
%                     tcomps(find(normalInts)+1) = tcomps(normalInts)-dts(normalInts);
%                     
%                     dts = dts + tcomps;
%                     
%                     assert(sum(dts) ==0, 'Non-zero net time perturbation');
% 
%                     [~,maxDts] = getDeltaT(obj);
% 
%                     if any(dts> maxDts)
%                         dts = dts./max(dts./maxDts);
%                     end
%         %             if any(obj.ts + dt < 0)
%                     [dts] = validateTimePerturbations(obj, tis, dts);
% 
%                     obj.fullts(subInts(tis),ints(tis)) = obj.fullts(subInts(tis),ints(tis))+dts;
% 
%                     obj.eliminateZeroTimeIntervals;
%                     
% 
%         end 
            %%%%%%%%%%%%%%%%%%%%
        function [dt, recurs] = adjustUncontrolledTiming(obj, ti, dt)
            obj.previousts = obj.fullts;    %Save old times in case an undo is needed.
            prevDT = dt;
            recurs = 0;

            if numel(ti) == length(obj.ts)
                [tps] = obj.validateTimePerturbations2(dt);
                if sum(tps-dt) < eps*numel(dt)
                    %They're already formatted correctly
                    [~, ints, subInts] = getIntervalts(obj);
                    locs = [sub2ind(size(obj.fullts),subInts,ints)]';
                    obj.fullts(locs) = obj.fullts(locs) + dt;

                end
            

            elseif numel(ti) > 1  
                for i=length(dt):-1:1
                    if dt(i) ~= 0
                        obj.adjustUncontrolledTiming(ti(i), dt(i));
                    end
                end
            

            elseif numel(ti) == 1
                [ots, ints, subInts] = getIntervalts(obj);
                numSubInts = sum(ints == ints(ti));
                
                if (numSubInts==1)
    %                 warning('Cannot adjust timing from controlled intervals without inserting an additional state');
                    dt = 0;
                    return
                elseif ti == length(ints) || ints(ti+1) ~= ints(ti)
                    %selected last interval
                    dt = adjustUncontrolledTiming(obj, ti-1, -dt);
                    recurs = 1;
                    return
                    
                else
                     if dt > 0
                        dt = min([dt, ots(ti+1)]);
                    else
                        dt = max([dt, -ots(ti)]);
                    end
    
                    obj.fullts(subInts(ti),ints(ti)) = obj.fullts(subInts(ti),ints(ti))+dt;
                    obj.fullts(subInts(ti)+1,ints(ti)) = obj.fullts(subInts(ti)+1,ints(ti))-dt;
                end
            end
            
            obj.eliminateZeroTimeIntervals;
 
%             try
                assert(all(abs((sum(obj.fullts,1) - obj.controlledts))<obj.timingThreshold), ...
                    'Timing mismatch due to uncontrolled switching');
%             catch
%                 warning('Timing mismatch due to uncontrolled switching');
%             end
        end
        
        function [that, dts] = getDeltaT(obj, ind)
            if nargin <2
                ind = 1:length(obj.swind);
            end
            
            fastestResFreq = zeros(length(ind),1);
            fastestFreq = fastestResFreq;
            for i = 1:length(ind)
                eigA = obj.topology.eigAs(:,:,obj.swind(i));%eigs(obj.As(:,:,i));
                fastestResFreq(i) =  max([1; abs(eigA(imag(eigA)>0))]);
                fastestFreq(i) = max(abs(eigA));
            end


            %% that = good approx for small-signal perturbations
%             that = min( [2*pi./fastestFreq/4, obj.ts'/100], [],2);             
            if size(fastestFreq) == size(obj.ts)
                that = min( [2*pi./fastestFreq/4, obj.ts(ind)/100], [],2); 
            else
                that = min( [2*pi./fastestFreq/4, obj.ts(ind)'/100], [],2); 
            end
            that = reshape(that, 1, length(that));
            
            %% dts = allowable large-signal perturbation to avoid introducing error between discrete samples
%             dts =  min( [2*pi./fastestResFreq/4, obj.ts'], [],2); 
            if size(fastestResFreq) == size(obj.ts)
                %Added the fastestFreq portion due to two dueling
                %exponentials causing a change in slope
                dts = min( [2*pi./fastestResFreq/4, obj.ts(ind)], [],2); 
            else
                dts = min( [2*pi./fastestResFreq/4, obj.ts(ind)'], [],2); 
            end
            dts = reshape(dts, 1, length(dts));
            
                             
        end
        
        
        function [tps] = validateTimePerturbations(obj, tis, dts, lim)
            assert(length(tis) == length(dts), 'vectors tis and dts must have the same number of elements');
            scaleFull = 1;          %scales the time perturbations to keep equal change vector
            individualLimit = 0;    %tries to only reduce the perturbations that violate (needs work)
            
            if nargin == 3
                lim = 0;
            end
            
            odts = dts;

            if lim == 1
                % Scale all values to make sure we're below 
                % a reasonable approximation to deltaT
                [~,dtLims] = getDeltaT(obj);
        
                tr = dts./dtLims;
                tr = max(tr(tr>1));
                
                if ~isempty(tr)
                    dts = dts/tr;
                end
            end
            
            [ots, ints, ~] = getIntervalts(obj);
            
            % Remove any zero intervals
            tis(dts == 0) = []; dts(dts == 0) = [];
            
            %% ???
            [tis, inds] = sort(tis);
            dts = dts(inds);
            
            % deal with controlled switching timing
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
            
            dts(dts > 0) = min([dts(dts > 0), ots(tis(dts > 0)+1)]);
            dts(dts < 0) = max([dts(dts < 0), -ots(tis(dts < 0))]);
            
            if(scaleFull)
                tps = zeros(1,length(ots));
                tps(tis)= dts;
                tps(tis+1) = tps(tis+1)-dts;
                    
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
            elseif(individualLimit)
                 nonNegTimes = 0;
                while nonNegTimes == 0

                    tps = zeros(1,length(ots));
                    tps(tis)= dts;
                    tps(tis+1) = tps(tis+1)-dts;

                    %% because multiple pertubations, still possible to go below zero
                    errLocs = ots + tps < 0;
                    errInds = find(errLocs);

                    chLoc = errLocs | circshift(errLocs,-1);
                    chLoc = intersect(find(chLoc), tis);
                    dtsLoc = ismember(tis,chLoc);
                    dts(dtsLoc) = -(ots(chLoc) + tps(chLoc) + eps) + dts(dtsLoc);

                    tps = zeros(1,length(ots));
                    tps(tis)= dts;
                    tps(tis+1) = tps(tis+1)-dts;

                    if(any(ots + tps< 0))
                        nonNegTimes = 0;
                    else
                        nonNegTimes = 1;
                    end
                end
            end
            
            assert(all(ots + tps >= 0), 'something went wrong, a time is negative');
            assert(sum(tps) == 0, 'something went wrong, time perturbations sum to non-zero value');
        end
        
        function eliminateZeroTimeIntervals(obj)
%             fts = obj.fullts;
%             dts = diff(fts,1);
%             [r,c] = find(fts(2:end,:) == dts & fts(2:end,:) ~= 0);
%             [r,c] = find(obj.fullts == 0 & obj.fullswind ~= 0);
            [r,c] = find(obj.fullts < obj.timingThreshold*0.9 & obj.fullswind ~= 0);
            
            for i = length(r):-1:1
                if obj.fullts(r(i),c(i)) > eps/100
                    % if it isn't actually zero, set it to zero and keep
                    % from altering the overall controlled interval length.
                    obj.adjustUncontrolledTiming(find(obj.ts == obj.fullts(r(i),c(i))), -obj.fullts(r(i),c(i)));
                    %NOTE: adjustUncontrolledTiming will re-call this
                    %function, so we can end here.
                else
                    obj.fullts(r(i),c(i)) = 0;
                    obj.fullswind(r(i),c(i)) = 0;
                    
                    if r(i) < size(obj.fullts,1) && obj.fullswind(r(i)+1,c(i)) ~= 0
                	    obj.fullts(:,c(i)) = ...
                            [obj.fullts(1:r(i)-1,c(i));
                             obj.fullts(r(i)+1:end,c(i));
                             0];     
                        obj.fullswind(:,c(i)) = ...
                            [obj.fullswind(1:r(i)-1,c(i));
                             obj.fullswind(r(i)+1:end,c(i));
                             0];   
                    end
                end
            end

            allZeroRows = all(obj.fullts ==0,2);
            allZeroCols = all(obj.fullts ==0,1);

            obj.fullts(allZeroRows,:) = [];
            obj.fullts(:,allZeroCols) = [];

            obj.fullswind(allZeroRows,:) = [];
            obj.fullswind(:,allZeroCols) = [];
            
%             if ~isempty(r)
%                 obj.fullts(r,c) = -1;
%                 inds = find(obj.ts == -1);
%                 obj.swind(inds) = [];
%                 
%                 for i = length(r):-1:1
%                     obj.fullts(:,c(i)) = ...
%                         [obj.fullts(1:r(i)-1,c(i));
%                          obj.fullts(r(i)+1:end,c(i));
%                          0];   
%                 end
%             end
        end

        function eliminateRedundtantTimeIntervals(obj)
            if size(obj.fullswind,1) == 1
                return
            end
            locs = diff(obj.fullswind) == 0 & obj.fullswind(2:end,:) ~= 0;
            locs = [zeros(1,size(locs,2)) == 1; locs];
            dts = obj.fullts(locs);
            
            obj.fullts(locs) = 0;
            obj.fullts(circshift(locs,-1)) = obj.fullts(circshift(locs,-1)) + dts;
            obj.eliminateZeroTimeIntervals;

        end
        
%         function undoLatestTimeChange(obj)
%             error("This isn't implemented")
%            obj.fullts = obj.previousts; 
%         end
        
        function [ts, ints, subInts] = getIntervalts(obj)
            % ts -- times
            % ints -- which controlled switcihng interval
            % subints -- which subinterval within parent controlled sw int.
            ts = obj.ts;
            [r,c] = find(obj.fullts ~= 0);
            ints = c;
            subInts = r;
            
            %Make sure they are tall vectors
            ints = reshape(ints,length(ints),1);
            subInts = reshape(subInts,length(subInts),1);
        end
        
        function [TF,lastInt,Ihc] = checkForSymmetry(obj) 
            lastInt = [];
            Ihc = [];
            TF = 0;
            if mod(length(obj.ts),2) == 0
                halfPoint = length(obj.ts)/2;
%                 Mat1 = reshape(cat(2, obj.As(:,:,1:halfPoint), obj.Bs(:,:,1:halfPoint)), size(obj.As(:,:,1:halfPoint),1),...
%                     (size(obj.As(:,:,1:halfPoint),2) + size(obj.Bs(:,:,1:halfPoint),2))*size(obj.As(:,:,1:halfPoint),3) );
%                 Mat2 = reshape(cat(2, obj.As(:,:,halfPoint+1:end), obj.Bs(:,:,halfPoint+1:end)), size(obj.As(:,:,halfPoint+1:end),1),...
%                     (size(obj.As(:,:,halfPoint+1:end),2) + size(obj.Bs(:,:,halfPoint+1:end),2))*size(obj.As(:,:,halfPoint+1:end),3) );
%                 if all(abs(Mat1) == abs(Mat2))
%                     TF = 1;
%                     lastInt = halfPoint;
%                     Ihc = obj.As(:,:,halfPoint+1)\obj.As(:,:,1);
%                     error('Ongoing Development Here')
%                 else
%                     TF = 0;
%                 end
                Atil = [obj.As, obj.Bs];
                Atil = [Atil; zeros(size(Atil,2)-size(Atil,1), size(Atil,2), size(Atil,3))];
                EA = eye(size(Atil,1));
                EA2 = EA;
                for i = 1:halfPoint
                    EA = expm(Atil(:,:,halfPoint-i+1))*EA;
                    EA2 = expm(Atil(:,:,end-i+1))*EA2;
                end
                if norm(abs(EA) - abs(EA2),2) < eps
                    TF=1;
                    lastInt = halfPoint;

                    signs = EA./EA2;
                    Ihc = diag(signs(1:size(obj.As,1),size(obj.As,1)+1));
                end
            end
        end
        
        
%         function [tnew] = adjustTiming(obj, ti, dt)
%             tnew = obj.ts;
%             
% %             if( dt + tnew(ti) < 0)
% %                 dt = -tnew(ti);
% %             end
% %             if( dt + tnew(ti) > obj.tsmax(ti))
% %                 dt = obj.tsmax(ti) - tnew(ti);
% %             end
% %             if( tnew(obj.tcomp(ti)) - dt > obj.tsmax(obj.tcomp(ti)))
% %                 dt = -obj.tsmax(obj.tcomp(ti)) + tnew(obj.tcomp(ti));
% %             end
%             
%             if dt > 0
%                 dt = min([dt, obj.ts(obj.tcomp(ti)), obj.tsmax(ti)-obj.ts(ti)]);
%             else
%                 dt = max([dt, -obj.ts(ti), -(obj.tsmax(ti)-obj.ts(obj.tcomp(ti)))]);
%             end
%            
%             tnew(ti) = tnew(ti) + dt;
%             tnew(obj.tcomp(ti)) = tnew(obj.tcomp(ti)) - dt;
%             
%             assert(~(sum(tnew < 0) && sum(tnew>obj.tsmax)), 'Timing adjustment failed');
%             
%             obj.ts = tnew;
%             obj.refreshCache();
%         end
        

        function eigs2tis(obj)
            for i = length(obj.ts):-1:1
                eigA = obj.topology.eigAs(:,:,obj.swind(i));%eigs(obj.As(:,:,i));
                fastestResFreq =  max([1; abs(eigA(imag(eigA)>0))]);
                dt =  min( [2*pi./fastestResFreq/8, obj.ts(i)/2'], [],2);
                
                ratio = ceil(obj.ts(i)/dt - eps);
                
                for j = 1:1:ratio-1
                    obj.addUncontrolledSwitching(i, 1, dt, [], [], 1);
                end
            end
        end

%         function refreshCache(obj)
%             if obj.caching
%                 if(size(obj.swseq) == size(obj.As,3))
%                     for i = 1:size(obj.As,3)
%                         obj.expAts(:,:,i) = expm(obj.As(:,:,obj.swseq(i))*obj.ts(i));
%                     end 
%                 end
%                 if isempty(obj.tsmax)
%                     obj.tsmax = sum(obj.ts)*ones(size(obj.ts));
%                 end
%                 if isempty(obj.tcomp)
%                     obj.tcomp = circshift(1:length(obj.ts),-1);
%                 end
%             end
%         end
        
        function setSwitchingPattern(obj, swind, ts)  
        % SETSWITCHINGPATTERN set switching pattern and controlled timing
        %
        %   setSwitchingPattern(obj, swind, ts)  
        %
        %   setSwitchingPattern(obj, swvec, ts)  
        %
        %
            if size(ts,2) == 1
                ts = ts';
            end
            if size(swind,2) == 1
                swind = swind';
            end
            if numel(swind) ~= length(swind) && all(swind == 1 | swind == 0, 'all')
                %2D swind, assume swseq
                swseqIn = swind;
                swind = zeros(1,size(swseqIn,1));
                for i = 1:size(swseqIn,1)
                    [~, idx] = intersect(obj.topology.swseq, swseqIn(i,:), 'rows');
                    if ~isempty(idx)
                        swind(i) = idx;
                    else
                        %If the switching state doesn't already exist,
                        %parse it
                        obj.topology.loadCircuit([],swseqIn(i,:));
                        [~, idx] = intersect(obj.topology.swseq, swseqIn(i,:), 'rows');
                        swind(i) = idx;
                    end
                end
            end
            
            assert(~any(sum(ts,1)==0), 'Cannot have zero-time intervals');
            
            assert(all(size(swind) == size(ts)), 'Dimensions of swind and ts must be the same');
            assert(all(ts>=0, 'all'), 'times in vector ts must all be positive');
            assert(all(swind >= 0 & swind <= size(obj.topology.As,3), 'all'), 'swind must be a vector of values indexing topology state matrices');          

            obj.controlledts = sum(ts,1);
            obj.fullts = ts;
            obj.fullswind = swind;

            obj.limitTimeIntervals();
        end
        

        %% Getters
        function res = get.As(obj)
            res = obj.topology.As(:,:,obj.swseq);
        end
        
        function res = get.Bs(obj)
            res = obj.topology.Bs(:,:,obj.swseq);
        end
        
        function res = get.Cs(obj)
            res = obj.topology.Cs(:,:,obj.swseq);
        end
        
        function res = get.Ds(obj)
            res = obj.topology.Ds(:,:,obj.swseq);
        end
        
        function res = get.Is(obj)
            res = obj.topology.Is(:,:,obj.swseq);      
        end
        
        function res = get.ts(obj)
            res = obj.fullts(obj.fullswind ~=0)';
            if size(res,1) > 1
                res = res';
            end
        end
        
        function res = get.swind(obj)
             res = obj.fullswind(obj.fullswind ~=0)';
             res = reshape(res,1,length(res));
        end
        
        function res = get.swseq(obj)
            if isempty(obj.swind)
                res = 1:size(obj.topology.As,3);
            else
                res = obj.swind;
            end
        end
        
        function res = get.swvec(obj)
            if isempty(obj.swind)
                res = obj.topology.swseq(1:size(obj.topology.As,3),:);
            else
                res = obj.topology.swseq(obj.swind,:);
            end
        end
        
        function expAts = get.expAts(obj)
            expAts = obj.expAts;
        end
        
        %% Setters
        function set.swseq(obj,newSeq)
%             obj.swind = newSeq;
            [~, index]=ismember(newSeq,obj.topology.swseq,'rows');
            obj.swind = index;
%             obj.refreshCache();
        end
        
        function set.ts(obj, newts)
            error('direct setting of ts not allowed.  Use SMPSconverter.setSwitchingPattern() instead');
%             obj.ts = newts;
%             obj.refreshCache();
        end
        
        function set.u(obj, newU)
            % enforce it being a tall vector
            if size(newU,1) ~= size(obj.Bs,2)
                if size(newU,1) == 1
                    newU = reshape(newU,[size(obj.Bs,2), 1, size(newU,3)]);
                end
            end
            assert(size(newU,1) == size(obj.Bs,2) && size(newU,2) ==1, 'invalid dimensions for u.  The input vector must be either (Ni x 1 x Nt) or (Ni x 1 x Nt) for Nt time intervals and Ni independent inputs');

            obj.u = newU;%reshape(newU,max(size(newU,1:2)),1, size(newU,3));
        end

        function set.topology(obj,top)
            arguments
                obj SMPSconverter
                top SMPStopology
            end
            obj.topology = top;
            top.linkConverter(obj);
        end
        
    end

    methods (Hidden)

        function limitTimeIntervals(obj)
            assert(sum((obj.fullts >0)+0,'all') <= obj.maximumTimeSteps, ...
                'Number of timing subintervals greater than setting in SMPSconverter.maximumTimeSteps.')
        end

    end
    
end

