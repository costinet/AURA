function [Xf,ts,swinds,us] = timeSteppingPeriod(obj, Xs, ts, origSwind, us)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 

    debug = 0;

    if nargin == 1
        if isempty(obj.Xs)
            obj.steadyState;
        end
        Xs = [obj.Xs(:,1); 1];
        origSwind = obj.converter.swind;
        ts = obj.ts;
        us = obj.fullu;
    end
    [~, deltaTs] = getDeltaT(obj.converter, 1:length(ts));
    if size(deltaTs,2) < size(deltaTs,1)
        deltaTs = deltaTs';
    end
    deltaTs = deltaTs/10;
%     nMinSteps = ceil(sum(obj.converter.fullts,1)./deltaTs);
    
    if size(ts,2) < size(ts,1)
        ts = ts';
    end
    nMinSteps = ceil(ts./deltaTs);
    tVecLength = sum(nMinSteps,'all');
    
    
    subInts = zeros(1,tVecLength);
    tvec = zeros(tVecLength,1);
    
    if debug == 1
        states = zeros(length(Xs),tVecLength);
    end
    
    intus = zeros([size(us,[1:2]),tVecLength]);
    for i = 1:length(origSwind)
        if i == 1
            intus(:,:,1:nMinSteps(1)) = repmat(us(:,:,i),1,nMinSteps(i));
            subInts(1:nMinSteps(1)) = repmat(origSwind(i),1,nMinSteps(i));
            tvec(1:nMinSteps(1)) = linspace(0,ts(1), nMinSteps(1));
        else
            intus(:,:,sum(nMinSteps(1:i-1))+1:sum(nMinSteps(1:i))) = repmat(obj.fullu(:,:,i),1,nMinSteps(i));
            subInts(sum(nMinSteps(1:i-1))+1:sum(nMinSteps(1:i))) = repmat(origSwind(i),1,nMinSteps(i));
            tvec(sum(nMinSteps(1:i-1))+1:sum(nMinSteps(1:i))) = sum(ts(1:i-1)) + linspace(ts(i)/nMinSteps(i),ts(i),nMinSteps(i));
        end
    end
% 
%     xs = zeros(size(obj.Xs,1),tVecLength);
%     xs(:,1) = obj.Xs(:,1);
%     for i = 2:length(obj.converter.ts)
% %         xs(:,i) = xs(:,i-1) + expm( ...
% %             IS it in violation?
% 
%     end



    
    
    %Make sure there is no violation at time zero
%     assert( all( ...
%         obj.converter.topology.Cbnd(:,:,1)*obj.Xs(:,1) + ...
%         obj.converter.topology.Dbnd(:,:,1)*obj.u + ...
%         - obj.converter.topology.bndHyst(:,1) + ...
%         obj.converter.topology.bndHyst(:,2) >= 0 ), ...
%     'Initial states are invalid')
    
    %% NEW Attempt
    try
        i = 1;
        j = 1;
        while i < length(tvec) % or 100 iterations without progress
            deltaT = diff(tvec(i:i+1));
            A = obj.converter.topology.As(:,:,subInts(i));
            B = obj.converter.topology.Bs(:,:,subInts(i));
            Cbnd = obj.converter.topology.Cbnd(:,:,subInts(i));
            Dbnd = obj.converter.topology.Dbnd(:,:,subInts(i));
            u = intus(:,:,i);
            
            M = [A, B*u; ...
                zeros(1,size(A,2) + 1)];
            
            oldXs = Xs;
            Xs = expm(M*deltaT)*Xs;
            
            if debug == 1
                states(:,i) = Xs;
            end
            
            if any(isnan(Xs))
                error('Time stepping resulted in NaN states');
            end

            
%             violateMarginEnd = Cbnd*Xs(1:end-1) + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1) + obj.converter.topology.bndHyst(:,2);
            [violateMarginEnd,targetValEnd] = obj.checkStateValidity(Xs(1:end-1), u, subInts(i));
            [violateMarginStart,targetValStart] = obj.checkStateValidity(oldXs(1:end-1), u, subInts(i));
            
%             if i == 378
%                 x=1;
%             end
            
            if any(violateMarginEnd < 0)
                % There is a violation, but based on how deltaTs was
                % formulated, linear projection to find it should work well
                
                %% Find time at which we hit the error bound
                
%                 violateMarginStart = Cbnd*oldXs(1:end-1) + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1) - obj.converter.topology.bndHyst(:,2);
%                 %Note the swapped error bounds
                
                if any(violateMarginStart < 0) && i == 1
                    %Initial state vector is an error
                    for ij = find(violateMarginStart < 0)'
                        jk = find(Cbnd(ij,:) ~=0,1,'first');
                        Xs(jk) = Xs(jk) - violateMarginStart(ij)./Cbnd(ij,jk)*1.5;
                    end
                    oldXs = Xs;
                    [violateMarginEnd,targetValEnd] = obj.checkStateValidity(Xs(1:end-1), u, subInts(i));
                    [violateMarginStart,targetValStart] = obj.checkStateValidity(oldXs(1:end-1), u, subInts(i));
                    continue
                elseif any(violateMarginStart < 0) && i>1 && (subInts(i) ~=subInts(i-1))
                    % Error at the beginning of a new subinterval
                    tcross = 0;
                    errAt = 0; %beginning
                elseif any((targetValStart) <= 0 & -violateMarginEnd+2*targetValEnd < 0)
                    % beginning point was within the hysteresis bounds               
                    tcross = 0;
                    errAt = 2; %end
                else
                    %crossing occured within the interval
                    crossSlope = ((2*violateMarginStart-violateMarginStart)-violateMarginEnd)/deltaT.*(violateMarginEnd < 0);
                    tcross = targetValStart.*(violateMarginEnd < 0)./ ...
                        (crossSlope);
                    tcross = tcross(violateMarginEnd < 0);
                    tcross = max(min(tcross, deltaT),0);
                    tcross = min(tcross);
                    
                   errAt = 1;
                end
                
                if errAt == 0
                    %% Error at beginning of period
%                     violateMarginStart = Cbnd*oldXs(1:end-1) + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1) + obj.converter.topology.bndHyst(:,2);
%                     FullPeriodViolateStart = zeros(size(violateMarginEnd,1),length(obj.ts));
%                     FullPeriodViolateEnd = zeros(size(violateMarginEnd,1),length(obj.ts));
                    currentInt = find(tvec(i) < cumsum(ts),1,'first');
%                     FullPeriodViolateStart(:, currentInt) = violateMarginStart;
%                     [tLocs,insertAt,adjType] = obj.findRequiredUncontrolledSwitching(FullPeriodViolateStart,FullPeriodViolateEnd)
                    
                    newSwVec = obj.converter.swvec(currentInt,:);
                    modSwitch = obj.topology.switchRef((violateMarginEnd < 0),1);
                    newState = ~obj.topology.switchRef((violateMarginEnd < 0),2);

                    newSwVec(modSwitch) = newState;
                    if ~ismember(newSwVec, obj.topology.swseq, 'rows')
                        obj.topology.loadCircuit([], [obj.converter.swvec; newSwVec], 0);
                    end
                    [~, newSwInd] = ismember(newSwVec, obj.topology.swseq, 'rows'); 
                    
                    %% adjust subints and tvec
%                     subInts(i) = newSwInd;
                    tvec(i) = tvec(i-1) + tcross;
                    controlledActions = find(diff(subInts) ~=0);
                    nextControlledSwitching = controlledActions(find(find(diff(subInts) ~=0) > i,1,'first'));
                    if(isempty(nextControlledSwitching))
                        nextControlledSwitching = length(subInts);
                    end
                    subInts(i:nextControlledSwitching) = newSwInd;
                else
                    %% Error at end of period
                    crossXs = expm(M*tcross)*oldXs;
%                     violateMarginCross = Cbnd*crossXs(1:end-1) + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1);
                    [violateMarginCross,targetValCross] = obj.checkStateValidity(crossXs(1:end-1), u, subInts(i));
                    k=1;
%                     while ~all(abs(violateMarginCross(violateMarginEnd < 0)) < abs(obj.converter.topology.bndHyst((violateMarginEnd < 0),2)))
                    if(errAt == 1) 
                        %We know we have an error, so fixing it should be
                        %when at least one constraint is within the
                        %hysteresis bands of violation and all others are
                        %not in violation.
%                         while ~any(abs(targetValCross(violateMarginEnd < 0)) < obj.converter.topology.bndHyst((violateMarginEnd < 0),2)) || ... 
%                             ~all((violateMarginCross(violateMarginEnd < 0)) > 0) 

                        if any(violateMarginCross < 0 & violateMarginEnd < 0)
                            firstViolationCondition = violateMarginCross < 0 & violateMarginEnd < 0;
                            firstViolationTime = tcross;
                        else
                            firstViolationCondition = violateMarginEnd < 0;
                            firstViolationTime = deltaT;
                        end

                        while ~(...
                                any(abs(targetValCross(violateMarginEnd < 0)) < obj.converter.topology.bndHyst((violateMarginEnd < 0),2)) ...
                                && ...
                                all((violateMarginCross(violateMarginEnd < 0)) > 0) ...
                                )
                            tcrossOld = tcross;
                            tcross = tcross + 1./crossSlope.*(targetValCross);
                            
                            tcross = tcross(firstViolationCondition);

%                             if any(violateMarginCross < 0)
%                                 tcross = tcross(violateMarginCross < 0 & violateMarginEnd < 0);
%                             else
%                                 tcross = tcross(violateMarginEnd < 0);
%                             end
                            tcross(tcross > deltaT) = deltaT;
                            tcross(tcross<0) = 0;
                            tcross = min(tcross(tcross <= deltaT & tcross >= 0));

                            if tcross == deltaT && tcrossOld == deltaT
                                tcross = 0;
                            end

                            crossXs = expm(M*tcross)*oldXs;
                            violateMarginCrossOld = violateMarginCross;
                            [violateMarginCross,targetValCross] = obj.checkStateValidity(crossXs(1:end-1), u, subInts(i));

                            if tcross < firstViolationTime && any(violateMarginCross < 0 & violateMarginEnd < 0)
                                firstViolationCondition = violateMarginCross < 0 & violateMarginEnd < 0;
                                firstViolationTime = tcross;
                            end
                            if tcross ~= tcrossOld
                                crossSlope = -(violateMarginCross-violateMarginCrossOld)./(tcross-tcrossOld).*(violateMarginEnd < 0);
                            end
                            k=k+1;
                            if k>100
                                error('Exceeded max number of iterations in solving for zero crossing');
                                
                                if(0)
                                    tx = linspace(0,deltaT, 1000);
                                    for kk = 1:length(tx)
                                        crossXsx(:,kk) =  expm(M*tx(kk))*oldXs;
                                        [violateMarginCrossx(:,kk),~] = obj.checkStateValidity(crossXsx(1:end-1,kk), u, subInts(i));
                                    end
                                    plot(tx, violateMarginCrossx((violateMarginEnd < 0),:))
                                end
                                break
                            end
                        end
                    end

                    assert(all(violateMarginCross(violateMarginEnd < 0) > 0), ...
                        'Linear projection didnt find correct uncontrolled switching time' )

                    %% Adjust timings and add new interval
%                     FullPeriodViolateStart = zeros(size(violateMarginEnd,1),length(obj.ts));
%                     FullPeriodViolateEnd = zeros(size(violateMarginEnd,1),length(obj.ts));
                    currentInt = find(tvec(i) < cumsum(ts),1,'first');
%                     FullPeriodViolateEnd(:, currentInt) = violateMarginEnd;
%                     [tLocs,insertAt,adjType] = obj.findRequiredUncontrolledSwitching(FullPeriodViolateStart,FullPeriodViolateEnd);

                    %% add new interval if it doesn't exist
                    newSwVec = obj.converter.swvec(currentInt,:);
                    modSwitch = obj.topology.switchRef((violateMarginEnd < 0),1);
                    newState = ~obj.topology.switchRef((violateMarginEnd < 0),2);

                    newSwVec(modSwitch) = newState;
                    if ~ismember(newSwVec, obj.topology.swseq, 'rows')
                        obj.topology.loadCircuit([], [obj.converter.swvec; newSwVec], 0);
                    end
                    [~, newSwInd] = ismember(newSwVec, obj.topology.swseq, 'rows'); 

                    %% adjust subints and tvec 

%                     tvec(i+1) = tvec(i) + (diff(tvec(i:i+1))-tcross);
                    tvec(i+1) = tvec(i) + tcross;
                    controlledActions = find(diff(subInts) ~=0);
                    nextControlledSwitching = controlledActions(find(find(diff(subInts) ~=0) > i,1,'first'));
                    if(isempty(nextControlledSwitching))
                        nextControlledSwitching = length(subInts);
                    end
                    subInts(i+1:nextControlledSwitching) = newSwInd;
                end
                Xs = oldXs;
                j=j+1;
                if j > 1000
                    error(['exceeded maximum number of iterations at t = ' num2str(tvec(i))])
                end
            else
                i=i+1;
                j=1;
            end

        end
    catch e 
        disp(e.message);
        rethrow(e);
    end
    
    Xf = Xs;
    swinds = [subInts(diff(subInts)~=0) subInts(end)];
    us = cat(3,intus(:,:,diff(subInts)~=0), intus(:,:,end));
    ts = [tvec(diff(subInts)~=0); tvec(end)];
    ts = [ts(1); diff(ts)];
    
    zeroLocs = (ts == 0);
    ts(zeroLocs) = [];
    swinds(zeroLocs) = [];
    
    if debug == 1
        NS = size(states,1)-1;
        for i = 1:NS
            subplot(NS,1,i)
            plot(tvec(1:end-1),states(i,1:end-1),'LineWidth',3);
            ylabel(obj.stateNames{i})
        end
    end
    
    return;

% %     %% OLD Algorithm
% %     
% %     i=1;
% %     xs = zeros(size(obj.Xs,1),tVecLength);
% %     xs(:,1) = obj.Xs(:,1);
% %     tvec = zeros(tVecLength,1);
% %     tind = 2;
% %     
% %     try 
% %     while i < length(obj.ts)
% % %         [~, deltaT] = getDeltaT(obj.converter, i);
% %         deltaT = deltaTs(i);
% %         A = obj.converter.As(:,:,origSwind(i));
% %         B = obj.converter.Bs(:,:,origSwind(i));
% %         Cbnd = obj.converter.topology.Cbnd(:,:,origSwind(i));
% %         Dbnd = obj.converter.topology.Dbnd(:,:,origSwind(i));
% %         
% %         
% %         while max(tvec) < sum(obj.ts(1:i))
% %             if max(tvec) + deltaT > sum(obj.ts(1:i))
% %                 %Make sure last step lines up on an interface
% %                 deltaT = (max(tvec) + deltaT) - sum(obj.ts(1:i));
% %             end
% %             M = [A, B*obj.u];
% %             M = [M; zeros(diff(size(M)),size(M,2))];
% %             xtilde = [xs(:,tind-1); 1];
% %             
% %             nextX = expm(M*deltaT)*xtilde;
% %             nextX = nextX(1:end-1);
% %             
% % %             nextX =  expm(A*deltaT)*xs(:,tind-1) + A\(expm(A*deltaT)-eye(size(A)))*B*obj.u;
% %             violateMarginEnd = Cbnd*nextX + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1) + obj.converter.topology.bndHyst(:,2);
% %             
% %             if any(violateMarginEnd < 0)
% %                %go back and find error 
% %                xdotStart = A*xs(:,tind-1) + B*obj.u;
% %                vErrDotStart = Cbnd*xdotStart;
% %                targetValStart= Cbnd*xs(:,tind-1) + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1);
% %                
% %                xdotEnd = A*nextX + B*obj.u;
% %                vErrDotEnd = Cbnd*xdotEnd;
% %                targetValEnd = Cbnd*nextX + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1);
% %                
% %                projErrTime = [targetValStart./vErrDotStart, -targetValEnd./vErrDotEnd];
% %                newDeltaT = mean(projErrTime(violateMarginEnd<0,:),'all');
% %                
% %                warning('then edit the sequence');
% %                
% % %                conv.addUncontrolledSwitching(i,(-1)^(j+1),dt,switchRef(tLocs(:,i),1),~switchRef(tLocs(:,i),2));
% %                
% %                tvec(tind) = tvec(tind-1)+deltaT;
% %                tind = tind+1
% %             else 
% %                 xs(:,tind) = nextX;
% %                 tvec(tind) = tvec(tind-1)+deltaT;
% %                 tind = tind+1
% %             end
% %         end
% %         i = i+1
% %     end
% %     catch e 
% %         disp(e.message);
% %         asdfsa=1
% %     end
% %     
% %     Xf = xs;
% %     ts = tvec;
% %     swinds = 1;
% % 
% % end
% % 
