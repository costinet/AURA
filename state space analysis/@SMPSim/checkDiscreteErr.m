function [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = checkDiscreteErr(obj)
%checkDiscreteErr find violations of constraint matrices at each discrete
%switching interface
%   
%   [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = checkDiscreteErr(obj)
%   using input SMPSim and its current steady-state solution, find the
%   errors according to constraint matrices at each switching interface.
%   Outputs *.Start refer to values at the beginning of any switching
%   interval and *.End refer to values at the end.  Output violateMargin is the
%   distance between current value of the signal and the (value +
%   tolerance) specified by the constraint.  Output targetVal is the same
%   without tolerance included.
%
%   checkDiscreteErr calls checkStateValidity recursively, and makes a few
%   corrections for common false positives.
%
%   See Also SMPSim.checkStateValidity, SMPStopology.constraints,
%   circuitParser

    Xss = obj.Xs;
%     us = obj.u;
    Cbnd = obj.topology.Cbnd(:,:,obj.converter.swind);
    Dbnd = obj.topology.Dbnd(:,:,obj.converter.swind);
    
    violateMarginStart = zeros(size(obj.topology.Cbnd,1) ,length(obj.converter.swind));
    targetValStart = violateMarginStart;
    violateMarginEnd = zeros(size(obj.topology.Cbnd,1) ,length(obj.converter.swind));
    targetValEnd = violateMarginEnd;

    for i = 1:length(obj.converter.swind)
        us = obj.fullu(:,:,i);
        [violateMarginStart(:,i), targetValStart(:,i)] = obj.checkStateValidity( Xss(:,i), us, obj.converter.swind(i));
        [violateMarginEnd(:,i), targetValEnd(:,i)] = obj.checkStateValidity( Xss(:,i+1), us, obj.converter.swind(i));

        
        %% Check for and correct hard switching issue:
        if any(violateMarginStart(:,i) < 0 & obj.topology.constraints.switchRef(:,2) == 0 )
            % error at start             &         its a device that is off          
            prevInt = subsref(circshift(obj.converter.swind,1),struct('type','()','subs',{{i}}));
            prevSwStates = obj.topology.swseq(prevInt, obj.topology.constraints.switchRef(:,1));
            if any(violateMarginStart(:,i) < 0 & obj.topology.constraints.switchRef(:,2) == 0 & prevSwStates')
                % error at start             &         its a device that is off       &       It was on
                falsePos = violateMarginStart(:,i) < 0 & obj.topology.constraints.switchRef(:,2) == 0 & prevSwStates';
                swind = obj.converter.swind(i);
                xdot = obj.As(:,:,i)*Xss(:,i) + obj.Bs(:,:,i)*us;
                deltaViolate = Cbnd(:,:,i)*xdot + Dbnd(:,:,i)*us;

                trueFalsePos = deltaViolate > 0 & falsePos;
                if(any(trueFalsePos))
                    violateMarginStart(trueFalsePos,i) = 0;
%                     warning('This is under test -- just deltaViolate > 0 is maybe too generous')
                end
                %the point of this whole thing is a diode can be on during
                %DT, turn off because the opposite FET came on, then try to
                %turn back on because INSTATNTANEOUSLY at the start of the
                %inverval, the I*R from when it was on makes it look like
                %if has a positive voltage.
                % One could conceivably address this with the hystersis on
                % the constraints, but it requires alteration based on what
                % the actual I*R is.
                % It cannot be neglected because maybe the FET
                % turning on or whatever made the switching transition
                % won't be taking over this current.
            end

        end

        %% More generally, check for any negligible error
        %If there is an error at the start of the interval, and it'll
        %resolve in less than 1.5 minimum time steps
        if any(violateMarginStart(:,i) < 0)
            swind = obj.converter.swind(i);
            xdot = obj.As(:,:,i)*Xss(:,i) + obj.Bs(:,:,i)*us;
            deltaViolate = Cbnd(:,:,i)*xdot + Dbnd(:,:,i)*us;

            negligibleError = violateMarginStart(:,i) < 0 & ...
                (violateMarginStart(:,i) + deltaViolate*obj.converter.timingThreshold*1.05) > 0;
            if(any(negligibleError))
                violateMarginStart(negligibleError,i) = 0;
                warning('check this, it is untested')
            end

        end

        %Same at the end of the interval
        if any(violateMarginEnd(:,i) < 0)
            swind = obj.converter.swind(i);
            nextXs = circshift(Xss,-1,2);
            xdot = obj.As(:,:,i)*nextXs(:,i) + obj.Bs(:,:,i)*us;
            deltaViolate = Cbnd(:,:,i)*xdot + Dbnd(:,:,i)*us;

            negligibleError = violateMarginEnd(:,i) < 0 & ...
                (violateMarginEnd(:,i) - deltaViolate*obj.converter.timingThreshold*1.05) > 0;

            if(any(negligibleError))
                violateMarginEnd(negligibleError,i) = 0;
                warning('check this, it is untested')
            end

        end

%         %% what if there is not an error, but there is about to be one
%         swind = obj.converter.swind(i);
%         xdot = obj.As(:,:,i)*Xss(:,i) + obj.Bs(:,:,i)*us;
%         deltaViolate = Cbnd(:,:,i)*xdot + Dbnd(:,:,i)*us;
% 
%         actuallyAnError = violateMarginStart(:,i) > 0 & ...
%             (violateMarginStart(:,i) + deltaViolate*obj.converter.ts(i)/250) < 0;
%         if(any(actuallyAnError))
%             violateMarginStart(actuallyAnError,i) = -1;
%             warning('check this, it is untested')
%         end


        
%         violateMarginStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1) + obj.topology.bndHyst(:,2);
%         targetValStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1);
%         violateMarginEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1) + obj.topology.bndHyst(:,2);
%         targetValEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1);
    end

end

