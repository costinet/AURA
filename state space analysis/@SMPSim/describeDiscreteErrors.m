function describeDiscreteErrors(obj)
%describeDiscreteErrors prints to command window a text description of the
%discrete time errors found at the current steady-state candidate solution
%   
%   For debugging purposes only.
%   
%   See Also SMPSim.describeAlteredTimes,
%   SMPSim.describeInsertedIntervals, SMPSim.describeSwitchState

    [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = obj.checkDiscreteErr();
    errBefore = min(violateMarginStart,0);
    errAfter = min(violateMarginEnd,0);

    for ti = 1:size(errBefore,2)
        vioLocs = find(errBefore(:,ti) < 0);
        fulli = find(obj.converter.fullts(:) > 0, ti);
        [tii, tij] = ind2sub(size(obj.converter.fullts),fulli(end));
        for si = vioLocs'
            if ~isempty(si)
                switchRef = obj.converter.topology.constraints.switchRef(si,1);
                switchName = obj.switchNames(switchRef);
                violationType = obj.converter.topology.constraints.switchRef(si,2);
                
    
                disp(['at beginning of interval ' num2str(ti) ' (' num2str(tij) ',' num2str(tii) ') switch ' switchName{:} ' is in violation']);
                if violationType == 0
                    disp([' --- it is off, but its voltage is ' num2str(targetValStart(si,ti) + obj.converter.topology.constraints.bndHys(si,1)) ' > ' num2str(obj.converter.topology.constraints.bndHys(si,1)) ' +/-' num2str(obj.converter.topology.constraints.bndHys(si,2))]);
                else
                    disp([' --- it is on, but its current is ' num2str(targetValStart(si,ti) + obj.converter.topology.constraints.bndHys(si,1)) ' > ' num2str(obj.converter.topology.constraints.bndHys(si,1)) ' +/-' num2str(obj.converter.topology.constraints.bndHys(si,2))]);
                end
            end
        end

        vioLocs = find(errAfter(:,ti) < 0);
        for si = vioLocs'
            if ~isempty(si)
                switchRef = obj.converter.topology.constraints.switchRef(si,1);
                switchName = obj.switchNames(switchRef);
                violationType = obj.converter.topology.constraints.switchRef(si,2);
                
    
                disp(['at end of interval ' num2str(ti) ' (' num2str(tij) ',' num2str(tii) ') switch ' switchName{:} ' is in violation']);
                if violationType == 0
                    disp([' --- it is off, but its voltage is ' num2str(targetValEnd(si,ti) + obj.converter.topology.constraints.bndHys(si,1)) ' > ' num2str(obj.converter.topology.constraints.bndHys(si,1)) ' +/-' num2str(obj.converter.topology.constraints.bndHys(si,2))]);
                else
                    disp([' --- it is on, but its current is ' num2str(targetValEnd(si,ti) + obj.converter.topology.constraints.bndHys(si,1)) ' > ' num2str(obj.converter.topology.constraints.bndHys(si,1))  ' +/-' num2str(obj.converter.topology.constraints.bndHys(si,2))]);
                end
            end
        end

    end

end