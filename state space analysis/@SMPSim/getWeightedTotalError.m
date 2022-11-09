function [weightTotalErr] = getWeightedTotalError(obj, errBefore,errAfter)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    err= cat(3,errBefore,errAfter);
    for z = 1:2
        errLocs = find(err(:,:,z) < 0);
        [r,c] = ind2sub(size(err(:,:,z)),errLocs);
    
        currentVioRows = find(obj.converter.topology.constraints.switchRef(:,2) == 1);
        [~,currentErrLocs] = intersect(r, currentVioRows);
        sw = r(currentErrLocs);
        col = c(currentErrLocs);
        ti = obj.converter.swseq(col);
        
        rDs = zeros(length(sw),1);
        for j = 1:length(sw)
            rDs(j) = -1/min(obj.converter.topology.constraints.Dbnd(sw(j),:,ti(j)));
            err(sw(j),col(j),z) = err(sw(j),col(j),z).*rDs(j);
        end
    end
    weightTotalErr = sum(err,'all');
    
%     err(sim.converter.topology.constraints.switchRef(:,2) == 1,:) = ...
%         err(sim.converter.topology.constraints.switchRef(:,2) == 1,:).*rDs
end