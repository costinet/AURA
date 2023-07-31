function [violateMargin,targetVal] = checkStateValidity(obj, X, u, swind)
%checkStateValidity use constraints to calculate validity of state vector
%   
%   [violateMargin,targetVal] = checkStateValidity(obj, X, u, swind)
%   obj is a SMPSim object, X and u are the state and input vectors, and
%   swind is an index into the SMPStopology class to get the desired
%   switching interval constraint matrices.
%       violateMargin is how far from closest hystersis edge
%       targetVal is how far from nominal value
%
%   See Also SMPSim, SMPStopology, circuitParser
   
%     swind = obj.converter.swind(swindi);
    Cbnd = obj.topology.Cbnd(:,:,swind);
    Dbnd = obj.topology.Dbnd(:,:,swind);

    violateMargin = Cbnd*X + Dbnd*u - obj.topology.bndHyst(:,1) + obj.topology.bndHyst(:,2);
    targetVal = Cbnd*X + Dbnd*u - obj.topology.bndHyst(:,1);

% warning('been messing here')
%     targetVal = - obj.topology.bndHyst(:,1);

end
