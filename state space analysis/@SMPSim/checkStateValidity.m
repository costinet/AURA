function [violateMargin,targetVal] = checkStateValidity(obj, X, u, swind)
%checkStateValidity Summary of this function goes here
%   violateMargin is how far from closes hystersis edge
%   targetVal is how far from nominal bound
   
%     swind = obj.converter.swind(swindi);
    Cbnd = obj.topology.Cbnd(:,:,swind);
    Dbnd = obj.topology.Dbnd(:,:,swind);

    violateMargin = Cbnd*X + Dbnd*u - obj.topology.bndHyst(:,1) + obj.topology.bndHyst(:,2);
    targetVal = Cbnd*X + Dbnd*u - obj.topology.bndHyst(:,1);

% warning('been messing here')
%     targetVal = - obj.topology.bndHyst(:,1);

end
