function [violateMargin,targetVal] = checkStateValidity(obj, X, u, swind)
%checkStateValidity Summary of this function goes here
%   violateMargin is how far from closes hystersis edge
%   targetVal is how far from nominal bound
   
    Cbnd = obj.topology.Cbnd(:,:,swind);
    Dbnd = obj.topology.Dbnd(:,:,swind);

% %     if swind == 5
% %         x=1;
% % %         obj.converter.topology.Cs(:,:,5)
% % %         obj.converter.topology.Ds(:,:,5)
% %         obj.converter.swseq
% %         obj.converter.topology.Cbnd(:,:,5)
% %         obj.converter.topology.Dbnd(:,:,5)
% % warning('test2')
% %     end
    
%     violateMargin = zeros(size(Cbnd,1) ,length(swind));
%     targetVal = violateMargin;


    violateMargin = Cbnd*X + Dbnd*u - obj.topology.bndHyst(:,1) + obj.topology.bndHyst(:,2);
    targetVal = Cbnd*X + Dbnd*u - obj.topology.bndHyst(:,1);

% warning('been messing here')
%     targetVal = - obj.topology.bndHyst(:,1);

end
