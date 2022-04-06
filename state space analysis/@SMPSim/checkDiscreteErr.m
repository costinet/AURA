function [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = checkDiscreteErr(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Xss = obj.Xs;
    us = obj.u;
    Cbnd = obj.topology.Cbnd(:,:,obj.converter.swind);
    Dbnd = obj.topology.Dbnd(:,:,obj.converter.swind);
    
    violateMarginStart = zeros(size(obj.topology.Cbnd,1) ,length(obj.converter.swind));
    targetValStart = violateMarginStart;
    violateMarginEnd = zeros(size(obj.topology.Cbnd,1) ,length(obj.converter.swind));
    targetValEnd = violateMarginEnd;

    for i = 1:length(obj.converter.swind)
        violateMarginStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1) + obj.topology.bndHyst(:,2);
        targetValStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1);
        violateMarginEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1) + obj.topology.bndHyst(:,2);
        targetValEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us - obj.topology.bndHyst(:,1);
    end

end

