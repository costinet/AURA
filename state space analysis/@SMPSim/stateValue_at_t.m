function [ x, xdot ] = stateValue_at_t(obj, x0, t, si)
%for switching interval si, calculate the state values at t, starting from
%initial states x0

    As = obj.converter.topology.As;
    Bs = obj.converter.topology.Bs;
    u = obj.converter.u;

    if(numel(t) == 1)
        [fresp, ~] = obj.forcedResponse(As(:,:,si), expm(As(:,:,si)*t), Bs(:,:,si), u, t, 1);
        x = expm(As(:,:,si)*t)*x0 + fresp;
        xdot = As(:,:,si)*x + Bs(:,:,si)*u;
    else
        % Probably should thrown an error (lsim is likely better) but
        % incorporated here for debugging.
        x = zeros(length(x0),length(t));
        for i = 1:length(t)
            [fresp, ~] = obj.forcedResponse(As(:,:,si), expm(As(:,:,si)*t(i)), Bs(:,:,si), u, t(i), 1);
            x(:,i) = expm(As(:,:,si)*t(i))*x0 + fresp;
            xdot(:,i) = As(:,:,si)*x(:,i) + Bs(:,:,si)*u;
        end
    end
end

