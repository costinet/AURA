function [Xs] = perturbedSteadyState(obj, dts)
% Function set up to accomodate DC code


    As = obj.As;
    Bs = obj.Bs;
    ts = obj.ts + dts;
    u = obj.u;
    
[Xs] = obj.SS_Soln_Aug(0,As,Bs,ts,u);

end

