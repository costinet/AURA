function [ x, fval, history ] = fminconHistory(fun, x0, xmin, xmax)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    history = [];
    options = optimset('Display','iter-detailed');
%     options.PlotFcns = @optimplotx;
    options.OutputFcn =  @outfun;
    options.Algorithm =  'interior-point';
    A = []; b=[]; Aeq=[]; beq=[];
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq, xmin, xmax, [], options);
    
     function stop = outfun(x,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
           history = [history; x, optimvalues.fval];
%            disp([x, optimvalues.fval]);
        end
    end


end

