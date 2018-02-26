function [ x, fval, history ] = fminHistory(fun, x0)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    history = [];
    options = optimset('Display','iter');
%     options.PlotFcns = @optimplotx;
    options.OutputFcn =  @outfun;
    [x,fval,exitflag,output] = fminsearch(fun,x0, options);
    
     function stop = outfun(x,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
           history = [history; x, optimvalues.fval];
%            disp([x, optimvalues.fval]);
        end
    end


end

