function plotModWf(ts,swvec,varargin)
%PLOTMODWF plot modulation pattern
%   plotModWf(ts,swvec)
%   plotModWf(ts,swvec,fn) plots to figure number specified by fn

    p = inputParser;
    p.addRequired('ts',@(x)numel(x)>1 && isnumeric(x) && size(x,1) == 1);
    p.addRequired('swvec',@(x)all(x==1 | x==0 ,'all'));
    
    p.addOptional('fn',get(gcf,'Number'))

    parse(p,ts,swvec,varargin{:});

    fn = p.Results.fn;

    figure(fn);
    numSw = size(swvec,2);


    ted = repmat([0 cumsum(ts)],2,1);
    ted = ted(2:end-1);

    for i = 1:numSw
        subplot(numSw,1,i);
        
        state =  repmat(swvec(:,i),1,2)';
        state = state(:);

        plot(ted,state, 'LineWidth', 3);
    end

    linkaxes


end

