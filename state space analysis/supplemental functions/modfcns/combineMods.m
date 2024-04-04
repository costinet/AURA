function [ts,swvec] = combineMods(ts1,swvec1,switchCol1,ts2,swvec2,switchCol2, varargin)
%COMBINEMODS combine two modulation functions
%   [ts,swvec] = combineMods(ts1,swvec1,switchCol1,ts2,swvec2,switchCol2)
%   ts is a vector of time durations, swvec is a binary matrix of switching
%   states and switchCol is the indexes of switches in the combined
%   modulation
%
%   if switchCol is 2-D, each row is a copy of the associated swvec.    


    p = inputParser;
    p.addRequired('ts1',@(x)numel(x)>1 && isnumeric(x) && size(x,1) == 1);
    p.addRequired('swvec1',@(x)all(x==1 | x==0 ,'all'));
%     p.addOptional('switchCol1',1:size(swvec1,2),@(x)isnumeric(x) && length(x) == size(swvec1,1));
    p.addRequired('switchCol1',@(x)isnumeric(x) && size(x,2) == size(swvec1,2));

    p.addRequired('ts2',@(x)numel(x)>1 && isnumeric(x) && size(x,1) == 1);
    p.addRequired('swvec2',@(x)all(x==1 | x==0 ,'all'));
%     p.addOptional('switchCol2',1:size(swvec2,2),@(x)isnumeric(x) && length(x) == size(swvec2,1));
    p.addRequired('switchCol2',@(x)isnumeric(x) && size(x,2) == size(swvec2,2));
    
%     p.addOptional('fn',get(gcf,'Number'))

    p.parse(ts1,swvec1,switchCol1,ts2,swvec2,switchCol2, varargin{:});

    ts1 = p.Results.ts1;
    ts2 = p.Results.ts2;

    assert(abs(sum(ts1)-sum(ts2))<10*max(eps(ts1)),'to combine modulations, they must have the same overall period');

    cumts1 = cumsum (ts1);
    cumts2 = cumsum(ts2);

    [cumts,~] = sort([cumts1,cumts2]);

    ts1loc = sum((cumts1' < cumts-eps),1);
    ts2loc = sum((cumts2' < cumts-eps),1);

    switchCol1 = switchCol1';
    switchCol2 = switchCol2';

    for i = 1:length(cumts)
        swvec(i,switchCol1(:)) = repmat(swvec1(ts1loc(i)+1,:),1,size(switchCol1,2));
        swvec(i,switchCol2(:)) = repmat(swvec2(ts2loc(i)+1,:),1,size(switchCol2,2));
    end

    ts = diff([0 cumts]);

    zeroints = ts <= 10*max(eps(ts));

    ts(zeroints) = [];
    swvec(zeroints,:) = [];

end