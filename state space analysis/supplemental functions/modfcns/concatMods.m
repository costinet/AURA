function [ts,swvec] = concatMods(ts1,swvec1,ts2,swvec2)
%CONCATMODS concatenate two modulation patterns sequentially
%
%   [ts,swvec] = concatMods(ts1,swvec1,ts2,swvec2)
     p = inputParser;
    p.addRequired('ts1',@(x)numel(x)>1 && isnumeric(x) && size(x,1) == 1);
    p.addRequired('swvec1',@(x)all(x==1 | x==0 ,'all'));

    p.addRequired('ts2',@(x)numel(x)>1 && isnumeric(x) && size(x,1) == 1);
    p.addRequired('swvec2',@(x)all(x==1 | x==0 ,'all'));

    p.parse(ts1,swvec1,ts2,swvec2);

    assert(size(swvec1,2)==size(swvec2,2),'Both swvec must have the same number of columns, corresponding to the same number of switches.')
    ts = [ts1, ts2];
    swvec = [swvec1; swvec2];
end