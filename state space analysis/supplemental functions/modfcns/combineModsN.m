function [tsO,swvecO] = combineModsN(ts,swvec,switchCol)
%COMBINEMODS combine two modulation functions
%   [ts,swvec] = combineMods(ts1,swvec1,switchCol1,ts2,swvec2,switchCol2, ...)
%   ts is a vector of time durations, swvec is a binary matrix of switching
%   states and switchCol is the indexes of switches in the combined
%   modulation
%
%   if switchCol is 2-D, each row is a copy of the associated swvec. 

%     p = inputParser;
%     p.addRequired('ts1',@(x)numel(x)>1 && isnumeric(x) && size(x,1) == 1);
%     p.addRequired('swvec1',@(x)all(x==1 | x==0 ,'all'));
% %     p.addOptional('switchCol1',1:size(swvec1,2),@(x)isnumeric(x) && length(x) == size(swvec1,1));
%     p.addRequired('switchCol1',@(x)isnumeric(x) && size(x,2) == size(swvec1,2));
% 
%     p.addRequired('ts2',@(x)numel(x)>1 && isnumeric(x) && size(x,1) == 1);
%     p.addRequired('swvec2',@(x)all(x==1 | x==0 ,'all'));
% %     p.addOptional('switchCol2',1:size(swvec2,2),@(x)isnumeric(x) && length(x) == size(swvec2,1));
%     p.addRequired('switchCol2',@(x)isnumeric(x) && size(x,2) == size(swvec2,2));
%     
% %     p.addOptional('fn',get(gcf,'Number'))
% 
%     p.parse(ts1,swvec1,switchCol1,ts2,swvec2,switchCol2, varargin{:});
% 
%     ts1 = p.Results.ts1;
%     ts2 = p.Results.ts2;

    arguments (Repeating)
        ts (1,:) {mustBeNumeric,mustBeReal}
        swvec {mustBeNumeric,mustBeReal,mustMatchts(ts,swvec)}
        switchCol {mustBeNumeric,mustBeReal,mustMatchSwvec(swvec,switchCol)}
    end

    numMods = numel(swvec);

    sumts = zeros(numMods,1);
    for j = 1:numMods
        sumts(j) = sum(ts{j});
        cumts{j} = cumsum(ts{j});
    end

    assert(all(sumts - sumts(1) < 10*max(eps(sumts))),'to combine modulations, they must have the same overall period');

    if numMods == 1
        tsO = ts;
        swvecO = swvec;
        return
    end

    [tsO,swvecO] = combineMods(ts{1},swvec{1},switchCol{1}, ts{2},swvec{2},switchCol{2});

    for j = 2:numMods-1
        [tsO,swvecO] = combineMods(tsO,swvecO,1:size(swvecO,2), ts{j+1},swvec{j+1},switchCol{j+1});
    end


% 
%     for j = 1:numMods
%         
%     
%         [cumts,~] = sort([cumts1,cumts2]);
%     
%         ts1loc = sum((cumts1' < cumts-eps),1);
%         ts2loc = sum((cumts2' < cumts-eps),1);
%     
%         switchCol1 = switchCol1';
%         switchCol2 = switchCol2';
%     
%         for i = 1:length(cumts)
%             swvec(i,switchCol1(:)) = repmat(swvec1(ts1loc(i)+1,:),1,size(switchCol1,2));
%             swvec(i,switchCol2(:)) = repmat(swvec2(ts2loc(i)+1,:),1,size(switchCol2,2));
%         end
%     
%         ts = diff([0 cumts]);
%     
%         zeroints = ts <= 10*max(eps(ts));
%     
%         ts(zeroints) = [];
%         swvec(zeroints,:) = [];
%     end

end


function mustMatchts(ts,swvec)
    if ~isequal(size(swvec,1),numel(ts))
        eid = 'Size:notEqual';
        msg = 'The number of rows of swvec must match the number of columns of ts';
        error(eid,msg)
    end
end

function mustMatchSwvec(swvec,switchCol)
    if ~isequal(size(swvec,2),size(switchCol,2))
        eid = 'Size:notEqual';
        msg = 'The number of columns of switchCol must match the number of columns of swvec';
        error(eid,msg)
    end
end