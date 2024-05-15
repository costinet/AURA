function plotModWf(varargin)
%plot modulation pattern
%
%   plotModWf(obj) obj is either an SMPSim object or SMPSconverter object
%
%   plotModWf(ts,swvec) ts and swvec are timing array and binary switching
%   state vector as specified by SMPSim objects
%
%   plotModWf(___,fn) plots to figure number specified by fn
%
%   plotModWf(___,fn, ls) applies linestyle (character array)
%   specified by ls
%
%   Example:
%       [ts1,swvec1] = dutyMod(.5,1e-6, 'dt', 100e-9, 'phase', -.25e-6, 'phaseUnits', 'time');
%       plotModWf(ts1,swvec1,1);
%
%   See also @SMPSim, dutyMod, phaseShiftMod, concatMods, combineMods,
%   combineModsN

    p = inputParser;
    p.addRequired('in1',@(x)(numel(x)>1 && isnumeric(x) && size(x,1) == 1) || (isa(x,'SMPSim')||isa(x,'SMPSconverter')));

    % Fake a swvec and ignore it to skip it in the input
    if isa(varargin{1},'SMPSim')||isa(varargin{1},'SMPSconverter') && numel(varargin) >= 2
        varargin = [varargin(1), {[0 0 0]}, varargin(2:end)];
    end

    p.addOptional('swvec',[],@(x)all(x==1 | x==0 ,'all'));
    
    p.addOptional('fn',get(gcf,'Number'),@isscalar)
    p.addOptional('linestyle','-',@ischar)

    parse(p,varargin{:});

    fn = p.Results.fn;
    linestyle = p.Results.linestyle;
    swvec = p.Results.swvec;

    in1 = p.Results.in1;

    if isa(in1,'SMPSim')||isa(in1,'SMPSconverter')
        sim = in1;
        ts = sim.ts;
        swvec = sim.swvec;
    else
        ts = in1;
        sim = [];
    end
    

    figure(fn);
    numSw = size(swvec,2);


    ted = repmat([0 cumsum(ts)],2,1);
    ted = ted(2:end-1);

    for i = 1:numSw
        subplot(numSw,1,i);
        
        state =  repmat(swvec(:,i),1,2)';
        state = state(:);

        plot(ted,state, linestyle,  'LineWidth', 3);

        if ~isempty(sim) &&  isa(sim,'SMPSim')
            ylabel(sim.switchNames{i});
        end
        if i< numSw
            set(gca, 'XTickLabel', [])
        end
    end

    linkaxes


end

