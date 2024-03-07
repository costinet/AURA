function [ts,swvec] = dutyMod(D,Ts, varargin)
%DUTYMOD generate duty-cycle specified modulation function
%
%   [ts,swvec] = dutyMod(D,Ts)
%
%   [ts,swvec] = dutyMod(D,Ts, 'paramName', paramVal)
%       paramNames include 'dt', 'phase', and 'phaseUnits'
%       phaseUnits may be 'rad', 'deg', or 'time'

    p = inputParser;
    p.addRequired('D',@(x)isnumeric(x)&&x<=1&&x>=0);
    p.addRequired('Ts',@(x)isnumeric(x));
    
    p.addParameter('dt',0)
    p.addParameter('phase',0)
    p.addParameter('phaseUnits','rad', ...
        @(x)strcmp(x,'rad') || strcmp(x,'deg') || strcmp(x,'time'));
    
    
    parse(p,D,Ts,varargin{:});

    D = p.Results.D;
    Ts = p.Results.Ts;
    phi = p.Results.phase;
    dt = p.Results.dt;

    if dt > 0
        swvec = [0 0
                 1 0 
                 0 0; 
                 0 1;
                 0 0];
        ts = [dt/2 (D*Ts-dt), dt,  ((1-D)*Ts-dt), dt/2];
    else
        swvec = [1 0; 
                 0 1];
        ts = [(D*Ts),  ((1-D)*Ts)];
    end

    if phi ~= 0
        [ts,swvec] = phaseShiftMod(ts,swvec, phi, 'phaseUnits', p.Results.phaseUnits);
%         if strcmp(p.Results.phaseUnits, 'rad')
%             tphi = mod(phi/(2*pi)*Ts,Ts);
%         elseif strcmp(p.Results.phaseUnits, 'deg')
%             tphi = mod(phi/(360)*Ts,Ts);
%         elseif strcmp(p.Results.phaseUnits, 'time')
%             tphi = mod(phi,Ts);
%         end
% 
%         
%         cumts = [cumsum(ts)] + tphi;
%         rshift = 0;
% 
%         if tphi > 0
%             brk = find(cumts > Ts, 1);
%             rshift = sum(cumts > Ts);
% 
%             if ~isempty(brk)
%                 cumts = [cumts(1:brk-1), Ts, cumts(brk:end) - Ts];
%                 swvec = [swvec(1:brk-1,:); swvec(brk,:); swvec(brk:end,:)];
% 
%                 swvec = circshift(swvec,rshift,1);    
%                 ts = diff([0 circshift(cumts,rshift)]); 
%             end
% %         elseif tphi < 0 %Shouldn't be reached because of mod?
% %             brk = find(cumts > 0, 1);
% %             rshift =  -sum(cumts < 0);
% % 
% %             if ~isempty(brk)
% %                 cumts = [cumts(1:brk-1) + Ts, 0, cumts(brk:end)];
% %                 swvec = [swvec(:,1:brk-1), swvec(:,brk-1), swvec(:,brk:end)];
% % 
% %                 swvec = circshift(swvec,rshift,2);    
% %                 ts = diff([0 circshift(cumts,rshift)]); 
% %             end
%         end
    end
end