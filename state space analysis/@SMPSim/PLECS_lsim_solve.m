function [Xss] = PLECS_lsim_solve(obj,IC,assign)
% This function takes the current steady state solution to the
% converter and updates it based on lsim-ing through the converter for
% 1 period. In this simulation there are no presumed diode conduction
% states. Diode conductions are assessed after every lsim to find if
% there was a violation. If there was deemed to be a violation then
% the converter state matricies are changed at the point of violation
% and then lsimed for the rest of the period. The steady state of the
% converter is chosen to be the final value of the state space
% solution.

try

if nargin==2
    assign = false;
end

As = obj.As;
Bs = obj.Bs;
Cs = obj.Cs;
Ds = obj.Ds;
u = obj.u;
eigA = obj.eigA;
ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
ts = obj.ts;
Xs = obj.Xs;
i = 1;
tol = 0.05;

%{

    for i = 1:nsub
        Aas(:,:,i) = [As(:,:,i), Bs(:,:,i)*u; 
                        zeros(1,ns+1)];
    end


    tsteps = 1000;
    t = linspace(0,Ts, tsteps);
    dt = Ts/tsteps;
    
    xplot = zeros(ns+1,tsteps);
    xplot(:,1) = X0;
    for i = 1:tsteps-1
        j = find(t(i)<=cumsum(ti),1,'first');
        xplot(:,i+1) = expm(Aas(:,:,j)*(dt))*xplot(:,i);
    end
%}



while i <= size(ONorOFF,2)
    
    Samples = min(pi./abs(imag(eigA(:,i))))/4\ts(i);
    
    % Run lsim for given time period
    C = zeros(1,size(As,1), 1);
    SS = ss(As(:,:,i), Bs(:,:,i), C, 0);
    tsim = linspace(0, ts(i), Samples);
    [~, ~, x] = lsim(SS, u*ones(size(tsim)), tsim, IC);
    
    % Find statemetn to determine if and where there are violation in
    % the lsim
    Infx = x;
    negInfx = x;
    Infx(:,[(ONorOFF(:,i)~=-1)']) = Inf;
    negInfx(:,[(ONorOFF(:,i)~=1)']) = -Inf;
    [Xoff,Yoff] = find(Infx<-1-tol);
    
    
    [Xon,Yon] = find(negInfx>-1+tol);
    
    % If there is no violation then continue on to the next interation
    if isempty(Xoff) && isempty(Xon)
        i = i+1;
        IC = x(end,:)';
        continue
    end
    
    if ~isempty(min(Xon)<min(Xoff)) || isempty(Xoff)
        min_Xon = min(Xon);
        
        
        if min_Xon~=1
            
            min_index=find(Xon==min_Xon);
            
            new_state = ONorOFF(:,i); % Find state that needs to be copied
            new_state(Yon(min_index)) = -1; % make correction in state
            ONorOFF = [ONorOFF(:,1:i), new_state ,ONorOFF(:,i+1:end) ]; % Place new state space in the ONorOFF matrix
            [A,B,C,D,eigA_n] = obj.add_state_matrix(new_state); % Calculate the needed
            
            ts = [ts(1:i-1) tsim(min_Xon) ts(i)-tsim(min_Xon) ts(i+1:end)];
            
            As(:,:,1:i) = As(:,:,1:i);
            As(:,:,i+2:size(As,3)+1) = As(:,:,i+1:end);
            As(:,:,i+1) =  A;
            Bs(:,:,1:i) = Bs(:,:,1:i);
            Bs(:,:,i+2:size(Bs,3)+1) = Bs(:,:,i+1:end) ;
            Bs(:,:,i+1) =  B;
            Cs(:,:,1:i) = Cs(:,:,1:i);
            Cs(:,:,i+2:size(Cs,3)+1) = Cs(:,:,i+1:end) ;
            Cs(:,:,i+1) =  C;
            Ds(:,:,1:i) = Ds(:,:,1:i);
            Ds(:,:,i+2:size(Ds,3)+1) = Ds(:,:,i+1:end) ;
            Ds(:,:,i+1) =  D;
            eigA(:,1:i) = eigA(:,1:i);
            eigA(:,i+2:size(eigA,2)+1) = eigA(:,i+1:end) ;
            eigA(:,i+1) = eigA_n;
            
            i = i+1;
            IC = x(min_Xon,:)';
            continue
        end
        
        
        
    end
    
    if ~isempty(min(Xoff)<=min(Xon)) || isempty(Xon)
        min_Xoff = min(Xoff);
        
        if min_Xoff~=1
            
            min_index=find(Xoff==min_Xoff);
            
            new_state = ONorOFF(:,i); % Find state that needs to be copied
            new_state(Yoff(min_index)) = 1; % make correction in state
            ONorOFF = [ONorOFF(:,1:i), new_state ,ONorOFF(:,i+1:end) ]; % Place new state space in the ONorOFF matrix
            [A,B,C,D,eigA_n] = obj.add_state_matrix(new_state); % Calculate the needed
            
            ts = [ts(1:i-1) tsim(min_Xoff) ts(i)-tsim(min_Xoff) ts(i+1:end)];
            
            As(:,:,1:i) = As(:,:,1:i);
            As(:,:,i+2:size(As,3)+1) = As(:,:,i+1:end);
            As(:,:,i+1) =  A;
            Bs(:,:,1:i) = Bs(:,:,1:i);
            Bs(:,:,i+2:size(Bs,3)+1) = Bs(:,:,i+1:end) ;
            Bs(:,:,i+1) =  B;
            Cs(:,:,1:i) = Cs(:,:,1:i);
            Cs(:,:,i+2:size(Cs,3)+1) = Cs(:,:,i+1:end) ;
            Cs(:,:,i+1) =  C;
            Ds(:,:,1:i) = Ds(:,:,1:i);
            Ds(:,:,i+2:size(Ds,3)+1) = Ds(:,:,i+1:end) ;
            Ds(:,:,i+1) =  D;
            eigA(:,1:i) = eigA(:,1:i);
            eigA(:,i+2:size(eigA,2)+1) = eigA(:,i+1:end) ;
            eigA(:,i+1) = eigA_n;
            
            i = i+1;
            IC = x(min_Xoff,:)';
            continue
        end
        
    end
    
    
    
end
Xss = x(end,:)';


if assign
    obj.As = As;
    obj.Bs = Bs;
    obj.Cs = Cs;
    obj.Ds = Ds;
    obj.ts = ts;
    obj.eigA = eigA;
    obj.Converter.Topology.Parser.ONorOFF = ONorOFF;
end

catch ME
    rethrow(ME)
end

end
