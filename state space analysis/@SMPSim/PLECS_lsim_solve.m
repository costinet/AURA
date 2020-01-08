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
    
    As = obj.As_OG;
    Bs = obj.Bs_OG;
    Cs = obj.Cs_OG;
    Ds = obj.Ds_OG;
    u = obj.u_OG;
    eigA = obj.eigA_OG;
    ONorOFF = obj.ONorOFF_OG;
    ts = obj.ts_OG;
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
    for i = 1:size(As,3)
        Aas(:,:,i) = [As(:,:,i), Bs(:,:,i)*u;
            zeros(1,size(As,1)+1)];
    end
    
    i = 1;
    x(1,:) = [IC',1];
    while i <= size(ONorOFF,2)
        
        Samples = min(pi./abs(imag(eigA(:,i))))/4\ts(i);
        
        Samples = max(Samples,10);
        
        % Run lsim for given time period
        dt = ts(i)/Samples;
        
        noviolation = true;
        notend = true;
        count = 0;
        while noviolation && notend 
            x(end+1,:) = (expm(Aas(:,:,i)*(dt))*x(end,:)')';
            noviolation =  ~(sum(sum(x(:,[(ONorOFF(:,i)==-1)'])<-1-tol)) | sum(sum(x(:,[(ONorOFF(:,i)==1)'])>-1+tol)));
            count = count+1;
            if count == floor(Samples)
                notend = false;
            end
        end
        
        % Need to check each state as I go to determine if there is a
        % violation
        
        
        Infx = x(:,1:end-1);
        negInfx = x(:,1:end-1);
        Infx(:,[(ONorOFF(:,i)~=-1)']) = Inf;
        negInfx(:,[(ONorOFF(:,i)~=1)']) = -Inf;
        [Xoff,Yoff] = find(Infx<-1-tol);
        [Xon,Yon] = find(negInfx>-1+tol);
        
        
        
        %{
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
        %}
        % If there is no violation then continue on to the next interation
        if isempty(Xoff) && isempty(Xon)
            i = i+1;
            IC = x(end,:);
            clear x
            x = IC;
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
                
                ts = [ts(1:i-1) count*dt ts(i)-count*dt ts(i+1:end)];
                
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
                Aas(:,:,1:i) = Aas(:,:,1:i);
                Aas(:,:,i+2:size(Aas,3)+1) = Aas(:,:,i+1:end);
                Aas(:,:,i+1) =  [A, B*u;
                    zeros(1,size(A,1)+1)];
                
                
                
                
                
                i = i+1;
                IC = x(end,:);
                clear x
                x = IC;
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
                
                ts = [ts(1:i-1) count*dt ts(i)-count*dt ts(i+1:end)];
                
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
                Aas(:,:,1:i) = Aas(:,:,1:i);
                Aas(:,:,i+2:size(Aas,3)+1) = Aas(:,:,i+1:end);
                Aas(:,:,i+1) =  [A, B*u;
                                 zeros(1,size(A,1)+1)];
                
                
                i = i+1;
                IC = x(end,:);
                clear x
                x = IC;
                continue
                
            
            else % If the inital condition needed to have a diode on
                
                min_index=find(Xoff==min_Xoff);
                
                new_state = ONorOFF(:,i); % Find state that needs to be copied
                new_state(Yoff(min_index)) = 1; % make correction in state
                ONorOFF = [new_state ,ONorOFF(:,1:end) ]; % Place new state space in the ONorOFF matrix
                [A,B,C,D,eigA_n] = obj.add_state_matrix(new_state); % Calculate the needed
                
                ts = [count*dt ts(i)-count*dt ts(i+1:end)];
                
              
                As(:,:,i+1:size(As,3)+1) = As(:,:,1:end);
                As(:,:,i) =  A;
              
                Bs(:,:,i+1:size(Bs,3)+1) = Bs(:,:,1:end) ;
                Bs(:,:,i) =  B;
              
                Cs(:,:,i+1:size(Cs,3)+1) = Cs(:,:,1:end) ;
                Cs(:,:,i) =  C;
              
                Ds(:,:,i+1:size(Ds,3)+1) = Ds(:,:,1:end) ;
                Ds(:,:,i) =  D;
              
                eigA(:,i+1:size(eigA,2)+1) = eigA(:,1:end) ;
                eigA(:,i) = eigA_n;
               
                Aas(:,:,i+1:size(Aas,3)+1) = Aas(:,:,1:end);
                Aas(:,:,i) =  [A, B*u;
                                 zeros(1,size(A,1)+1)];
                
                
                i = i+1;
                IC = x(end,:);
                clear x
                x = IC;
                
                
                
            end
            
        end
        
        
        
    end
    Xss = x(end,1:end-1)';
    
    
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
    ME.stack.line
    rethrow(ME)
end

end
