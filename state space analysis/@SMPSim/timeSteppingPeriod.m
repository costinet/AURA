function [Xf,ts,swinds] = timeSteppingPeriod(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 
    [~, deltaTs] = getDeltaT(obj.converter, 1:length(obj.ts));
    tVecLength = sum(ceil(obj.converter.ts./deltaTs),'all');
% 
%     xs = zeros(size(obj.Xs,1),tVecLength);
%     xs(:,1) = obj.Xs(:,1);
%     for i = 2:length(obj.converter.ts)
% %         xs(:,i) = xs(:,i-1) + expm( ...
% %             IS it in violation?
% 
%     end



    i=1;
    xs = zeros(size(obj.Xs,1),tVecLength);
    xs(:,1) = obj.Xs(:,1);
    tvec = zeros(tVecLength,1);
    tind = 2;
    try
    while i < length(obj.ts)
        [~, deltaT] = getDeltaT(obj.converter, i);
        A = obj.converter.As(:,:,obj.converter.swind(i));
        B = obj.converter.Bs(:,:,obj.converter.swind(i));
        Cbnd = obj.converter.topology.Cbnd(:,:,obj.converter.swind(i));
        Dbnd = obj.converter.topology.Dbnd(:,:,obj.converter.swind(i));
        
        
        while max(tvec) < sum(obj.ts(1:i))
            if max(tvec) + deltaT > sum(obj.ts(1:i))
                %Make sure last step lines up on an interface
                deltaT = (max(tvec) + deltaT) - sum(obj.ts(1:i));
            end
            M = [A, B*obj.u];
            M = [M; zeros(diff(size(M)),size(M,2))];
            xtilde = [xs(:,tind-1); 1];
            
            nextX = expm(M*deltaT)*xtilde;
            nextX = nextX(1:end-1);
            
%             nextX =  expm(A*deltaT)*xs(:,tind-1) + A\(expm(A*deltaT)-eye(size(A)))*B*obj.u;
            violateMarginEnd = Cbnd*nextX + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1) + obj.converter.topology.bndHyst(:,2);
            
            if any(violateMarginEnd < 0)
               %go back and find error 
               xdotStart = A*xs(:,tind-1) + B*obj.u;
               vErrDotStart = Cbnd*xdotStart;
               targetValStart= Cbnd*xs(:,tind-1) + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1);
               
               xdotEnd = A*nextX + B*obj.u;
               vErrDotEnd = Cbnd*xdotEnd;
               targetValEnd = Cbnd*nextX + Dbnd*obj.u - obj.converter.topology.bndHyst(:,1);
               
               projErrTime = [targetValStart./vErrDotStart, -targetValEnd./vErrDotEnd];
               newDeltaT = mean(projErrTime(violateMarginEnd<0,:),'all');
               
               warning('then edit the sequence');
               
               tvec(tind) = tvec(tind-1)+deltaT;
               tind = tind+1
            else 
                xs(:,tind) = nextX;
                tvec(tind) = tvec(tind-1)+deltaT;
                tind = tind+1
            end
        end
        i = i+1
    end
    catch e 
        disp(e.message);
        asdfsa=1
    end
    
    Xf = xs;
    ts = tvec;
    swinds = 1;

end

