function [ Xs] = SS_Soln( As, Bs, ts, u, HCM)
%HCM is the half-cycle matrix which negates any states which are odd with
%respect to the switching period.  For full cycle, just use HCM = eye(n)

n = size(As,3);
ns = size(Bs,1);

if (nargin < 5)
    HCM = eye(ns);
end


expAs = [];
for i=1:n
    expAs(:,:,i) = expm(As(:,:,i)*ts(i));
end

cumProdExp = [];
for i=1:n
    cumProdExp(:,:,i) = eye(ns);
    for k = 1:i
        cumProdExp(:,:,i) = expAs(:,:,k)*cumProdExp(:,:,i);
    end
end

cumProdExpRev = [];
for i=1:n
    cumProdExpRev(:,:,i) = eye(ns);
    if(i<n)
        for k = i+1:n
            cumProdExpRev(:,:,i) = expAs(:,:,k)*cumProdExpRev(:,:,i);
        end
    end
end

RHSsum = zeros(ns,1);
for i=1:n
   RHSsum = RHSsum +  cumProdExpRev(:,:,i)*As(:,:,i)^-1*(expAs(:,:,i)-eye(ns))*(Bs(:,:,i)*u);
end


Xss = (HCM - cumProdExp(:,:,n))^-1*RHSsum;

Xs(:,1) = Xss;
for(i=1:n)
    Xs(:,i+1) = expm(As(:,:,i)*ts(i))*Xs(:,i) + As(:,:,i)^-1*(expAs(:,:,i)-eye(ns))*(Bs(:,:,i)*u);
end

if sum(isnan(Xss))
    ME = MException('resultisNaN:noSuchVariable', ...
                       'Waveform reconstruction resulted in NaN');
    throw(ME);
end




% J2 = Xs{2}(2)/u(1)*Ro;
% Mo2 = Xs{2}(3)/u(1);
% r2 = sqrt(J^2 + Mo^2)
% 
% ZVS2 = J2 > (1-M);
% 
% if(ZVS2)
%     tzvs = min((pi - acos(Mo/r) - acos((1-Mo)/r)/wo,ts{2});
%     ts{3} = ts{3} + tzvs;
% end


% %% Check and revise
% diodeConduct = 0;
% multires = 0;
%     Ro = 20.6236;
%     wo = 1.7186e+07;
% if((Xs{3}(1) < -2) || Xs{2}(2) < 0)
%     diodeConduct = 1;
% end
% if(ts{2}>2*pi*wo/4*1.2) 
%     multires = 1;
% end
% 
% ts2 = ts{2};
% if(diodeConduct || multires)
%     J = Xs{2}(2)/u(1)*Ro;
%     if(J<0) 
%         ts{2} = 0; 
%         ts{3} = ts{3} + ts2-ts{2};
%     else
%         Mo = Xs{2}(3)/u(1);
%         r = sqrt(J^2 + Mo^2);
% 
%         theta1 = pi - atan2(J,Mo);
%         theta2 = pi - acos(Mo/r) - acos((1-Mo)/r);
% 
%         ts{2} = min([theta1/wo, theta2/wo, ts{2}(1)]);
%         ts{3} = ts{3} + ts2-ts{2};
%     end
% end
% 
% diodeConduct2 = 0;
% multires2 = 0;
% 
% if((Xs{end}(1) > u(1)+2) || Xs{4}(2) > 0)
%     diodeConduct2 = 1;
% end
% if(ts{4}>2*pi*wo/4*1.2) 
%     multires2 = 1;
% end
% 
% ts4 = ts{4};
% if(diodeConduct2 || multires2)
%     J = Xs{4}(2)/u(1)*Ro;
%     if(J>0) 
%         ts{4} = 0; 
%         ts{1} = ts{1} + ts4-ts{4};
%     else
%         Mo = Xs{4}(3)/u(1);
%         r = sqrt(J^2 + (1-Mo)^2);
% 
%         theta3 = pi - acos((1-Mo)/r) - acos(Mo/r);
% 
%         ts{4} = min([theta3/wo, ts{4}(1)]);
%     end
% end
% 
% if(diodeConduct || multires || diodeConduct2 || multires2)
%     if ((ts{4}-ts4) > 100e-12) || ((ts{2}-ts2) > 100e-12) 
%         [ Xs, ts ] = SS_Soln( As, Bs, ts, u, HCM);
%         x = 1
%     end
% end



