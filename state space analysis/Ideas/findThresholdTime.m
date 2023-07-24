function t = findThresholdTime(Am, xt, row, x0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% Just a test case
if nargin == 0
    load('testSyncBuck', 'As', 'Bs', 'u', 'ts', 'Xss', 'xs', 't');
    A = As(:,:,2);
    B = Bs(:,:,2);
    tover = ts(2);
    
    xt = [0 48];
    row = [3 4];
    x0 = Xss(:,2);
    
    Am = [A, B*u; zeros(1,5)];
    
    xresult = [xs(:,1123); 1];
    tresult = t(1123)-ts(1);
end
    
    


%% basis for solution space
% finds basis for solutions for y in expm(Am*t0)*y == [~ ... ~ xt ~ ... ~ 1]
M = expm(Am*1e-12);
Xb = M([row end],:)\[xt'; 1];
Zb = null(M([row end],:));

coeffs = 100*ones(size(Zb,2),1);

%coeffs aren't correct -- need to be found

Xf = M*(Zb*coeffs + Xb)


%% let's try an assumption -- other states don't change
Xf = [x0; 1];
Xf([row end]) = [xt'; 1];

Xb = M\Xf;

%% if we did have it...
coeffs = Zb\(expm(-Am*1e-16)*xresult - Xb);
 M*(Zb*coeffs + Xb)
 xresult




%% Eigenvalue decomp
[V,D] = eigs(M);
% M = V*D/V
% M^N = V*D^N/V

%% Need to solve this somehow
% all possible things I can multiple M by to zero the state in "row" ==
% all possible state combinations from the system
Xf == M^N*x0;
(Zb*coeffs + Xb) == M^(N-1)*x0

M*(Zb*coeffs + Xb) ==  V*D^N/V*x0;

t = Am\logm(M^N)+1e-6;

end

