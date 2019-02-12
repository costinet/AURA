function [m] = binary_search(~,A,n,T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L = 1;
R = n;
while L <= R
    m = ceil((L+R)/2);
    if A(m) < T
        L = m + 1;
    elseif A(m) > T
        R = m - 1;
    else
        return 
    end
end
m = NaN;
end

