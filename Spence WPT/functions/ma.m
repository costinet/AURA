function [m,a] = ma(num, print)

% return the magnitud and angle (deg.) of a complex number
m = abs(num); 
a = angle(num)*180/pi; 

if(nargin < 2)  % check number of args
    print = 1; 
end; 

if(print)
    dispMag = ['mag:   ', num2str(m)]; 
    dispAng = ['angle: ', num2str(a)]; 
    disp(dispMag); 
    disp(dispAng); 
end; 
end