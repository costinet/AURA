function [arrSat] = sat(arr,min,max)

if(nargin < 2)  % check number of args
    min = 0;    % decide on graph limits
    max = 1; 
    disp('[Too few input arguments]: min & max set to [0, 1] in function: sat'); 
elseif(nargin < 3)
    max = min + 1; 
    disp(['[Too few input arguments]: max set to ' num2str(max) ' in function: sat']); 
end 

arr(arr>max) = max; 
arr(arr<min) = min; 

arrSat = arr; 

end