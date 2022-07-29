function  [ron,Coss,w,l]=Select_FET(FET_selection)


switch FET_selection
    
    %% EPC 2023
    case 1
        
        ron = 1.45e-3;
        Coss = 1530e-12;
        w = 6.05;
        l = 2.3;
        
        %% EPC 2014C
    case 2
        
        ron = 16e-3;
        Coss = 1530e-12;
        w = 1.7;
        l = 1.1;
        
        %% EPC 2015C
    case 3
        
        ron = 4e-3;
        Coss = 710e-12;
        w = 4.1;
        l = 1.6;
        
        %% EPC 2055
    case 4
        
        ron = 3.6e-3;
        Coss = 408e-12;
        w = 2.5;
        l = 1.5;
        
        %% EPC 2030
    case 5
        
        ron = 2.4e-3;
        Coss = 1120e-12;
        w = 4.6;
        l = 2.6;
        
        %% EPC 2024
    case 6
        
        ron = 1.5e-3;
        Coss = 1620e-12;
        w = 6.05;
        l = 2.3;
        
        %% EPC 2031
    case 7
        
        ron = 2.6e-3;
        Coss = 980e-12;
        w = 4.6;
        l = 2.6;
        
        
        %% EPC 2020
    case 8
        
        ron = 2.2e-3;
        Coss = 1020e-12;
        w = 6.05;
        l = 2.3;
        
        %% EPC 2065
    case 9
        
        ron = 3.6e-3;
        Coss = 534e-12;
        w = 3.5;
        l = 1.95;
        
        %% EPC 2029
    case 10
        
        ron = 3.2e-3;
        Coss = 820e-12;
        w = 4.6;
        l = 2.6;
        
        
        %% EPC 2021
    case 11
        
        ron = 2.2e-3;
        Coss = 1100e-12;
        w = 6.05;
        l = 2.3;
        
        %% EPC 2059
    case 12
        
        ron = 9e-3;
        Coss = 267e-12;
        w = 2.8;
        l = 1.4;
        
        %% EPC 2052
    case 13
        
        ron = 13.5e-3;
        Coss = 195e-12;
        w = 1.5;
        l = 1.5;
        
        %%  EPC 2204
    case 14
        
        ron = 6e-3;
        Coss = 304e-12;
        w = 2.5;
        l = 1.5;
        
        %%  EPC 2218
    case 15
        
        ron = 3.2e-3;
        Coss = 562e-12;
        w = 3.5;
        l = 1.95;
        
        %%  EPC 2302
    case 16
        
        ron = 1.8e-3;
        Coss = 1000e-12;
        w = 3;
        l = 5;
        
        
        %%  EPC 2069
    case 17
        
        ron = 2.25e-3;
        Coss = 1044e-12;
        w = 3.25;
        l = 3.25;
        
        
        %%  EPC 2067
    case 18
        
        ron = 1.55e-3;
        Coss = 1071e-12;
        w = 2.85;
        l = 3.25;
        
        %%  EPC 2066
    case 19
        
        ron = 1.1e-3;
        Coss = 1670e-12;
        w = 6.05;
        l = 2.3;
        
        
end

end