% Kgfe method

%% Variables
%{
rho; % Wire effective ressitivity (Ohm-cm)

I % Total rms winding currents referred to primary (A)

lambda % Applied primary volt-seconds

P_tot % Allowed total power dissipation (W)

K_u % Winding fill factor (-)

beta % Core loss exponent (-)

Kfe % Core loss coefficent (W/(cm^3T^beta))

A_c % Core cross sectional area (cm^2)
W_A % Core window area (cm^2)
MLT % Mean length per turn (cm)
l_m % Magnetic path length (cm)
DELTA_B % Peak ac flux density (Tesla)
A_w % Wire areas (cm^2)

u_0 = 4*pi*10^-7 % Permeability of free space (H/m)

%}
%% Set of Equations
%%{
% Core size ineequality


Kgfe >= ((rho*lambda^2*I^2*Kfe^(2/beta))/(2*Ku*Ptot^((beta+2)/beta)))*10^8;


delta_B  = (10^8*((rho*lambda*I)/(2*Ku))*((MLT)/(W_A*A_c^3*l_m))*((1)/(beta*Kfe)))^(1/(beta+2));

n = (lembda/(2*DELTA_B*A_c))*10^4;

l_g = (u_0*A_c*n^2/L)*10^-4;

A_L = (L/n^2)*10^9;

B_max = DELTA_B+((L*I_dc)/(n*A_c))*10^4;

A_w <= K_u*W_A/n;

P_cu = (rho*n*MLT/A_w)*I^2;

P_fe = K_fe*(DELTA_B)^beta*A_c*L_m;

P_tot = P_cu+P_fe;







%}




