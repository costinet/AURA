
% This code is for the EPC 9118 converter components to run in
% simulink


Vg = 48;
Cin = 4*4.7e-6+4*1e-6;
Cout = 4*47e-6;
L = 1.2e-6;
fs = 400e3;
Q1_Cds = 430e-12;
Q1_Ron = 7e-3;

RL = 1e-3;

Q2_Cds = 1100e-12;
Q2_Ron = 2.2e-3;

%Rload = 5/5;
%Rload = 5/10;
%Rload = 5/15;
Rload = 5/20;


% Duty cycle is in pu
% Phase delay is in seconds

% these are in seconds
ON_time = 265e-9;
Primary_Dead = 37.5e-9;
Secondary_Dead = 31.25e-9;

Q1_dc = ON_time*fs;
Q1_pd = 0;


Q2_dc = (1/fs-ON_time-Secondary_Dead-Primary_Dead)*fs ;
Q2_pd = ON_time+Primary_Dead;