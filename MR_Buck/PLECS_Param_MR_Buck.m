%{

This code forms the parameters for the Active Clamp flyback model
found in the reference listed below:

R. Perrin, N. Quentin, B. Allard, C. Martin and M. Ali,
"High-Temperature GaN Active-Clamp Flyback Converter With Resonant
Operation Mode," in IEEE Journal of Emerging and Selected Topics in
Power Electronics, vol. 4, no. 3, pp. 1077-1085, Sept. 2016.

Modifications have been made to account for a 5 to 20 V 3W power
supply for a gate driver based on the topology.

%}

clear

 Vg = 12;
        L1 = 1e-6; %L
        C1 = 1*10^-6; %Cout
        L2 = 0.01e-9; % Resonate inductor
        fs = 1e6;
        Ts = 1/fs;
        V = 3;
        Io = 0.1; % was 1
        M1_C = 1.5e-9; % CHS
        M2_C = 1.5e-9; % LHS

        M1_R_ON = 0.002; % ronHS
        M2_R_ON = 0.002; % ronLS

duty = 0.25;





ModelPath = 'Flyback_PLECS_Model_EVAL/Circuit';

% Switch vector is in the order of [Rec_D Snub_D FET_D Fet]
swVec(1,:) = [1 0 0 0];

ssOrder = plecs('get', ModelPath, 'StateSpaceOrder');
plecs('set', ModelPath, 'SwitchVector', swVec(1,:));
names = plecs('get', ModelPath, 'Topology');
As(:,:,1) = names.A;
Bs(:,:,1) = names.B;
Cs(:,:,1) = names.C;
Ds(:,:,1) = names.D;
Is(:,:,1) = names.I;

eigA = eig(As(:,:,1));



