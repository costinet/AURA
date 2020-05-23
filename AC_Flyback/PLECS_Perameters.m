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


        M1_R = 0.054;
        M2_R = 0.0305;
        M3_R = 0.2834;
        
        M1_C = 166e-12;
        M2_C = 272e-12;
        M3_C = 28.47e-12;




clear

P = 3;
Lr = 50e-9;
Lm = 5e-6;
Co = 10e-6;
Rds1 = 0.016;
Rds2 = 0.016;
Rds3 = 0.016;
Cds1 = 150e-12;
Cds2 = 150e-12;
Cds3 = 0;
duty = .53;
np = 1;
ns = 4;
Vo = 20;
Ro = 133.3333;
Vg = 5; 
fs = 1000e3;
deadtime = 0.02/fs;
alpha = duty;
n = ns/np;
Cr = 5e-9;
Vf1 = 1;
Vf2 = 1;
Vf3 = 0.2;



if Cds1 == Cds2
    Cds = Cds1;
end

i_Lr_to = 2*Cds*(Vg+(Vo/(n)))/deadtime;
Lm = (Vg*alpha^2/fs^2-2*Lr*((P/Vg/fs)+i_Lr_to*alpha/fs))/(2*(((P/Vg/fs)+i_Lr_to*alpha/fs)));

fs_check  = fs<=(1-alpha)/(2*pi()*sqrt(Lr*Cr));

Cr = ((1-alpha)/fs/2/pi())^2/Lr;

fs_diff = (1-alpha)/(2*pi()*sqrt(Lr*Cr))-fs;
Gain = (n*alpha)/(1-alpha)*(Lm/(Lm+Lr));

%{
ModelPath = 'Flyback_PLECS_Model/Circuit2';

swVec(1,:) = [0 0 0];

ssOrder = plecs('get', ModelPath, 'StateSpaceOrder');
plecs('set', ModelPath, 'SwitchVector', swVec(1,:));
names = plecs('get', ModelPath, 'Topology');
As(:,:,1) = names.A;
Bs(:,:,1) = names.B;
Cs(:,:,1) = names.C;
Ds(:,:,1) = names.D;
Is(:,:,1) = names.I;

eigA = eig(As(:,:,1));

%}


