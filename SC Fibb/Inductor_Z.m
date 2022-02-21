

%{

         |-R1-|
--R2--L2-|-C1-|--
         |-L1-|

%}



R1 = 35;
L1 = 300e-9;
C1 = 375e-11;




R2 = 0.01;
L2 = 900e-9;
omega = 5000000*2*pi;

%{
for f = 1:1:10000000

w = 2*pi*f;

Z(f) = R2+L2+1/(1/R1+1/(1j*w*L1)+1j*w*C1);

absZ(f) = abs(Z(f));
angZ(f) = angle(Z(f));

end

%}

s = tf('s');

Z = R2+s*L2+1/(1/R1+1/(s*L1)+s*C1);
%Z =  1/(1/(R1+s*L1)+s*C1);

plotoptions = bodeoptions;
plotoptions.FreqUnits = 'Hz';
plotoptions.FreqScale = 'linear';
plotoptions.MagUnits = 'abs';
plotoptions.Xlim = [1,10000000];
bodeplot(Z,plotoptions)
w = 1/omega;