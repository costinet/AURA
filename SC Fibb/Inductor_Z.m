
%% OG USB Cable

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



%% 10ft 60W USB Cable

%{

         |-R1-|
--R2--L2-|-C1-|--
         |-L1-|

%}

% Biggest cable
% This is close but not exact yet 5.26.2022 JAB


R1 = 10;
L1 = 1522e-9;
C1 = 15e-08;

ideal_fs = 600000;

C_test = 1/L1/(ideal_fs*2*pi())^2;

fs = 1/2/pi()/sqrt(L1*C1);


R2 = 0.1;
L2 = 128e-9;


s = tf('s');

Z = R2+s*L2+1/(1/R1+1/(s*L1)+s*C1);
%Z =  1/(1/(R1+s*L1)+s*C1);

plotoptions = bodeoptions;
plotoptions.FreqUnits = 'Hz';
plotoptions.FreqScale = 'log';
plotoptions.MagUnits = 'abs';
plotoptions.MagScale =  'log';
plotoptions.Xlim = [10000,10000000];
bodeplot(Z,plotoptions)
[mag,phs,RadianFrequency] = bode(Z);
Magnitude = squeeze(mag);
Phase = squeeze(phs);
RadianFrequency = RadianFrequency./2./pi;
T = table(RadianFrequency,Magnitude,Phase);




%% 10ft 60W USB Cable

%{

         |-R1-|
--R2--L2-|-C1-|--
         |-L1-|

%}

% Biggest cable
% This is close but not exact yet 5.26.2022 JAB


R1 = 10;
L1 = 1522e-9;
C1 = 15e-08;

ideal_fs = 600000;

C_test = 1/L1/(ideal_fs*2*pi())^2;

fs = 1/2/pi()/sqrt(L1*C1);


R2 = 0.1;
L2 = 128e-9;


s = tf('s');

Z = R2+s*L2+1/(1/R1+1/(s*L1)+s*C1);
%Z =  1/(1/(R1+s*L1)+s*C1);

plotoptions = bodeoptions;
plotoptions.FreqUnits = 'Hz';
plotoptions.FreqScale = 'log';
plotoptions.MagUnits = 'abs';
plotoptions.MagScale =  'log';
plotoptions.Xlim = [10000,10000000];
bodeplot(Z,plotoptions)


