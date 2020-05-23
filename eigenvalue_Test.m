
clear

load('Test2.mat')
number_selected = 2;

%sysss = ss(obj.As(:,:,6),[obj.Bs(:,:,6) eye(size(obj.As(:,:,6)))],eye(size(obj.As(:,:,6))),zeros(size([obj.Bs(:,:,6) eye(size(obj.As(:,:,6)))])));


sysss = ss(obj.As(:,:,number_selected),[obj.Bs(:,:,number_selected) eye(size(obj.As(:,:,number_selected)))],eye(size(obj.As(:,:,number_selected))),zeros(size([obj.Bs(:,:,number_selected) eye(size(obj.As(:,:,number_selected)))])));



sys = tf(sysss);
sys.InputName = [obj.Converter.Topology.Parser.ConstantNames ; obj.Converter.Topology.stateLabels];
sys.OutputName = obj.Converter.Topology.stateLabels;


figure(3)
h = bodeplot(sys(:,size(obj.Bs,2)+1:end));
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');

[Z,P,K]=zpkdata(sys);





n = 1000;
t = linspace(0,obj.ts(end),n);
%u = [5 0.2 0.2 0.2 0 0 0 0 0 0 0]';
u = [5 0.2 0.2 0 0 0]';

%lsim(sysss,u*ones(size(t)),t,obj.Xs(:,end-1)')
figure(1)
lsim(sysss,u*ones(size(t)),t,[7.83 -0.2 -0.000277])

Poles = P{end,end};

Poles(2) = real(Poles(2))*20+imag(Poles(2))*1i;
Poles(3) = real(Poles(3))*20+imag(Poles(3))*1i;
%Pole = num2cell(Poles);


for i = size(obj.Bs,2)+1:1:size(obj.Bs,2)+size(obj.As,2)
   for j = 1:1:size(obj.As,2)
       P{j,i} = Poles;
   end 
end

 P{1,1} = Poles;
 P{2,1} = Poles;
 P{3,1} = Poles;

 
new_sys = zpk(Z,P,K);

new_sys.InputName = [obj.Converter.Topology.Parser.ConstantNames ; obj.Converter.Topology.stateLabels];
new_sys.OutputName = obj.Converter.Topology.stateLabels;

figure(4)
h = bodeplot(new_sys(:,size(obj.Bs,2)+1:end));
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');


%% Nonlinear (bounded) least squares regression to find zero for Buck implementation

As = obj.As(:,:,number_selected);
Bs = obj.Bs(:,:,number_selected);
Cs = eye(size(obj.As(:,:,number_selected)));
ti = obj.ts(2);
u = [5 0.2 0.2]';
X0 = [7.83 -0.2 -0.000277]';
Stick=X0;


new_sys_ss = ss(new_sys);
As = new_sys_ss.A;
Bs = new_sys_ss.B;
Cs = new_sys_ss.C;
u = [5 0.2 0.2 0 0 0]';
X0 = zeros(size(As,1),1);

% Stick = [    3.8091;  1.1909;   -0.0091];
   
   
   
[X,RESNORM,RESIDUAL,EXITFLAG] = ...
    lsqnonlin(@(x) IC_Error(x, As, Bs, ti,u,Cs,Stick),X0);

figure(2)
lsim(new_sys_ss,u*ones(size(t)),t,X')

return

sys = ss(sys);


obj.Xs(:,end-1);

[I,J]=find(sys.C~=0);

X0 = sum(sys.C)';

i = 1;
%X0(1:i*7) = obj.Xs(i,end-1)/sys.C(I(i),J(i));

for i  = 1:1:7
    
    X0(i*7) = obj.Xs(i,end-1)/sys.C(I(i),J(i));
    
end

[Y,T,X] = lsim(sys,u*ones(size(t)),t,X0);

sysss_old = ss(obj.As(:,:,6),[obj.Bs(:,:,6)],eye(size(obj.As(:,:,6))),zeros(size([obj.Bs(:,:,6)])));

tf(sysss_old)

%% Small RLC

clear


syms R L C C1 C2 R1 R2 R3
Asyms = [0 -1/C2
    1/L R3/L];
eig(Asyms)





A33 = [-1/R1/C1 0 -1/C1
    0 -1/R2/C2 -1/C2
    1/L 1/L -R3/L];


ShortedL = 0.47382e-6;
ShortedR = 326.9e-3;
L = 1.0393e-6-ShortedL;  % Resonate indcutor
L = 1.0393e-6-ShortedL;  % Resonate indcutor
C1 = 220e-9; % Resonate Cap
R1 = 390;
R3 = 0.63521-ShortedR; % Inductor series resistance
R2 = 0.2834;
C2 = 28.47e-12;

eig(eval(A33))




function [err] = IC_Error(X0, As, Bs, ti,u,Cs,Stick)
% This function just accounts for error and runs the function below

%Xi = AugSS(X0,As,Bs,ti,u);

X_Try = Cs*X0;
err = (X_Try-Stick)./Stick;

end

function Xi = AugSS(X0,As,Bs,ti,u)

ns = size(As,2);
pi = size(Bs,2);
nsub = length(ti);

X0 = [X0; 1 ; zeros(ns,1)];

Aas = [As, Bs*u;
    zeros(1,ns+1)];

Aais = [Aas, zeros(ns+1, ns);
    eye(ns), zeros(ns, ns+1)];

Xi = expm(Aais*ti)*X0;

Xi(ns+1:end,:) = [];
 
end


