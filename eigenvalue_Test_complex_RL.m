clear


load('Test3.mat')
number_selected = 2;
 u = [5 0.2 0.2 0 0 0 0]';

%load('Test.mat')
%number_selected = 6;
%u = [5 0.2 0.2 0.2 0 0 0 0 0 0 0]';

%sysss = ss(obj.As(:,:,6),[obj.Bs(:,:,6) eye(size(obj.As(:,:,6)))],eye(size(obj.As(:,:,6))),zeros(size([obj.Bs(:,:,6) eye(size(obj.As(:,:,6)))])));


ShortedR = 326.9e-3;
R3 = 0.63521-ShortedR; % Inductor series resistance


M1_R = 0.054;
M2_R = 0.0305;
M3_R = 0.2834;
R9 = 0.001; % input resistance 
L20 = 1e-3;
R20 = 6.2-R3;


Ress_out = [1/10000000 0 0 0
    0 1/10000000 0 0
    0 0 R3 0
    0 0 R20 -R20
    -1 -1 0 0
    -1/R9 -1/R9 0 0];

Bx = [obj.Bs(:,:,number_selected) eye(size(obj.As(:,:,number_selected)))];
Cx = [eye(size(obj.As(:,:,number_selected)));Ress_out];
Dx = zeros(size(Cx,1), size(Bx,2));
Dx(end,1) = 1/R9;
Dx(end-1,1) = 1;


sysss = ss(obj.As(:,:,number_selected),Bx,Cx,Dx);



sys = tf(sysss);
sys.InputName = [obj.Converter.Topology.Parser.ConstantNames ; obj.Converter.Topology.stateLabels];
sys.OutputName = {obj.Converter.Topology.stateLabels{:} , 'M1_R [I]', 'M3_R [I]', 'L_R [V]','L_Rac [V]'  , 'Rg [V]', 'Rg [A]' };


figure(3)
h = bodeplot(sys(:,size(obj.Bs,2)+1:end));
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');

[Z,P,K]=zpkdata(sys);


n = 1000000;
t = linspace(0,obj.ts(end),n);

%lsim(sysss,u*ones(size(t)),t,obj.Xs(:,end-1)')
% Need code that will find certain poles and change them 

figure(1)
Y1= lsim(sysss,u*ones(size(t)),t,[7 -2 -0.000277 -0.000277]);
lsim(sysss,u*ones(size(t)),t,[7 -2 -0.000277 -0.000277]);


sysss_plot = ss( sysss.A, sysss.B, sysss.C(1:4,:), sysss.D(1:4,:));

figure(5)
lsim(sysss_plot,u*ones(size(t)),t,[7 -2 -0.000277 -0.000277]);


%{
Poles = P{end,end};
Poles(2) = real(Poles(2))*20+imag(Poles(2))*1i;
Poles(3) = real(Poles(3))*20+imag(Poles(3))*1i;
%Pole = num2cell(Poles);


for i = size(obj.Bs,2)+1:1:size(obj.Bs,2)+size(obj.As,2)
   for j = 1:1:size(Cx,1)
       P{j,i} = Poles;
   end 
end

 P{1,1} = Poles;
 P{2,1} = Poles;
 P{3,1} = Poles;
 P{4,1} = Poles;
 P{5,1} = Poles;
 P{6,1} = Poles;
 P{7,1} = Poles;
 P{8,1} = Poles;


 
new_sys = zpk(Z,P,K);

new_sys.InputName = [obj.Converter.Topology.Parser.ConstantNames ; obj.Converter.Topology.stateLabels];
new_sys.OutputName = {obj.Converter.Topology.stateLabels{:} , 'M1_R [I]', 'M3_R [I]', 'L_R [V]', 'Rg [V]', 'Rg [A]'};

figure(4)
h = bodeplot(new_sys(:,size(obj.Bs,2)+1:end));
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');




%% Nonlinear (bounded) least squares regression to find zero for Buck implementation

As = obj.As(:,:,number_selected);
Bs = obj.Bs(:,:,number_selected);
Cs = [eye(size(obj.As(:,:,number_selected)));Ress_out];
Ds = Dx ; 
u = [5 0.2 0.2]';
X0 = [7 -2 -0.000277]';
ti = 5e-9;
Xi = AugSS(X0,As,Bs,ti,u);
Xi2 = AugSS(X0,As,Bs,ti*2,u);
Stick0 = Cs*X0+Dx(:,1:3)*u;
Stick1 = Cs*Xi+Dx(:,1:3)*u;
Stick2 = Cs*Xi2+Dx(:,1:3)*u;

Stick=[Stick0 Stick1 Stick2];


new_sys_ss = ss(new_sys);
As = new_sys_ss.A;
Bs = new_sys_ss.B;
Cs = new_sys_ss.C;
Ds = new_sys_ss.D;
u = [5 0.2 0.2 0 0 0]';
X0 = zeros(size(As,1),1);

%Stick = [3.8091;  1.1909;   -0.0091];


options = optimoptions('lsqnonlin','Display','iter');


   
[X,RESNORM,RESIDUAL,EXITFLAG] = ...
    lsqnonlin(@(x) IC_Error(x, As, Bs, ti,u,Cs,Ds,Stick),X0,[],[],options);







figure(2)
lsim(new_sys_ss,u*ones(size(t)),t,X');
Y2 = lsim(new_sys_ss,u*ones(size(t)),t,X');
%}
%{
R_before = Y1(:,6)./Y1(:,3); 
R_after = Y2(:,6)./Y2(:,3);
R_g_after = Y2(:,7)./Y2(:,8);
R_g_before = Y1(:,7)./Y1(:,8);

P_R_g_after = mean(Y2(:,7).*Y2(:,8));
P_R_g_before = mean(Y1(:,7).*Y1(:,8));


P_before = mean(Y1(:,6).*Y1(:,3));
P_after = mean(Y2(:,6).*Y2(:,3));


P_source_before = -mean(5.*Y1(:,8));
P_source_after = -mean(5.*Y2(:,8));

P_M1_before = mean(Y1(:,1).*Y1(:,4));
P_M1_after = mean(Y2(:,1).*Y2(:,4));
P_M2_before = mean(Y1(:,2).*Y1(:,5));
P_M2_after = mean(Y2(:,2).*Y2(:,5));


L = 5.6548e-07;
M1 = 1.6600e-10;
D2 = 2.8470e-11;

EM1_start = 0.5*M1*7^2;
ED2_start = 0.5*D2*(-2)^2;
EL_start = 0.5*L*(-0.000277)^2;

E_start = EM1_start + ED2_start + EL_start;

EM1_before = 0.5*M1*Y1(end,1)^2;
ED2_before = 0.5*D2*Y1(end,2)^2;
EL_before = 0.5*L*Y1(end,3)^2;

E_before = EM1_before + ED2_before + EL_before;

EM1_after = 0.5*M1*Y2(end,1)^2;
ED2_after = 0.5*D2*Y2(end,2)^2;
EL_after = 0.5*L*Y2(end,3)^2;

E_after = EM1_after + ED2_after + EL_after;

Power_change_before = (E_start - E_before) / obj.ts(end);
Power_change_after = (E_start - E_after) / obj.ts(end);

Add_before =  P_source_before + P_M1_before + P_M2_before + P_before + P_R_g_before ;
Add_after = P_source_after + P_M1_after + P_M2_after + P_after + P_R_g_after ;


%}


%% The stuff for RL in parallel



R_before = Y1(:,7)./Y1(:,3); 
R_g_before = Y1(:,9)./Y1(:,10);



P_R_g_before = mean(Y1(:,9).*Y1(:,10));
P_before = mean(Y1(:,7).*Y1(:,3));

Pac_before = mean(Y1(:,8).*(Y1(:,3)-Y1(:,4)));


P_source_before = -mean(5.*Y1(:,10));


P_M1_before = mean(Y1(:,1).*Y1(:,5));

P_M2_before = mean(Y1(:,2).*Y1(:,6));




L = 5.6548e-07;
M1 = 1.6600e-10;
D2 = 2.8470e-11;
L20 = 1e-3;


EM1_start =  0.5* M1*   7          ^2;
ED2_start =  0.5* D2* (-2)         ^2;
EL_start =   0.5* L*  (-0.000277)  ^2;
ELac_start = 0.5* L20*(-0.000277)  ^2;

E_start = EM1_start + ED2_start + EL_start+ELac_start;

EM1_before =  0.5* M1  * Y1(end,1) ^2;
ED2_before =  0.5* D2  * Y1(end,2) ^2;
EL_before  =  0.5* L   * Y1(end,3) ^2;
ELac_before = 0.5* L20 * Y1(end,4) ^2;



E_before = EM1_before + ED2_before + EL_before +ELac_before;

Power_change_before = (E_start - E_before) / obj.ts(end);

Add_before =  P_source_before + P_M1_before + P_M2_before + P_before + P_R_g_before + Pac_before;





return


%% This is all old stuff


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





%}



function [err] = IC_Error(X0, As, Bs, ti,u,Cs,Ds,Stick)
% This function just accounts for error and runs the function below


Xi = AugSS(X0,As,Bs,ti,u);
Xi2 = AugSS(X0,As,Bs,ti*2,u);
X_Try = Cs*X0+Ds*u;
%X_Try = Cs*X0;
err = [ ((X_Try-Stick(:,1))./Stick(:,1))' (((Cs*Xi+Ds*u)-Stick(:,2))./Stick(:,2))' (((Cs*Xi2+Ds*u)-Stick(:,3))./Stick(:,3))'];

if size(X_Try,1)>6
err(7) = (X_Try(7)-Stick(7,1));
err(8) = (X_Try(8)-Stick(8,1));
end


end

function Xi = AugSS(X0,As,Bs,ti,u)

ns = size(As,2);
pi = size(Bs,2);
nsub = length(ti);

if isempty(u)
    Xi = expm(As*ti)*X0;
    return
end



X0 = [X0; 1 ; zeros(ns,1)];

Aas = [As, Bs*u;
    zeros(1,ns+1)];

Aais = [Aas, zeros(ns+1, ns);
    eye(ns), zeros(ns, ns+1)];

Xi = expm(Aais*ti)*X0;

Xi(ns+1:end,:) = [];
 
end
