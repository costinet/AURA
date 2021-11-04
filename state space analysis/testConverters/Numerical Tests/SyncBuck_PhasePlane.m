L = 100e-6;
C = 1e-6;
R = 10;
Vg = 50;

A = [0, -1/L;
     1/C, -1/R/C];
B = [1/L; 0];

ilrange = linspace(-1, 10, 10);
vorange = linspace(0, 70, 15);

[x1 x2] = meshgrid(ilrange, vorange);

for i = 1:length(x1(:))
%     [i1, i2] = ind2sub(size(x1), i);
    xdot_ON(i,:) = A*[x1(i); x2(i)] + B*Vg;
    xdot_OFF(i,:) = A*[x1(i) + diff(ilrange(1:2))/2; x2(i)+diff(vorange(1:2))/2] + B*0;
end

Vb = Vg;
Ib = Vg/R;

hold off
quiver(x1(:)/Ib, x2(:)/Vb, xdot_ON(:,1)/Ib, xdot_ON(:,2)/Vb, 0.5);
hold on
quiver(x1(:)/Ib + diff(ilrange(1:2))/2/Ib, x2(:)/Vb+diff(vorange(1:2))/2/Vb, xdot_OFF(:,1)/Ib, xdot_OFF(:,2)/Vb, 0.5, 'r');
    
% quiver(x1(:)/Ib + diff(ilrange(1:2))/4/Ib, x2(:)/Vb+diff(vorange(1:2))/4/Vb, xdot_OFF(:,1)/Ib + xdot_ON(:,1)/Ib, xdot_OFF(:,2)/Vb + xdot_ON(:,2)/Vb, 0.5, 'm');
   

sys = ss(A, B, [1,0; 0 1], 0);

x0s = [0 0; 5 70; 5 0; 1, 70];
for i = 1:size(x0s,1)
    Y_ON = lsim(sys,Vg*ones(1,1000),linspace(0,1e-3,1000), x0s(i,:));
    plot(Y_ON(:,1)/Ib, Y_ON(:,2)/Vb, 'b', 'LineWidth', 3)
end

x0s = [5 50; 10 25; 7 0; 3 0];
for i = 1:size(x0s,1)
    Y_OFF = lsim(sys,zeros(1,1000),linspace(0,1e-3,1000), x0s(i,:));
    plot(Y_OFF(:,1)/Ib, Y_OFF(:,2)/Vb, 'r', 'LineWidth', 3)
end

xlim([0 10/Ib]);
ylim([0 70/Vb]);
xlabel('i_l(t)/V_g*R')
ylabel('v_o(t)/V_g')
legend('FET ON', 'FET OFF')

%% Solve Steady State
D = 0.5; Ts = 1/100e3;
X0 = (eye(2) - expm(A*Ts))\(expm(A*(1-D)*Ts)*(A\(expm(A*D*Ts)-eye(2))*B*Vg) ...
    + (A\(expm(A*(1-D)*Ts)-eye(2))*B*0));
ys1 = lsim(sys,Vg*ones(1,1000),linspace(0,D*Ts,1000), X0);
ys2 = lsim(sys,Vg*zeros(1,1000),linspace(0,(1-D)*Ts,1000), ys1(end,:));
ys = [ys1; ys2];
plot(ys(:,1)/Ib, ys(:,2)/Vb,'k','LineWidth',3)

%% Solve dynamic
X0 = [10; 0];
ys = X0';
for i = 1:20
    if(mod(i,2) == 1)
        ys1 = lsim(sys,Vg*ones(1,1000),linspace(0,Ts,1000), X0);
        ys = [ys; ys1(ys1(:,1) < 0.628*Ib,:)];
        X0 = ys(end,:);
    else
        ys1 = lsim(sys,zeros(1,1000),linspace(0,Ts,1000), X0);
        ys = [ys; ys1(ys1(:,1) > 0.372*Ib,:)];
        X0 = ys(end,:);   
    end
end
plot(ys(:,1)/Ib, ys(:,2)/Vb,':g','LineWidth',2)
