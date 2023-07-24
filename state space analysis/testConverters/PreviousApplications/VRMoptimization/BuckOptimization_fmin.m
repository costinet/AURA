% clear all;
clc;

debug = 0;

HSFET.ron = .05;
HSFET.Cp = 3.4874e-10;
LSFET = HSFET;

inductor.Rl = .01;
inductor.L = 230e-9;

Vg = 5;
Cout = 100e-6;
fs = 500e3;
V = 1.8;
Io = .01;

x0 = [Io, fs/1e6];
history = [];
fun = @(x) -BuckSimulation( inductor, HSFET, LSFET, x(2)*1e6, Vg, V, x(1), Cout);

%% Optimization over fs and Io
nonopt= 1;
while nonopt
    % xmin = [0 100e3/1e6];
    % xmax = [5 10e6/1e6];
    [ x, fval, history_iter ] = fminHistory(fun, x0);
    history = [history; history_iter];
    % [ x, fval, history ] = fminconHistory(fun, x0, xmin, xmax);

    Iopt = x(1); dI = max(Iopt/100, .1);
    fsopt = x(2)*1e6; dfs = max(fsopt/1e3, 10e3);
    [eta_opt, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fsopt, Vg, V, Iopt, Cout);
    [eta_dI, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fsopt, Vg, V, Iopt+dI, Cout);
    [eta_dfs, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fsopt+dfs, Vg, V, Iopt, Cout);
    [eta_dIm, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fsopt, Vg, V, Iopt-dI, Cout);
    [eta_dfsm, ~, ~, ~, ~, ~, exitflag ] = BuckSimulation( inductor, HSFET, LSFET, fsopt-dfs, Vg, V, Iopt, Cout);

    partials = ([eta_dI, eta_dfs; eta_dIm, eta_dfsm] - eta_opt)./[dI dfs/1e6];
    if(sum(sum(sign(partials))) == -4)
        disp('local optimum found')
        nonopt = 0;
    else
        nonopt = 1;
        x0 = x + (partials(1,:) - partials(2,:) )/2.*x;
    end
end


load('bruteForceIofs_smallRes.mat');
% %% Brute force for checking
% Imax = 2;
% Iorange = logspace(-3,log10(Imax),100);
% fsmax = 10e6;
% fsrange = logspace(5, log10(fsmax), 100);
% 
% etas = zeros(length(Iorange), length(fsrange));
% % i = 26; j = 50;
% for i = 1:length(Iorange)
%     for j = 1:length(fsrange)
%         eta_iter = BuckSimulation( inductor, HSFET, LSFET, fsrange(j), Vg, V, Iorange(i), Cout);
%         etas(i,j) = eta_iter;
%     end
%     disp(i)
% end

%% i, j that didnt work
% broken = [];
%  row = 13;
% i = broken(row,1);
% j = broken(row,2);
% eta_iter = BuckSimulation( inductor, HSFET, LSFET, fsrange(j), Vg, V, Iorange(i), Cout);
   eta_iter = BuckSimulation( inductor, HSFET, LSFET, 3.6487*1e6, Vg, V, 0.2797 , Cout); 

% F = scatteredInterpolant(Ios', fss'/1e6, etas');
figure(1)
loglog(0,0); hold on;
[X, Y] = meshgrid(Iorange, fsrange/1e6);
contourf(X,Y, etas', 'Levellist', [.9894 .98 .97 .95 .93 .91 .89 .87 .85 .83 .81 .79 0]); caxis([.8 1]);
hold on;
scatter(history(:,1), history(:,2), [],-history(:,3));
colorbar; caxis([.8 1]);
xlabel('I_o');
ylabel('f_s [MHz]');
hold on;
plot(history(:,1), history(:,2), 'LineWidth',3)