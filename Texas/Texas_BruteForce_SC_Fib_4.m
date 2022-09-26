clear
sf = [1 1 1 1 1 1e-6 100 100];  
x = [3 3 3 3 3  1.2315e6    0.009    0.009];
x = x.*sf;
Coss_adj = [0 0 0 0 0 ];
Ron_adj = [0 0 0 0 0 ];

Tstart=tic;
Number_of_FETs = 15;
stick  = zeros(Number_of_FETs,Number_of_FETs,Number_of_FETs);
for i = 1:Number_of_FETs
    for j = 1:Number_of_FETs
        for k = 1:Number_of_FETs
            x(1:5) = [i j k 15 3];
            [stick(i,j,k),graph1]=Texas_SC_Fib_4_Sweep(x,Coss_adj,Ron_adj);
        end
    end
end


toc(Tstart)

    Number_of_FETs = 19;
    FETs = zeros(Number_of_FETs,2);
    for i = 1:Number_of_FETs
        [FETs(i,1),FETs(i,2),~,~]=Select_FET(i);
    end
    
F = scatteredInterpolant(FETs(:,1), FETs(:,2), stick', 'linear','none');


Ronrange = linspace(min(FETs(:,1)),max(FETs(:,1)));
CossRange = linspace(min(FETs(:,2)),max(FETs(:,2)));
[RonMesh, CossMesh] = meshgrid(Ronrange, CossRange);


%surf(Ronrange,CossRange,F(RonMesh,CossMesh))


% Create surf
surf(Ronrange,CossRange,F(RonMesh,CossMesh));

% Create zlabel
zlabel('Ploss (W)');

% Create ylabel
ylabel('Coss (F)');

% Create xlabel
xlabel('Ron (\Omega)');
axis1 = gca;
% Set the remaining axes properties
set(axis1,'FontName','Arial','FontSize',24);

hold on
scatter3(FETs(:,1),FETs(:,2),stick,100,'k','filled') ;





%{

Ronrange = linspace(min(FETs(:,1)),max(FETs(:,1));
CossRange = linspace(min(FETs(:,2)),max(FETs(:,2));
[RonMesh, CossMesh] = meshgrid(Ronrange, CossRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
c.LevelList = [0 1 2 3 4 5 6 7 8 9];






[VgMesh, PoutMesh] = meshgrid(Ronrange, CossRange);
[f,c] = contourf(Vgrange,PoutRange,F(VgMesh,PoutMesh),'ShowText','on');


% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create surf
surf(xdata1,ydata1,zdata1,cdata1,'Parent',axes1);

% Create zlabel
zlabel('Ploss (W)');

% Create ylabel
ylabel('Coss (F)');

% Create xlabel
xlabel('Ron (\Omega)');

view(axes1,[-3.57187499999999 13.7727272727273]);
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Arial','FontSize',24);

%}


