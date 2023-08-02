    load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')
    transDB = obj;
    clear obj
    
    load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductors.mat')
    indDB = obj;
    clear obj

%% SC Buck Hybrid Buck


%load('COMPEL_2023_SCBuckHybridBuck_D_Sweep_GA_saved.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions');

sav_x = [6.0000   13.0000   69.0000   67.0000   70.0000   14.0000    1.3404    0.8980
26.0000   34.0000   46.0000   46.0000   46.0000   50.0000    6.9742    3.7730
8.0000   28.0000   31.0000   50.0000   41.0000   50.0000    2.8200    1.4852
26.0000   35.0000   46.0000   48.0000   45.0000   51.0000    6.9697    1.7014
15.0000   29.0000   41.0000   49.0000   44.0000   50.0000    6.3618    2.4473
8.0000   28.0000   31.0000   50.0000   43.0000   50.0000    2.8071    1.4829
6.0000   13.0000   70.0000   68.0000   71.0000   15.0000    1.3345    0.8804
6.0000   18.0000   22.0000   61.0000   23.0000   60.0000    1.3244    0.6673
6.0000   18.0000   23.0000   61.0000   23.0000   60.0000    1.1890    0.6149
8.0000   28.0000   31.0000   50.0000   43.0000   50.0000    3.0216    1.7040
6.0000   18.0000   22.0000   61.0000   23.0000   60.0000    1.2722    0.6037
8.0000   27.0000   53.0000   49.0000   53.0000   35.0000    5.4564    2.2068
9.0000   28.0000   38.0000   49.0000   44.0000   50.0000    5.0299    2.2256
6.0000   13.0000   70.0000   68.0000   70.0000   15.0000    1.3381    0.8895
6.0000   13.0000   70.0000   68.0000   71.0000   15.0000    1.3349    0.8795
26.0000   35.0000   46.0000   48.0000   45.0000   51.0000    6.9714    1.5316
26.0000   35.0000   46.0000   48.0000   46.0000   50.0000    6.9693    1.6293
6.0000   13.0000   70.0000   67.0000   70.0000   15.0000    1.3381    0.8763];

stick = [16.1875    0.0751
    2.1136    0.3079
    9.9848    0.0838
    3.1354    0.1964
    5.0899    0.1100
   13.2899    0.0800
   22.4834    0.0680
   37.0785    0.0678
   42.9955    0.0663
   12.3462    0.0807
   38.6005    0.0677
    4.0663    0.1447
    7.4984    0.1040
   19.3645    0.0683
   22.4777    0.0680
    3.1346    0.1992
    2.3969    0.2552
   18.9733    0.0692
];


stick = stick(:,2);
sav_eff = 1-stick;


FETs = sav_x(:,1:5);
Inductor_Chosen = sav_x(:,end-2);
Iout = sav_x(:,end-1)*10;
fs = sav_x(:,end)*1e6;


% Number of FETs to parallel through
ParallelNumberM = 20;

% Number of Inductors to parallel through
ParallelNumberL = 20;

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];

parallel_fix = (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;

parallel_fix = parallel_fix + ones(size(parallel_fix));


for i = 1:size(FETs,1)
for j = 1:size(FETs,2)
    ron(i,j)= (1/parallel_fix(i,j))*transDB(FETs(i,j)).ron.typ*1e-3;
    Coss(i,j) = (parallel_fix(i,j))*transDB(FETs(i,j)).Coss.typ*1e-12;
    Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 3 * parallel_fix(i,j);
end
end

parallel_fixL = (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;

parallel_fixL = parallel_fixL + ones(size(parallel_fixL));

for j = 1:size(Inductor_Chosen,1)
RL(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).Rdc.typ*1e-3;
L(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).L.typ*1e-6;
Area_L(j) = (indDB(Inductor_Chosen(j)).width.typ) * (indDB(Inductor_Chosen(j)).length.typ) * 3 * parallel_fixL(j);
end


Area_L = Area_L';


scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff))

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
findfigs

% Color maps for each M1-3 (#||)
c = [];
c = sav_x(:,1)<=20;
c = c + (sav_x(:,1)>=21 & sav_x(:,1)<=40).*2;
c = c + (sav_x(:,1)>=41 & sav_x(:,1)<=60).*3;
c = c + (sav_x(:,1)>=61 & sav_x(:,1)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M1-3 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M1-3 (#)
c = sav_x(:,1);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M1-3 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs

% Color maps for each M4,8,9 (#||)
c = [];
c = sav_x(:,2)<=20;
c = c + (sav_x(:,2)>=21 & sav_x(:,2)<=40).*2;
c = c + (sav_x(:,2)>=41 & sav_x(:,2)<=60).*3;
c = c + (sav_x(:,2)>=61 & sav_x(:,2)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M4,8,9 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M4,8,9 (#)
c = sav_x(:,2);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M4,8,9 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs

% Color maps for each M5,10,13 (#||)
c = [];
c = sav_x(:,3)<=20;
c = c + (sav_x(:,3)>=21 & sav_x(:,3)<=40).*2;
c = c + (sav_x(:,3)>=41 & sav_x(:,3)<=60).*3;
c = c + (sav_x(:,3)>=61 & sav_x(:,3)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M5,10,13 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M5,10,13 (#)
c = sav_x(:,3);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M5,10,13 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M6,11,14 (#||)
c = [];
c = sav_x(:,4)<=20;
c = c + (sav_x(:,4)>=21 & sav_x(:,4)<=40).*2;
c = c + (sav_x(:,4)>=41 & sav_x(:,4)<=60).*3;
c = c + (sav_x(:,4)>=61 & sav_x(:,4)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M6,11,14 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M6,11,14 (#)
c = sav_x(:,4);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M6,11,14 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M7,12,15 (#||)
c = [];
c = sav_x(:,5)<=20;
c = c + (sav_x(:,5)>=21 & sav_x(:,5)<=40).*2;
c = c + (sav_x(:,5)>=41 & sav_x(:,5)<=60).*3;
c = c + (sav_x(:,5)>=61 & sav_x(:,5)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M7,12,15 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M7,12,15 (#)
c = sav_x(:,5);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front M7,12,15 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#||)
c = [];
c = sav_x(:,6)<=20;
c = c + (sav_x(:,6)>=21 & sav_x(:,6)<=40).*2;
c = c + (sav_x(:,6)>=41 & sav_x(:,6)<=60).*3;
c = c + (sav_x(:,6)>=61 & sav_x(:,6)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front L (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#)
c = sav_x(:,6);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front L (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each fs (#)
c = [];
c = sav_x(:,8);

figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('SC Buck Hybrid Buck Pareto Front fs','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs
set(gca,'ColorScale','log')

%{
hold on

scatter(10.12200, 0.0895,50,'r')

hold on

scatter(6.7582, 0.0889,50,'k')
%}
%% Buck Hybrid Buck


%load('COMPEL_2023_BuckHybridBuck_D_Sweep_GA_saved.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions');

sav_x = [7.0000   47.0000    9.0000    6.9939    2.0913
    8.0000   48.0000   35.0000    1.7978    0.5578
    8.0000   49.0000   13.0000    3.7423    0.6346
    8.0000   49.0000   34.0000    3.2075    0.7045
    6.0000   46.0000    9.0000    6.9939    2.0913
    8.0000   49.0000   38.0000    2.2510    0.4788
    8.0000   30.0000   38.0000    2.9308    0.5205
    8.0000   30.0000   38.0000    2.4605    0.4799
    8.0000   49.0000   35.0000    2.2776    0.6385
    8.0000   49.0000   35.0000    2.1475    0.5037
    8.0000   49.0000   14.0000    2.4514    0.7189
    8.0000   49.0000   34.0000    3.0783    0.6990
    8.0000   49.0000   35.0000    2.6541    0.6409
    8.0000   30.0000   38.0000    2.6345    0.4617
    8.0000   30.0000   38.0000    2.6853    0.5283
    8.0000   30.0000   35.0000    2.6689    0.5505
    7.0000   49.0000   29.0000    6.0766    1.1562
    7.0000   41.0000    9.0000    6.8805    1.1068
];

stick = [0.8665    0.1564
    7.1336    0.0610
    2.3999    0.0707
    4.2965    0.0663
    0.6700    0.2185
    9.2886    0.0562
    7.3798    0.0593
    8.7901    0.0572
    6.0506    0.0641
    6.4174    0.0613
    3.6637    0.0674
    4.4769    0.0660
    5.1924    0.0658
    8.2098    0.0574
    8.0544    0.0584
    5.4335    0.0647
    1.7724    0.0950
    1.0770    0.1150
];


stick = stick(:,2);
sav_eff = 1-stick;


FETs = sav_x(:,1:2);
Inductor_Chosen = sav_x(:,end-2);
Iout = sav_x(:,end-1)*10;
fs = sav_x(:,end)*1e6;


% Number of FETs to parallel through
ParallelNumberM = 20;

% Number of Inductors to parallel through
ParallelNumberL = 20;

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];
Area_FET = [];
Area_L = [];

parallel_fix = (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;

parallel_fix = parallel_fix + ones(size(parallel_fix));


for i = 1:size(FETs,1)
for j = 1:size(FETs,2)
    ron(i,j)= (1/parallel_fix(i,j))*transDB(FETs(i,j)).ron.typ*1e-3;
    Coss(i,j) = (parallel_fix(i,j))*transDB(FETs(i,j)).Coss.typ*1e-12;
    Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 3 * parallel_fix(i,j);
end
end

parallel_fixL = (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;

parallel_fixL = parallel_fixL + ones(size(parallel_fixL));

for j = 1:size(Inductor_Chosen,1)
RL(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).Rdc.typ*1e-3;
L(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).L.typ*1e-6;
Area_L(j) = (indDB(Inductor_Chosen(j)).width.typ) * (indDB(Inductor_Chosen(j)).length.typ) * 3 * parallel_fixL(j);
end


Area_L = Area_L';


scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff))

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
findfigs


% Color maps for each M1-3 (#||)
c = [];
c = sav_x(:,1)<=20;
c = c + (sav_x(:,1)>=21 & sav_x(:,1)<=40).*2;
c = c + (sav_x(:,1)>=41 & sav_x(:,1)<=60).*3;
c = c + (sav_x(:,1)>=61 & sav_x(:,1)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M1-3 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M1-3 (#)
c = sav_x(:,1);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M1-3 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs

% Color maps for each M4-6 (#||)
c = [];
c = sav_x(:,2)<=20;
c = c + (sav_x(:,2)>=21 & sav_x(:,2)<=40).*2;
c = c + (sav_x(:,2)>=41 & sav_x(:,2)<=60).*3;
c = c + (sav_x(:,2)>=61 & sav_x(:,2)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M4-6 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M4-6 (#)
c = sav_x(:,2);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M4-6 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#||)
c = [];
c = sav_x(:,3)<=20;
c = c + (sav_x(:,3)>=21 & sav_x(:,3)<=40).*2;
c = c + (sav_x(:,3)>=41 & sav_x(:,3)<=60).*3;
c = c + (sav_x(:,3)>=61 & sav_x(:,3)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front L (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#)
c = sav_x(:,3);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front L (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each fs (#||)
c = [];
c = sav_x(:,5);

figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front fs','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
set(gca,'ColorScale','log')
findfigs


%% Dickson2


%load('COMPEL_2023_Dickson2_Sweep_GA_saved.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions');


sav_x = [8.0000   46.0000   29.0000   18.0000    2.0061    0.5279
   13.0000    5.0000   67.0000   10.0000    6.7204    0.5890
    8.0000   46.0000   29.0000   15.0000    2.5047    0.5593
    8.0000   46.0000   29.0000   18.0000    2.1020    0.5310
    8.0000   46.0000   28.0000   14.0000    4.3230    0.7494
    8.0000   46.0000   29.0000   18.0000    2.1901    0.5316
    8.0000   46.0000   29.0000   15.0000    2.2548    0.5390
   11.0000    7.0000   61.0000    9.0000    6.1013    1.0042
    7.0000   46.0000   28.0000   14.0000    3.5280    0.5744
    8.0000   46.0000   29.0000   18.0000    2.3056    0.5291
    8.0000   46.0000   30.0000   14.0000    3.6436    0.6272
    7.0000   14.0000   54.0000    4.0000    6.4552    0.8307
    8.0000   46.0000   30.0000   18.0000    3.2196    0.5618
   13.0000    6.0000   66.0000    9.0000    6.1053    1.0047
    8.0000   46.0000   29.0000   15.0000    2.2569    0.5431
    8.0000   46.0000   29.0000   18.0000    2.6315    0.5392
    8.0000   46.0000   30.0000   14.0000    3.2591    0.5796
    8.0000   46.0000   30.0000   18.0000    3.0503    0.6089

];

stick = [17.2065    0.0395
    3.6690    0.0787
    9.5132    0.0452
   16.4215    0.0397
    5.2168    0.0596
   15.7607    0.0399
   10.5676    0.0439
    4.5025    0.0770
    6.3180    0.0578
   14.9716    0.0402
    7.5523    0.0461
    2.5037    0.1905
   11.8672    0.0429
    2.6516    0.1042
   10.5575    0.0439
   13.1173    0.0416
    8.4433    0.0454
   12.5262    0.0428
];


stick = stick(:,2);
sav_eff = 1-stick;


FETs = sav_x(:,1:3);
Inductor_Chosen = sav_x(:,end-2);
Iout = sav_x(:,end-1)*10;
fs = sav_x(:,end)*1e6;


% Number of FETs to parallel through
ParallelNumberM = 20;

% Number of Inductors to parallel through
ParallelNumberL = 20;

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];
Area_FET = [];
Area_L=[];

parallel_fix = (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;

parallel_fix = parallel_fix + ones(size(parallel_fix));


for i = 1:size(FETs,1)
    for j = 1:size(FETs,2)
        ron(i,j)= (1/parallel_fix(i,j))*transDB(FETs(i,j)).ron.typ*1e-3;
        Coss(i,j) = (parallel_fix(i,j))*transDB(FETs(i,j)).Coss.typ*1e-12;
        switch j
            case 1
                Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 6*parallel_fix(i,j);
            case 2
                Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 10*parallel_fix(i,j);
            case 3
                Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 6*parallel_fix(i,j);
            otherwise

        end
    end
end

parallel_fixL = (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;

parallel_fixL = parallel_fixL + ones(size(parallel_fixL));

for j = 1:size(Inductor_Chosen,1)
RL(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).Rdc.typ*1e-3;
L(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).L.typ*1e-6;
Area_L(j) = (indDB(Inductor_Chosen(j)).width.typ) * (indDB(Inductor_Chosen(j)).length.typ) * 9 * parallel_fixL(j);
end


Area_L = Area_L';


scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff))

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);



% Color maps for each M1-6 (#||)
c = [];
c = sav_x(:,1)<=20;
c = c + (sav_x(:,1)>=21 & sav_x(:,1)<=40).*2;
c = c + (sav_x(:,1)>=41 & sav_x(:,1)<=60).*3;
c = c + (sav_x(:,1)>=61 & sav_x(:,1)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front M1-6 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M1-6 (#)
c = sav_x(:,1);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front M1-6 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs

% Color maps for each M7-16 (#||)
c = [];
c = sav_x(:,2)<=20;
c = c + (sav_x(:,2)>=21 & sav_x(:,2)<=40).*2;
c = c + (sav_x(:,2)>=41 & sav_x(:,2)<=60).*3;
c = c + (sav_x(:,2)>=61 & sav_x(:,2)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front M7-16 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M7-16 (#)
c = sav_x(:,2);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front M7-16 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs



% Color maps for each M17-22 (#||)
c = [];
c = sav_x(:,3)<=20;
c = c + (sav_x(:,3)>=21 & sav_x(:,3)<=40).*2;
c = c + (sav_x(:,3)>=41 & sav_x(:,3)<=60).*3;
c = c + (sav_x(:,3)>=61 & sav_x(:,3)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front M17-12 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M17-22 (#)
c = sav_x(:,3);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front M17-22 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#||)
c = [];
c = sav_x(:,4)<=20;
c = c + (sav_x(:,4)>=21 & sav_x(:,4)<=40).*2;
c = c + (sav_x(:,4)>=41 & sav_x(:,4)<=60).*3;
c = c + (sav_x(:,4)>=61 & sav_x(:,4)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front L (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#)
c = sav_x(:,4);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front L (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each fs (#)
c = [];
c = sav_x(:,6);

figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dickson2 Pareto Front fs','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
set(gca,'ColorScale','log')
findfigs



%% LEGO


%load('COMPEL_2023_LEGO_GA_saved.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions');


sav_x = [   28.0000   48.0000   68.0000   52.0000    1.6905    0.3198
   28.0000   48.0000   49.0000   50.0000    3.3470    0.4593
   26.0000   14.0000   46.0000   13.0000    6.2559    1.0230
   28.0000   48.0000   68.0000   50.0000    4.6220    0.3858
   28.0000   48.0000   68.0000   51.0000    2.7808    0.3388
   28.0000   26.0000   48.0000   21.0000    6.2032    0.5894
   27.0000   15.0000   46.0000   14.0000    6.2559    1.0230
   28.0000   27.0000   61.0000   17.0000    5.4008    0.3405
   28.0000   48.0000   68.0000   51.0000    2.1076    0.3879
   28.0000   48.0000   68.0000   52.0000    1.7868    0.3225
   27.0000   46.0000   20.0000   13.0000    3.0520    0.5071
   28.0000   48.0000   68.0000   51.0000    2.4801    0.3642
   29.0000   48.0000   68.0000   50.0000    6.9382    0.3030
   28.0000   47.0000   20.0000   13.0000    3.0516    0.5071
   28.0000   16.0000   48.0000   17.0000    6.2076    0.6628
   28.0000   48.0000   68.0000   51.0000    3.7999    0.3872
   28.0000   48.0000   68.0000   52.0000    1.9408    0.3229
   27.0000   14.0000   46.0000   11.0000    6.4339    2.3469


];

stick = [16.4552    0.0720
    8.4007    0.0814
    1.5673    0.3082
    6.0186    0.0962
   10.0034    0.0791
    2.7393    0.1626
    1.9308    0.2495
    4.6453    0.1026
   13.1986    0.0742
   15.5679    0.0727
    6.4925    0.0928
   11.2162    0.0763
    4.1931    0.1241
    7.7909    0.0871
    3.8317    0.1379
    7.3206    0.0912
   14.3328    0.0741
    1.5137    0.3809
];


stick = stick(:,2);
sav_eff = 1-stick;



FETs = sav_x(:,1:3);
Inductor_Chosen = sav_x(:,end-2);
Iout = sav_x(:,end-1)*10;
fs = sav_x(:,end)*1e6;


% Number of FETs to parallel through
ParallelNumberM = 20;

% Number of Inductors to parallel through
ParallelNumberL = 20;

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];
Area_FET = [];
Area_L=[];

parallel_fix = (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;

parallel_fix = parallel_fix + ones(size(parallel_fix));


for i = 1:size(FETs,1)
    for j = 1:size(FETs,2)
        ron(i,j)= (1/parallel_fix(i,j))*transDB(FETs(i,j)).ron.typ*1e-3;
        Coss(i,j) = (parallel_fix(i,j))*transDB(FETs(i,j)).Coss.typ*1e-12;
        switch j
            case 1
                Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 6*parallel_fix(i,j);
            case 2
                Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 10*parallel_fix(i,j);
            case 3
                Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 6*parallel_fix(i,j);
            otherwise

        end
    end
end

parallel_fixL = (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;

parallel_fixL = parallel_fixL + ones(size(parallel_fixL));

for j = 1:size(Inductor_Chosen,1)
RL(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).Rdc.typ*1e-3;
L(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).L.typ*1e-6;
Area_L(j) = (indDB(Inductor_Chosen(j)).width.typ) * (indDB(Inductor_Chosen(j)).length.typ) * 3 * parallel_fixL(j);
end


Area_L = Area_L';


scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff))

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);



% Color maps for each M1-6 (#||)
c = [];
c = sav_x(:,1)<=20;
c = c + (sav_x(:,1)>=21 & sav_x(:,1)<=40).*2;
c = c + (sav_x(:,1)>=41 & sav_x(:,1)<=60).*3;
c = c + (sav_x(:,1)>=61 & sav_x(:,1)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front M1-6 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M1-6 (#)
c = sav_x(:,1);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front M1-6 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs

% Color maps for each M7-16 (#||)
c = [];
c = sav_x(:,2)<=20;
c = c + (sav_x(:,2)>=21 & sav_x(:,2)<=40).*2;
c = c + (sav_x(:,2)>=41 & sav_x(:,2)<=60).*3;
c = c + (sav_x(:,2)>=61 & sav_x(:,2)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front M7-16 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M7-16 (#)
c = sav_x(:,2);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front M7-16 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs



% Color maps for each M17-22 (#||)
c = [];
c = sav_x(:,3)<=20;
c = c + (sav_x(:,3)>=21 & sav_x(:,3)<=40).*2;
c = c + (sav_x(:,3)>=41 & sav_x(:,3)<=60).*3;
c = c + (sav_x(:,3)>=61 & sav_x(:,3)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front M17-22 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M17-22 (#)
c = sav_x(:,3);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front M17-22 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#||)
c = [];
c = sav_x(:,4)<=20;
c = c + (sav_x(:,4)>=21 & sav_x(:,4)<=40).*2;
c = c + (sav_x(:,4)>=41 & sav_x(:,4)<=60).*3;
c = c + (sav_x(:,4)>=61 & sav_x(:,4)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front L (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#)
c = sav_x(:,4);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front L (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each fs (#)
c = [];
c = sav_x(:,6);

figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('LEGO Pareto Front fs','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
set(gca,'ColorScale','log')
findfigs

%% Dual Phase Multi-Inductor Hybrid


load('COMPEL_2023_DualPhaseMultiIndHybrid_D_Sweep_GA_saved.mat','sav_x','sav_eff','sav_Ploss','sav_Pout','sav_conditions');


FETs = sav_x(:,1:2);
Inductor_Chosen = sav_x(:,end-2);
Iout = sav_x(:,end-1);
fs = sav_x(:,end);


% Number of FETs to parallel through
ParallelNumberM = 20;

% Number of Inductors to parallel through
ParallelNumberL = 20;

% FET Selection
ron = [];
Coss = [];
FET_w = [];
FET_l = [];
Area_FET = [];
Area_L=[];

parallel_fix = (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;
parallel_fix = parallel_fix + (FETs>ParallelNumberM);
FETs(FETs>ParallelNumberM) = FETs(FETs>ParallelNumberM)-ParallelNumberM;

parallel_fix = parallel_fix + ones(size(parallel_fix));


for i = 1:size(FETs,1)
for j = 1:size(FETs,2)
    ron(i,j)= (1/parallel_fix(i,j))*transDB(FETs(i,j)).ron.typ*1e-3;
    Coss(i,j) = (parallel_fix(i,j))*transDB(FETs(i,j)).Coss.typ*1e-12;
    Area_FET(i,j) = (transDB(FETs(i,j)).width.typ) * (transDB(FETs(i,j)).length.typ) * 3;
end
end

parallel_fixL = (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;
parallel_fixL = parallel_fixL + (Inductor_Chosen>ParallelNumberL);
Inductor_Chosen(Inductor_Chosen>ParallelNumberL) = Inductor_Chosen(Inductor_Chosen>ParallelNumberL)-ParallelNumberL;

parallel_fixL = parallel_fixL + ones(size(parallel_fixL));

for j = 1:size(Inductor_Chosen,1)
RL(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).Rdc.typ*1e-3;
L(j) = (1/parallel_fixL(j))*indDB(Inductor_Chosen(j)).L.typ*1e-6;
Area_L(j) = (indDB(Inductor_Chosen(j)).width.typ) * (indDB(Inductor_Chosen(j)).length.typ) * 3;
end


Area_L = Area_L';


scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff))

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Dual Phase Multi-Inductor Hybrid Pareto Front','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
findfigs


% Color maps for each M1-3 (#||)
c = [];
c = sav_x(:,1)<=20;
c = c + (sav_x(:,1)>=21 & sav_x(:,1)<=40).*2;
c = c + (sav_x(:,1)>=41 & sav_x(:,1)<=60).*3;
c = c + (sav_x(:,1)>=61 & sav_x(:,1)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M1-3 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M1-3 (#)
c = sav_x(:,1);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M1-3 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs

% Color maps for each M4-6 (#||)
c = [];
c = sav_x(:,2)<=20;
c = c + (sav_x(:,2)>=21 & sav_x(:,2)<=40).*2;
c = c + (sav_x(:,2)>=41 & sav_x(:,2)<=60).*3;
c = c + (sav_x(:,2)>=61 & sav_x(:,2)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M4-6 (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each M4-6 (#)
c = sav_x(:,2);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front M4-6 (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#||)
c = [];
c = sav_x(:,3)<=20;
c = c + (sav_x(:,3)>=21 & sav_x(:,3)<=40).*2;
c = c + (sav_x(:,3)>=41 & sav_x(:,3)<=60).*3;
c = c + (sav_x(:,3)>=61 & sav_x(:,3)<=80).*4;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front L (#||)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each L (#)
c = sav_x(:,3);
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
c(c>20) =  c(c>20) - 20;
figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front L (#)','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


% Color maps for each fs (#||)
c = [];
c = sav_x(:,5);

figure
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],c)

% Create ylabel
ylabel('(1-eff)','FontSize',20);

% Create xlabel
xlabel('mm/A','FontSize',20);


% Create title
title('Buck Hybrid Buck Pareto Front fs','FontSize',24);

axis1 = gca;
set(axis1,'FontName','Times New Roman','FontSize',14);
colorbar
findfigs


%% Final Comparison

figure
load('BuckHybridBuck_graph_5_11_2023.mat')
hold on
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],'magenta')
load('Dickson2_graph_5_11_2023.mat')
hold on
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],"cyan")
load('LEGO_graph_5_11_2023.mat')
hold on
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],"red")
load('DualPhaseMultiIndHybrid_graph_5_11_2023.mat')
hold on
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],"blue")
load('SCBuckHybridBuckGraph5_10_2023.mat')
hold on
scatter((sum(Area_FET,2)+Area_L)./Iout,(1-sav_eff),[],'yellow')



findfigs





