%% Test Cases adding and subtracting devices to transistor database

clear

topDB = topologyDB();

%% SCBuckHybridBuck
SCBuckHybridBuck = topology('SCBuckHybridBuck Feb28',1,0,'SCBuckHybridBuck.net');

topDB.add(SCBuckHybridBuck)
%% Dual_Phase_MultiInductor_Hybrid
Dual_Phase_MultiInductor_Hybrid = topology('Dual_Phase_MultiInductor_Hybrid Feb28',1,0,'Dual_Phase_MultiInductor_Hybrid.net');

topDB.add(Dual_Phase_MultiInductor_Hybrid)


%% BuckHybridBuck
BuckHybridBuck = topology('BuckHybridBuck Feb28',1,0,'BuckHybridBuck.net');


topDB.add(BuckHybridBuck)

%% 3LevelBuck
COMPEL2023_3LevelBuck = topology('3LevelBuck Feb28',1,0,'COMPEL2023_3LevelBuck.net');


topDB.add(COMPEL2023_3LevelBuck)

%% 4LevelBuck
COMPEL2023_4LevelBuck = topology('4LevelBuck Feb28',1,0,'COMPEL2023_4LevelBuck.net');

topDB.add(COMPEL2023_4LevelBuck)

%% Buck
COMPEL2023_Buck = topology('Buck',1,0,'Buck_Vout_D.net');

topDB.add(COMPEL2023_Buck)


topDB.saveDB