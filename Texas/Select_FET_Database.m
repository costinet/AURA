function  [ron,Coss,w,l]=Select_FET_Database(FET_selection)

FETs = Transistor.listDevices;
Chosen_FET = FETs(FET_selection);
Coss = Transistor.getTableParam(Chosen_FET{:},'Coss','Typ');
ron = Transistor.getTableParam(Chosen_FET{:},'Rds_on','Typ');
w = 1;
l = 1;