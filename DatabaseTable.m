%% Sorting Database based on Tables

clear
load('D:\GitHub\AURA\AURAdb\databases\@inductorDB\inductors.mat')

indDB = obj;
clear obj

load('D:\GitHub\AURA\AURAdb\databases\@transistorDB\transistors.mat')

transDB = obj;
clear obj

for i=1:length(indDB)
    NameI{i,1} = indDB(i).partNumber;
    L(i,1) = indDB(i).L.approx;
    Rdc(i,1) = indDB(i).Rdc.approx;
    Isat30(i,1) = indDB(i).Isat(3).approx;
    Length(i,1) = indDB(i).Length.approx;
    Width(i,1) = indDB(i).Width.approx;
    Height(i,1) = indDB(i).Height.approx;
end

IndTable = table(L,Rdc,Isat30,Length,Width,Height,...
    'RowNames',NameI);


UsableInd = NameI(IndTable.Isat30 >= 16);

IndTable(IndTable.Isat30 >= 16,:)



for i=1:length(transDB)
    NameT{i,1} = transDB(i).partNumber;
    Vds(i,1) = transDB(i).Vds.approx;
    Ids(i,1) = transDB(i).Ids.approx;
    Rds(i,1) = transDB(i).Rds.approx;
    Coss(i,1) = transDB(i).Coss.approx;
    Length(i,1) = transDB(i).Length.approx;
    Width(i,1) = transDB(i).Width.approx;
    Height(i,1) = transDB(i).Height.approx;
end


TransTable = table(Vds,Ids,Rds,Coss,Length,Width,Height,...
    'RowNames',NameT);



UsableTrans = NameT(TransTable.Vds >= 40 & TransTable.Ids >= 33);

TransTable(TransTable.Vds >= 40 & TransTable.Ids >= 33,:)



