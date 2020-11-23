% conn = sqlite('si_mos_clean.sql');
% sqlquery = 'SELECT * from si_mos_clean';
% results = fetch(conn,sqlquery)

%% Initial
% conn = sqlite('MOS.db','create');
% createTable = ['create table FETs (ron NUMERIC, vds NUMERIC, Id NUMERIC, Qg NUMERIC)'];
% exec(conn, createTable);

%% Insert FET
% conn = sqlite('MOS.db');
% insertFET = ['INSERT into FETs (ron, vds, Id, Qg) VALUES (0.16, 6, 1.8, 0)'];
% exec(conn, insertFET);
% close(conn);

%% Read Data
% conn = sqlite('MOS.db');
% getFETs = ['Select * from FETs'];
% results = fetch(conn,getFETs);
% close(conn);

%% Create and insert
% transistor(partno, ron, qg, Rjc, Id, CossVds, IsdVsd, RonT, Vbr)
EPC2034C = transistor('EPC2034C', 6e-3, 11.1e-9, 0.3, 48, [], [], [], 200);
Addplots(EPC2034C);
EPC2034C.datasheet();

%% Create nrdb

FETjson(1).PartNo = 'EPC2034C';
FETjson(1).PlotData = EPC2034C.PlotData;
FETjson(1).TableData = EPC2034C.TableData;

FETjson(2).PartNo = 'EPC2034C2';
FETjson(2).PlotData = EPC2034C.PlotData;
FETjson(2).TableData = EPC2034C.TableData;

nrdb = jsonencode(FETjson);

jsondecode(nrdb)
% FETs = transistorDB('transistors.db','transistors.json');
% FETs.addTransistor(EPC2034C)


function Addplots(transistor)
    transistor.PlotData(1).Type = 'CV';
    transistor.PlotData(1).Content = 'CrssVds';
    transistor.PlotData(1).Xlabel = 'VDS � Drain-to-Source Voltage [V]';
    transistor.PlotData(1).Ylabel = 'Capacitance [pF]';
    transistor.PlotData(1).TestCond = [];
    transistor.PlotData(1).Title = 'Capacitance (Linear Scale)';
    transistor.PlotData(1).Data = [0.11551847389541368, 25.732183976231294
0.47009622551925156, 25.70798888631885
0.8476374204010726, 25.33364214408111
1.2251866446148276, 24.933441980609615
1.6027499201594684, 24.485739379014948
1.9803192177030584, 24.023457759205332
2.3578905225796354, 23.562511514172215
2.735467849455164, 23.088671770956672
3.1130451763306874, 22.62436089745079
3.4906164812072635, 22.190259599191258
3.868175742085938, 21.80548908179793
4.2457289809656595, 21.447564081322902
4.623278205179415, 21.10875300147127
5.000821407394221, 20.79485398282984
5.378362602276042, 20.49204982812892
5.755895767825928, 20.219008519899038
6.133426926042833, 19.955864126763572
6.510960091592722, 19.68996709382987
6.888493257142608, 19.427612940907462
7.266026422692495, 19.168754461757782
7.643555573576418, 18.92521440924157
8.021088739126304, 18.673050016505368
8.398619897343208, 18.430025814317332
8.776153062893096, 18.18445944095511
9.153688235775965, 17.936537778810965
9.531221401325851, 17.697546766134508
9.908756574208724, 17.456263530612425
10.286291747091596, 17.218269880958562
10.663828927307451, 16.978194353661998
11.041366107523304, 16.741466216039605
11.418903287739157, 16.508038795205756
11.796440467955012, 16.277866069038065
12.173983670169818, 16.035805053671993
12.551522857718656, 15.807257462118496
12.929062045267491, 15.581967206347663
13.306607254815283, 15.345440232173706
13.684150457030087, 15.117244912954018
14.061693659244892, 14.892442986359512
14.4392368614597, 14.67098399073499
14.799596691135294, 14.513191991080276
];


    transistor.PlotData(2).Type = 'CV';
    transistor.PlotData(2).Content = 'CossVds';
    transistor.PlotData(2).Xlabel = 'VDS � Drain-to-Source Voltage [V]';
    transistor.PlotData(2).Ylabel = 'Capacitance [pF]';
    transistor.PlotData(2).TestCond = [];
    transistor.PlotData(2).Title = 'Capacitance (Linear Scale)';
transistor.PlotData(2).Data = [0.14023225575581233, 73.77495983048881
0.5319802940695668, 73.73503921813943
0.9094713056267963, 73.23340354619361
1.2869643245170141, 72.71236836645906
1.6644593507402128, 72.17239736733767
2.041962406295351, 71.5466079638582
2.4194594398515323, 70.9930212148981
2.796954466074733, 70.46581829939332
3.174453506963895, 69.89866458828648
3.55194652585411, 69.4013551433691
3.9294395447443233, 68.90758391603372
4.306932563634543, 68.41732573285202
4.684419560525812, 67.99451176421229
5.06190455008409, 67.59551102329827
5.43938552497641, 67.24102323581603
5.816870514534686, 66.84644406472147
6.194351489427007, 66.49588457195591
6.5718344716522985, 66.12641748343727
6.949311431878648, 65.8209149053699
7.326786384772014, 65.5373785104689
7.704269366997311, 65.17323711159473
8.081752349222613, 64.81111896969831
8.459227302115972, 64.53193246722279
8.836700247676356, 64.27410717778054
9.214185237234634, 63.896938260464786
9.591660190127996, 63.621689769227174
9.969133135688383, 63.36750118456197
10.346616117913678, 63.01541614471221
10.72409307814003, 62.724286323638495
11.101564016367426, 62.49328318233983
11.479046998592716, 62.14605550745473
11.856521951486084, 61.878349284269646
12.233992889713482, 61.650461588755896
12.611475871938778, 61.30791683285628
12.988952832165126, 61.02467562062717
13.366417748393575, 60.85717459366051
13.743900730618876, 60.51903753063281
14.121369661513285, 60.3150728288192
14.498842607073664, 60.07409521486478
14.790531843025915, 59.836194539448286
];


    transistor.PlotData(3).Type = 'CV';
    transistor.PlotData(3).Content = 'CissVds';
    transistor.PlotData(3).Xlabel = 'VDS � Drain-to-Source Voltage [V]';
    transistor.PlotData(3).Ylabel = 'Capacitance [pF]';
    transistor.PlotData(3).TestCond = [];
    transistor.PlotData(3).Title = 'Capacitance (Linear Scale)';
transistor.PlotData(3).Data = [0.15314035474737317, 91.66126436395534
0.5306153076407455, 91.26641563711546
0.9080822312021738, 90.98736112794163
1.285557184095532, 90.59541537038321
1.663040166320833, 90.09204550673681
2.0405171265471767, 89.67582226021716
2.4179940867735232, 89.26152195583927
2.795475061665844, 88.79341218129355
3.172960051224119, 88.27236075626857
3.550439018783453, 87.83698708033299
3.927907949677863, 87.54095387146076
4.3053889245701775, 87.08186718682461
4.682851833465645, 86.87008960547439
5.060330801024981, 86.44163215949287
5.437789695254475, 86.28552729876144
5.815268662813806, 85.8599530093887
6.19272554971031, 85.73178702915806
6.570202509936662, 85.33570807982903
6.947657389500194, 85.23505729248528
7.325134349726543, 84.84127322349201
7.70259123662305, 84.71462786013481
8.080060167517468, 84.42911780386532
8.457529098411882, 84.14456998982351
8.834985985308395, 84.01896461833385
9.212462945534739, 83.63079887044825
9.589915817765283, 83.55836573724119
9.96739478532461, 83.14624225334981
10.344851672221125, 83.02212711855901
10.722316588449573, 82.79424730745168
11.099789534009956, 82.46345834824862
11.477244413573477, 82.36619528924787
11.854721373799824, 81.9856652989037
12.232174246030375, 81.91465702568532
12.609651206256721, 81.5362131320762
12.987112107819202, 81.36344059340112
13.364568994715711, 81.24198670537666
13.742045954942059, 80.86665053368662
14.11949681983962, 80.82195998246931
14.496971772732989, 80.4738036677902
14.805807981586337, 80.27700807276624
];

    transistor.PlotData(4).Type = 'VQ';
    transistor.PlotData(4).Content = 'VdsQg';
    transistor.PlotData(4).Xlabel = 'QG � Gate Charge [nC]';
    transistor.PlotData(4).Ylabel = 'VGS � Gate-to-Source Voltage [V]';
    transistor.PlotData(4).TestCond = 'ID = 1.5 A, Vds = 6 V';
    transistor.PlotData(4).Title = 'Gate Charge';
transistor.PlotData(4).Data = [0.009457930340156765, 0.14201979175952986
0.02683637166777404, 0.31870239563868846
0.04147244202431391, 0.477885001473131
0.05610836951047362, 0.6361353780769144
0.07074443986701351, 0.7953179839113572
0.085380224482793, 0.9526361312844814
0.1000157233578122, 1.1080898201962874
0.11465165084397187, 1.2663401968000707
0.1302015042776414, 1.4291338514935035
0.1466655059016343, 1.5979209186353898
0.16312849155847928, 1.7600788001370318
0.1805038665125228, 1.9167533164474557
0.19970620669934772, 2.0767289178329658
0.21982194060543792, 2.2377750080124996
0.2399235433321068, 2.3066151524686465
0.2600145476742099, 2.3063008376322562
0.28010555201631293, 2.305986522795865
0.300196556358416, 2.305672207959474
0.320287560700519, 2.3053578931230834
0.3403785650426221, 2.3050435782866927
0.36047019281911136, 2.308797172820451
0.38057699083233226, 2.4115365620278433
0.4006906466238017, 2.559022954306878
0.42080419850954004, 2.7058313616908887
0.4409175425838163, 2.851283799284848
0.4610304710351685, 2.9940242972987092
0.48114381510944476, 3.139476734892669
0.5012571591837209, 3.2849291724866294
0.5213700876350731, 3.42766967050049
0.5414832238978874, 3.5717661383044
0.5615964640664326, 3.716540591003336
0.5817094964235159, 3.859959073912221
0.6018225287805992, 4.003377556821107
0.6219357689491443, 4.148152009520041
0.6420485934947653, 4.290214522638878
0.6621616258518485, 4.433633005547764
0.6822748660203938, 4.578407458246698
0.7023874827545529, 4.719114001575484
0.722501626902958, 4.869786922861136
0.7398688645951262, 4.973365805434236
];

    transistor.PlotData(5).Type = 'VQtest';
    transistor.PlotData(5).Content = 'VdsQg';
    transistor.PlotData(5).Xlabel = 'QG � Gate Charge [nC]';
    transistor.PlotData(5).Ylabel = 'VGS � Gate-to-Source Voltage [V]';
    transistor.PlotData(5).TestCond = 'ID = 1.5 A, Vds = 6 V';
    transistor.PlotData(5).Title = 'Gate Charge';
transistor.PlotData(5).Data = [0.009457930340156765, 0.14201979175952986
0.02683637166777404, 0.31870239563868846
0.04147244202431391, 0.477885001473131
0.05610836951047362, 0.6361353780769144
0.07074443986701351, 0.7953179839113572
0.085380224482793, 0.9526361312844814
0.1000157233578122, 1.1080898201962874
0.11465165084397187, 1.2663401968000707
0.1302015042776414, 1.4291338514935035
0.1466655059016343, 1.5979209186353898
0.16312849155847928, 1.7600788001370318
0.1805038665125228, 1.9167533164474557
0.19970620669934772, 2.0767289178329658
0.21982194060543792, 2.2377750080124996
0.2399235433321068, 2.3066151524686465
0.2600145476742099, 2.3063008376322562
0.28010555201631293, 2.305986522795865
0.300196556358416, 2.305672207959474
0.320287560700519, 2.3053578931230834
0.3403785650426221, 2.3050435782866927
0.36047019281911136, 2.308797172820451
0.38057699083233226, 2.4115365620278433
0.4006906466238017, 2.559022954306878
0.42080419850954004, 2.7058313616908887
0.4409175425838163, 2.851283799284848
0.4610304710351685, 2.9940242972987092
0.48114381510944476, 3.139476734892669
0.5012571591837209, 3.2849291724866294
0.5213700876350731, 3.42766967050049
0.5414832238978874, 3.5717661383044
0.5615964640664326, 3.716540591003336
0.5817094964235159, 3.859959073912221
0.6018225287805992, 4.003377556821107
0.6219357689491443, 4.148152009520041
0.6420485934947653, 4.290214522638878
0.6621616258518485, 4.433633005547764
0.6822748660203938, 4.578407458246698
0.7023874827545529, 4.719114001575484
0.722501626902958, 4.869786922861136
0.7398688645951262, 4.973365805434236
];
end
