%% Inductors Database

indDB = inductorDB;


%% XGL3515900ME
XGL3515900ME = inductor('XGL3515-900ME','Power','Composite');
XGL3515900ME.manufacturer = 'Coilcraft';
XGL3515900ME.material = 'Composite';
XGL3515900ME.addParameter('L', 0.09e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3515900ME.addParameter('Rdc', 2.2e-3 ,'typ')
XGL3515900ME.addParameter('Rdc', 2.6e-3 ,'max')
XGL3515900ME.addParameter('SRF', 350000000 ,'typ' )
XGL3515900ME.addParameter('Isat', 7.6 ,'typ', {'10% drop in L'} )
XGL3515900ME.addParameter('Isat', 12.6 ,'typ', {'20% drop in L'} )
XGL3515900ME.addParameter('Isat', 17.8 ,'typ', {'30% drop in L'} )
XGL3515900ME.addParameter('Irms', 19.5 ,'typ', {'20°C rise in L'} )
XGL3515900ME.addParameter('Irms', 26.6 ,'typ', {'40°C rise in L'} )
XGL3515900ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3515900ME.addParameter('Fspice', 650e6 ,'max')
XGL3515900ME.addParameter('R1spice', 5 ,'typ')
XGL3515900ME.addParameter('R2spice', 0.0026 ,'typ')
XGL3515900ME.addParameter('Cspice', 2.4e-12 ,'typ')
XGL3515900ME.addParameter('K1spice', 1e-6 ,'typ')
XGL3515900ME.addParameter('K2spice', 0.0013 ,'typ')
XGL3515900ME.addParameter('K3spice', 0.09 ,'typ')
XGL3515900ME.addParameter('K4spice', 1e-6 ,'typ')
XGL3515900ME.addParameter('K5spice', 1e-6 ,'typ')
XGL3515900ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3515900ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3515900ME.addParameter('Height', 1.5e-3 ,'typ')
XGL3515900ME_DP = digitizedPlot([]);

XGL3515900ME_DP.plotData = {[0, 0.09000000000000002
1.0517241379310331, 0.09000000000000002
2.707828518173345, 0.08837837837837839
4.363932898415657, 0.08675675675675679
7.539142590866728, 0.0818918918918919
9.74883504193849, 0.07702702702702703
12.510251630941287, 0.07216216216216216
15.409599254426842, 0.06729729729729733
17.34249767008388, 0.06405405405405405
20.794501397949674, 0.05756756756756756
23.554986020503257, 0.05432432432432434]};

XGL3515900ME_DP.logAxes = [0 0];
XGL3515900ME_DP.normAxes = [0 0];
XGL3515900ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3515900ME_DP.dataLabels = {'Y1'};
XGL3515900ME_DP.SIUnits = { '' , char(181)};


XGL3515900ME_graph = componentPlotData(inductor, XGL3515900ME_DP);
XGL3515900ME.addGraph(XGL3515900ME_graph);

XGL3515900ME_DP = digitizedPlot([]);


XGL3515900ME_DP.plotData = {[0.1, 0.09
21.153269486715, 0.08916790781594854
32.94528756210256, 0.08558215188848489
42.97770509076642, 0.08359894530188042
53.476744133487635, 0.08117852066213044
67.33148480489835, 0.0788287711850899
86.80206045830667, 0.07564962010327668
126.68488228157155, 0.07303438710811522
208.08515701516816, 0.07177500524117415]};

XGL3515900ME_DP.logAxes = [1 1];
XGL3515900ME_DP.normAxes = [0 0];
XGL3515900ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3515900ME_DP.dataLabels = {'Y1'};
XGL3515900ME_DP.SIUnits = { 'M' , char(181)};


XGL3515900ME_graph = componentPlotData(inductor, XGL3515900ME_DP);
XGL3515900ME.addGraph(XGL3515900ME_graph);

indDB.add(XGL3515900ME)


%% XGL3515181ME

XGL3515181ME = inductor('XGL3515-181ME','Power','Composite');
XGL3515181ME.manufacturer = 'Coilcraft';
XGL3515181ME.material = 'Composite';
XGL3515181ME.addParameter('L', 0.18e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3515181ME.addParameter('Rdc', 3.5e-3 ,'typ')
XGL3515181ME.addParameter('Rdc', 4.1e-3 ,'max')
XGL3515181ME.addParameter('SRF', 205000000 ,'typ' )
XGL3515181ME.addParameter('Isat', 5.4 ,'typ', {'10% drop in L'} )
XGL3515181ME.addParameter('Isat', 8.8 ,'typ', {'20% drop in L'} )
XGL3515181ME.addParameter('Isat', 12.4 ,'typ', {'30% drop in L'} )
XGL3515181ME.addParameter('Irms', 15.1 ,'typ', {'20°C rise in L'} )
XGL3515181ME.addParameter('Irms', 20.4 ,'typ', {'40°C rise in L'} )
XGL3515181ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3515181ME.addParameter('Fspice', 500e6 ,'max')
XGL3515181ME.addParameter('R1spice', 6 ,'typ')
XGL3515181ME.addParameter('R2spice', 4.10e-3 ,'typ')
XGL3515181ME.addParameter('Cspice', 3.1e-12 ,'typ')
XGL3515181ME.addParameter('K1spice', 1e-6 ,'typ')
XGL3515181ME.addParameter('K2spice', 0.0270 ,'typ')
XGL3515181ME.addParameter('K3spice', 0.18 ,'typ')
XGL3515181ME.addParameter('K4spice', 1e-6 ,'typ')
XGL3515181ME.addParameter('K5spice', 1e-6 ,'typ')
XGL3515181ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3515181ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3515181ME.addParameter('Height', 1.5e-3 ,'typ')

XGL3515181ME_DP = digitizedPlot([]);

XGL3515181ME_DP.plotData = {[0, 0.17852696584724065
0.7845958527295805, 0.17819084860242265
1.531104528142192, 0.17692756103794885
2.5823106220905623, 0.1744161746898346
3.6639864578925088, 0.17035671113662357
4.791366906474818, 0.16567686444903482
5.781633516716037, 0.16069232558042093
7.12230215827338, 0.15322049291260054
8.363944138806607, 0.14606147056728525
9.498942022852304, 0.13983433230119457
10.732966567922132, 0.13298497777200488
12.24121878967414, 0.12581678978308325
13.581887431231484, 0.11927317493704509
15.151079136690644, 0.11364992164684093
]};

XGL3515181ME_DP.logAxes = [0 0];
XGL3515181ME_DP.normAxes = [0 0];
XGL3515181ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3515181ME_DP.dataLabels = {'Y1'};
XGL3515181ME_DP.SIUnits = { '' , char(181)};


XGL3515181ME_graph = componentPlotData(inductor, XGL3515181ME_DP);
XGL3515181ME.addGraph(XGL3515181ME_graph);

XGL3515181ME_DP = digitizedPlot([]);


XGL3515181ME_DP.plotData = {[0.10, 0.18
11.272112480396759, 0.18
21.718560876054816, 0.17998470721797707
32.56253253746654, 0.17155072591380072
47.889791876375185, 0.16749918113614023
68.42440650228761, 0.16592307004536755
95.89830426477921, 0.16675428861834854
119.71592415723538, 0.16838928580233806
150.1774399723319, 0.17503063317997464]};

XGL3515181ME_DP.logAxes = [1 1];
XGL3515181ME_DP.normAxes = [0 0];
XGL3515181ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3515181ME_DP.dataLabels = {'Y1'};
XGL3515181ME_DP.SIUnits = { 'M' , char(181)};


XGL3515181ME_graph = componentPlotData(inductor, XGL3515181ME_DP);
XGL3515181ME.addGraph(XGL3515181ME_graph);

indDB.add(XGL3515181ME)

%% XGL3515331ME

XGL3515331ME = inductor('XGL3515-331ME','Power','Composite');
XGL3515331ME.manufacturer = 'Coilcraft';
XGL3515331ME.material = 'Composite';
XGL3515331ME.addParameter('L', 0.33e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3515331ME.addParameter('Rdc', 6.4e-3 ,'typ')
XGL3515331ME.addParameter('Rdc', 7.4e-3 ,'max')
XGL3515331ME.addParameter('SRF', 145000000 ,'typ' )
XGL3515331ME.addParameter('Isat', 3.9 ,'typ', {'10% drop in L'} )
XGL3515331ME.addParameter('Isat', 6.3 ,'typ', {'20% drop in L'} )
XGL3515331ME.addParameter('Isat', 8.7 ,'typ', {'30% drop in L'} )
XGL3515331ME.addParameter('Irms', 11.0 ,'typ', {'20°C rise in L'} )
XGL3515331ME.addParameter('Irms', 14.8 ,'typ', {'40°C rise in L'} )
XGL3515331ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3515331ME.addParameter('Fspice', 500e6 ,'max')
XGL3515331ME.addParameter('R1spice', 8 ,'typ')
XGL3515331ME.addParameter('R2spice', 7.40e-3 ,'typ')
XGL3515331ME.addParameter('Cspice', 3.3e-12 ,'typ')
XGL3515331ME.addParameter('K1spice', 1e-6 ,'typ')
XGL3515331ME.addParameter('K2spice', 0.0460 ,'typ')
XGL3515331ME.addParameter('K3spice', 0.33 ,'typ')
XGL3515331ME.addParameter('K4spice', 1e-6 ,'typ')
XGL3515331ME.addParameter('K5spice', 1e-6 ,'typ')
XGL3515331ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3515331ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3515331ME.addParameter('Height', 1.5e-3 ,'typ')

XGL3515331ME_DP = digitizedPlot([]);

XGL3515331ME_DP.plotData = {[0, 0.3253333333333333
0.5889508362899144, 0.32533333333333336
1.37658388241257, 0.32059259259259265
2.2635580334515963, 0.31288888888888894
3.1718195641155598, 0.30281481481481487
4.647744551444501, 0.2838518518518519
5.392802838317283, 0.27377777777777784
6.187531677648252, 0.2631111111111112
7.337050177394831, 0.24592592592592585
8.167257982767362, 0.2346666666666668
8.898124683223516, 0.2245925925925927
9.664470349721235, 0.2145185185185185
10.892042574759252, 0.20029629629629642
11.53775975671566, 0.1949629629629631]};

XGL3515331ME_DP.logAxes = [0 0];
XGL3515331ME_DP.normAxes = [0 0];
XGL3515331ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3515331ME_DP.dataLabels = {'Y1'};
XGL3515331ME_DP.SIUnits = { '' , char(181)};


XGL3515331ME_graph = componentPlotData(inductor, XGL3515331ME_DP);
XGL3515331ME.addGraph(XGL3515331ME_graph);

XGL3515331ME_DP = digitizedPlot([]);


XGL3515331ME_DP.plotData = {[0.10, 0.33
18.718287878265404, 0.3240485246911819
35.205609902204166, 0.32415567929687505
55.39575690739808, 0.32737362678297677
76.89629973365635, 0.3386669407183749
98.81782116155613, 0.35543886814666265]};

XGL3515331ME_DP.logAxes = [1 1];
XGL3515331ME_DP.normAxes = [0 0];
XGL3515331ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3515331ME_DP.dataLabels = {'Y1'};
XGL3515331ME_DP.SIUnits = { 'M' , char(181)};


XGL3515331ME_graph = componentPlotData(inductor, XGL3515331ME_DP);
XGL3515331ME.addGraph(XGL3515331ME_graph);

indDB.add(XGL3515331ME)


%% XGL3515561ME

XGL3515561ME = inductor('XGL3515-561ME','Power','Composite');
XGL3515561ME.manufacturer = 'Coilcraft';
XGL3515561ME.material = 'Composite';
XGL3515561ME.addParameter('L', 0.56e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3515561ME.addParameter('Rdc', 10.2e-3 ,'typ')
XGL3515561ME.addParameter('Rdc', 11.8e-3 ,'max')
XGL3515561ME.addParameter('SRF', 100000000 ,'typ' )
XGL3515561ME.addParameter('Isat', 3.0 ,'typ', {'10% drop in L'} )
XGL3515561ME.addParameter('Isat', 4.8 ,'typ', {'20% drop in L'} )
XGL3515561ME.addParameter('Isat', 6.6 ,'typ', {'30% drop in L'} )
XGL3515561ME.addParameter('Irms', 8.1 ,'typ', {'20°C rise in L'} )
XGL3515561ME.addParameter('Irms', 10.9 ,'typ', {'40°C rise in L'} )
XGL3515561ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3515561ME.addParameter('Fspice', 400e6 ,'max')
XGL3515561ME.addParameter('R1spice', 5 ,'typ')
XGL3515561ME.addParameter('R2spice', 1.18e-2 ,'typ')
XGL3515561ME.addParameter('Cspice', 4.4e-12 ,'typ')
XGL3515561ME.addParameter('K1spice', 1e-6 ,'typ')
XGL3515561ME.addParameter('K2spice', 0.0777 ,'typ')
XGL3515561ME.addParameter('K3spice', 0.56 ,'typ')
XGL3515561ME.addParameter('K4spice', 1e-6 ,'typ')
XGL3515561ME.addParameter('K5spice', 1e-6 ,'typ')
XGL3515561ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3515561ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3515561ME.addParameter('Height', 1.5e-3 ,'typ')

XGL3515561ME_DP = digitizedPlot([]);

XGL3515561ME_DP.plotData = {[0, 0.541553637484587
0.405405405405405, 0.5425400739827375
0.8868243243243239, 0.5385943279901357
1.3935810810810811, 0.5297163995067818
1.9104729729729724, 0.5198520345252775
2.7618243243243237, 0.4971639950678175
3.390202702702703, 0.47842170160295927
4.109797297297299, 0.45770653514180026
4.738175675675675, 0.4379778051787917
5.336148648648648, 0.41923551171393336
6.390202702702707, 0.3847102342786684
7.201013513513516, 0.36202219482120856
7.859797297297297, 0.3422934648581999
8.407094594594597, 0.32848335388409355
8.827702702702704, 0.3205918618988902]};

XGL3515561ME_DP.logAxes = [0 0];
XGL3515561ME_DP.normAxes = [0 0];
XGL3515561ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3515561ME_DP.dataLabels = {'Y1'};
XGL3515561ME_DP.SIUnits = { '' , char(181)};


XGL3515561ME_graph = componentPlotData(inductor, XGL3515561ME_DP);
XGL3515561ME.addGraph(XGL3515561ME_graph);

XGL3515561ME_DP = digitizedPlot([]);


XGL3515561ME_DP.plotData = {[0.1, 0.56
5.09501095775818, 0.56
9.49085439615902, 0.555813135206041
14.507702429916147, 0.5559365926865236
25.382628854611433, 0.5614866430456364
45.71733798582625, 0.6066934937480514
67.24971884921256, 0.6812418245889834]};

XGL3515561ME_DP.logAxes = [1 1];
XGL3515561ME_DP.normAxes = [0 0];
XGL3515561ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3515561ME_DP.dataLabels = {'Y1'};
XGL3515561ME_DP.SIUnits = { 'M' , char(181)};


XGL3515561ME_graph = componentPlotData(inductor, XGL3515561ME_DP);
XGL3515561ME.addGraph(XGL3515561ME_graph);

indDB.add(XGL3515561ME)


%% XGL4015101ME

XGL4015101ME = inductor('XGL4015-101ME','Power','Composite');
XGL4015101ME.manufacturer = 'Coilcraft';
XGL4015101ME.material = 'Composite';
XGL4015101ME.addParameter('L', 0.1e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4015101ME.addParameter('Rdc', 1.8e-3 ,'typ')
XGL4015101ME.addParameter('Rdc', 2.2e-3 ,'max')
XGL4015101ME.addParameter('SRF', 260000000 ,'typ' )
XGL4015101ME.addParameter('Isat', 8.9 ,'typ', {'10% drop in L'} )
XGL4015101ME.addParameter('Isat', 16.0 ,'typ', {'20% drop in L'} )
XGL4015101ME.addParameter('Isat', 24.5 ,'typ', {'30% drop in L'} )
XGL4015101ME.addParameter('Irms', 19.0 ,'typ', {'20°C rise in L'} )
XGL4015101ME.addParameter('Irms', 25.9 ,'typ', {'40°C rise in L'} )
XGL4015101ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4015101ME.addParameter('Fspice', 500e6 ,'max')
XGL4015101ME.addParameter('R1spice', 6.8 ,'typ')
XGL4015101ME.addParameter('R2spice', 2.2e-3 ,'typ')
XGL4015101ME.addParameter('Cspice', 3.5e-12 ,'typ')
XGL4015101ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4015101ME.addParameter('K2spice', 0.016 ,'typ')
XGL4015101ME.addParameter('K3spice', 0.1 ,'typ')
XGL4015101ME.addParameter('K4spice', 8.35e-4 ,'typ')
XGL4015101ME.addParameter('K5spice', 1.04e-4 ,'typ')
XGL4015101ME.addParameter('Length', 4.0e-3 ,'typ')
XGL4015101ME.addParameter('Width', 4.0e-3 ,'typ')
XGL4015101ME.addParameter('Height', 1.5e-3 ,'typ')

XGL4015101ME_DP = digitizedPlot([]);

XGL4015101ME_DP.plotData = {[0, 0.09465484680428346
2.9424877396344185, 0.09310146996508331
5.3945608559964295, 0.09115538941664668
9.40704413731609, 0.08608633785190609
11.814534106107892, 0.08296668005517385
14.890771288452962, 0.078285100800346
19.037004012483283, 0.07282543808609841
21.711992866696388, 0.06892455799035217
26.482389656709756, 0.06307620230634962 
29.201961658493087, 0.060348899458899
32.902362906821224, 0.05723429868151303
35.131520285332144, 0.055287346233189094
38.742755238519834, 0.053736933853605756
39.81275078020507, 0.05334998468362542]};

XGL4015101ME_DP.logAxes = [0 0];
XGL4015101ME_DP.normAxes = [0 0];
XGL4015101ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4015101ME_DP.dataLabels = {'Y1'};
XGL4015101ME_DP.SIUnits = { '' , char(181)};


XGL4015101ME_graph = componentPlotData(inductor, XGL4015101ME_DP);
XGL4015101ME.addGraph(XGL4015101ME_graph);

XGL4015101ME_DP = digitizedPlot([]);


XGL4015101ME_DP.plotData = {[0.1, 0.09226718495210182
9.188383380501373, 0.09295639980146427
20.016269587420464, 0.09155435431547598
37.611462080852526, 0.08973735411829478
55.78643359458444, 0.08796754488088265
113.42534180315404, 0.0866439431405688
174.14415619614482, 0.0874818990679973
]};

XGL4015101ME_DP.logAxes = [1 1];
XGL4015101ME_DP.normAxes = [0 0];
XGL4015101ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4015101ME_DP.dataLabels = {'Y1'};
XGL4015101ME_DP.SIUnits = { 'M' , char(181)};


XGL4015101ME_graph = componentPlotData(inductor, XGL4015101ME_DP);
XGL4015101ME.addGraph(XGL4015101ME_graph);

indDB.add(XGL4015101ME)


%% XGL4015221ME

XGL4015221ME = inductor('XGL4015-221ME','Power','Composite');
XGL4015221ME.manufacturer = 'Coilcraft';
XGL4015221ME.material = 'Composite';
XGL4015221ME.addParameter('L', 0.22e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4015221ME.addParameter('Rdc', 3.4e-3 ,'typ')
XGL4015221ME.addParameter('Rdc', 4.1e-3 ,'max')
XGL4015221ME.addParameter('SRF', 150000000 ,'typ' )
XGL4015221ME.addParameter('Isat', 5.9 ,'typ', {'10% drop in L'} )
XGL4015221ME.addParameter('Isat', 10.1 ,'typ', {'20% drop in L'} )
XGL4015221ME.addParameter('Isat', 15.2 ,'typ', {'30% drop in L'} )
XGL4015221ME.addParameter('Irms', 13.0 ,'typ', {'20°C rise in L'} )
XGL4015221ME.addParameter('Irms', 18.0 ,'typ', {'40°C rise in L'} )
XGL4015221ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4015221ME.addParameter('Fspice', 500e6 ,'max')
XGL4015221ME.addParameter('R1spice', 7.5 ,'typ')
XGL4015221ME.addParameter('R2spice', 4.1e-3 ,'typ')
XGL4015221ME.addParameter('Cspice', 4.8e-12 ,'typ')
XGL4015221ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4015221ME.addParameter('K2spice', 0.035 ,'typ')
XGL4015221ME.addParameter('K3spice', 0.22 ,'typ')
XGL4015221ME.addParameter('K4spice', 1.61e-4 ,'typ')
XGL4015221ME.addParameter('K5spice', 1.04e-4 ,'typ')
XGL4015221ME.addParameter('Length', 4.0e-3 ,'typ')
XGL4015221ME.addParameter('Width', 4.0e-3 ,'typ')
XGL4015221ME.addParameter('Height', 1.5e-3 ,'typ')

XGL4015221ME_DP = digitizedPlot([]);

XGL4015221ME_DP.plotData = {[0, 0.2106631989596879
1.9881055793757318, 0.20754226267880363
4.164776429805335, 0.199739921976593
6.7336649008982015, 0.18803641092327697
9.195458325136986, 0.17633289986996098
12.620457266140132, 0.15916775032509747
14.868339127351774, 0.1498049414824447
17.61590503032213, 0.13966189856957084
20.970250021905276, 0.12873862158647598
23.18289775196727, 0.12327698309492846
25.14565703936973, 0.11781534460338097
27.358304769431722, 0.11235370611183353
28.35758000938881, 0.11001300390117034
]};

XGL4015221ME_DP.logAxes = [0 0];
XGL4015221ME_DP.normAxes = [0 0];
XGL4015221ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4015221ME_DP.dataLabels = {'Y1'};
XGL4015221ME_DP.SIUnits = { '' , char(181)};


XGL4015221ME_graph = componentPlotData(inductor, XGL4015221ME_DP);
XGL4015221ME.addGraph(XGL4015221ME_graph);

XGL4015221ME_DP = digitizedPlot([]);


XGL4015221ME_DP.plotData = {[0.10, 0.21426586069182588
10.924661959123759, 0.2169125763344823
36.900152358087, 0.2135908019644648
58.93224644220127, 0.21886321734425257
77.28089016039117, 0.22539724687434878
99.85730737224378, 0.23909289267193135
]};

XGL4015221ME_DP.logAxes = [1 1];
XGL4015221ME_DP.normAxes = [0 0];
XGL4015221ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4015221ME_DP.dataLabels = {'Y1'};
XGL4015221ME_DP.SIUnits = { 'M' , char(181)};


XGL4015221ME_graph = componentPlotData(inductor, XGL4015221ME_DP);
XGL4015221ME.addGraph(XGL4015221ME_graph);

indDB.add(XGL4015221ME)


%% XGL4015471ME

XGL4015471ME = inductor('XGL4015-471ME','Power','Composite');
XGL4015471ME.manufacturer = 'Coilcraft';
XGL4015471ME.material = 'Composite';
XGL4015471ME.addParameter('L', 0.47e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4015471ME.addParameter('Rdc', 6.2e-3 ,'typ')
XGL4015471ME.addParameter('Rdc', 7.5e-3 ,'max')
XGL4015471ME.addParameter('SRF', 95000000 ,'typ' )
XGL4015471ME.addParameter('Isat', 4.4 ,'typ', {'10% drop in L'} )
XGL4015471ME.addParameter('Isat', 7.3 ,'typ', {'20% drop in L'} )
XGL4015471ME.addParameter('Isat', 10.5 ,'typ', {'30% drop in L'} )
XGL4015471ME.addParameter('Irms', 9.6 ,'typ', {'20°C rise in L'} )
XGL4015471ME.addParameter('Irms', 13.0 ,'typ', {'40°C rise in L'} )
XGL4015471ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4015471ME.addParameter('Fspice', 300e6 ,'max')
XGL4015471ME.addParameter('R1spice', 3.0 ,'typ')
XGL4015471ME.addParameter('R2spice', 7.5e-3 ,'typ')
XGL4015471ME.addParameter('Cspice', 5.2e-12 ,'typ')
XGL4015471ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4015471ME.addParameter('K2spice', 0.070 ,'typ')
XGL4015471ME.addParameter('K3spice', 0.47 ,'typ')
XGL4015471ME.addParameter('K4spice', 9.22e-5 ,'typ')
XGL4015471ME.addParameter('K5spice', 1.04e-4 ,'typ')
XGL4015471ME.addParameter('Length', 4.0e-3 ,'typ')
XGL4015471ME.addParameter('Width', 4.0e-3 ,'typ')
XGL4015471ME.addParameter('Height', 1.5e-3 ,'typ')

XGL4015471ME_DP = digitizedPlot([]);

XGL4015471ME_DP.plotData = {[0, 0.4703104818519697
1.230668667961964, 0.4640144270051577
3.3026955261047486, 0.4420584884560116
5.178113149173069, 0.41386022995511035
6.857014584604254, 0.3856696478677888
8.839622855474648, 0.357467202232208
10.822231126345045, 0.3292647565966272
13.251420572672082, 0.29791986638388
16.180786518356804, 0.2634304380266265
16.91313963570765, 0.25558933048298
]};

XGL4015471ME_DP.logAxes = [0 0];
XGL4015471ME_DP.normAxes = [0 0];
XGL4015471ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4015471ME_DP.dataLabels = {'Y1'};
XGL4015471ME_DP.SIUnits = { '' , char(181)};


XGL4015471ME_graph = componentPlotData(inductor, XGL4015471ME_DP);
XGL4015471ME.addGraph(XGL4015471ME_graph);

XGL4015471ME_DP = digitizedPlot([]);


XGL4015471ME_DP.plotData = {[0.1, 0.47
10.560060430288459, 0.4701552212123313
24.16735047207813, 0.4888406553636805
41.559793920873744, 0.528783836248639
63.498811446196115, 0.5950245683555894
]};

XGL4015471ME_DP.logAxes = [1 1];
XGL4015471ME_DP.normAxes = [0 0];
XGL4015471ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4015471ME_DP.dataLabels = {'Y1'};
XGL4015471ME_DP.SIUnits = { 'M' , char(181)};


XGL4015471ME_graph = componentPlotData(inductor, XGL4015471ME_DP);
XGL4015471ME.addGraph(XGL4015471ME_graph);

indDB.add(XGL4015471ME)

%% XGL4015681ME

XGL4015681ME = inductor('XGL4015-681ME','Power','Composite');
XGL4015681ME.manufacturer = 'Coilcraft';
XGL4015681ME.material = 'Composite';
XGL4015681ME.addParameter('L', 0.68e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4015681ME.addParameter('Rdc', 8.4e-3 ,'typ')
XGL4015681ME.addParameter('Rdc', 10.1e-3 ,'max')
XGL4015681ME.addParameter('SRF', 78000000 ,'typ' )
XGL4015681ME.addParameter('Isat', 3.6 ,'typ', {'10% drop in L'} )
XGL4015681ME.addParameter('Isat', 6.2 ,'typ', {'20% drop in L'} )
XGL4015681ME.addParameter('Isat', 9.0 ,'typ', {'30% drop in L'} )
XGL4015681ME.addParameter('Irms', 8.1 ,'typ', {'20°C rise in L'} )
XGL4015681ME.addParameter('Irms', 11.0 ,'typ', {'40°C rise in L'} )
XGL4015681ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4015681ME.addParameter('Fspice', 250e6 ,'max')
XGL4015681ME.addParameter('R1spice', 5.0 ,'typ')
XGL4015681ME.addParameter('R2spice', 1.01e-2 ,'typ')
XGL4015681ME.addParameter('Cspice', 6.0e-12 ,'typ')
XGL4015681ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4015681ME.addParameter('K2spice', 0.104 ,'typ')
XGL4015681ME.addParameter('K3spice', 0.68 ,'typ')
XGL4015681ME.addParameter('K4spice', 3.45e-3 ,'typ')
XGL4015681ME.addParameter('K5spice', 1.35e-3 ,'typ')
XGL4015681ME.addParameter('Length', 4.0e-3 ,'typ')
XGL4015681ME.addParameter('Width', 4.0e-3 ,'typ')
XGL4015681ME.addParameter('Height', 1.5e-3 ,'typ')

XGL4015681ME_DP = digitizedPlot([]);

XGL4015681ME_DP.plotData = {[0, 0.625
0.8352399553571445, 0.61875
2.01339285714286, 0.6000000000000001
3.155691964285716, 0.5750000000000001
4.886951264880953, 0.5354166666666668
6.67168898809524, 0.49166666666666664
8.42080543154762, 0.45208333333333334
10.27701822916667, 0.4104166666666666
11.75851004464286, 0.38125000000000003
13.73990885416667, 0.34791666666666654
16.203497023809522, 0.3166666666666667
]};

XGL4015681ME_DP.logAxes = [0 0];
XGL4015681ME_DP.normAxes = [0 0];
XGL4015681ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4015681ME_DP.dataLabels = {'Y1'};
XGL4015681ME_DP.SIUnits = { '' , char(181)};


XGL4015681ME_graph = componentPlotData(inductor, XGL4015681ME_DP);
XGL4015681ME.addGraph(XGL4015681ME_graph);

XGL4015681ME_DP = digitizedPlot([]);


XGL4015681ME_DP.plotData = {[0.10, 0.6149944704253409
12.863680316788962, 0.6318055490982214
21.90367613634016, 0.6570193853246803
34.469989243293334, 0.7177765143946596
48.19657812396218, 0.8238075668765414
60.46587352129199, 0.9549202341490749
]};

XGL4015681ME_DP.logAxes = [1 1];
XGL4015681ME_DP.normAxes = [0 0];
XGL4015681ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4015681ME_DP.dataLabels = {'Y1'};
XGL4015681ME_DP.SIUnits = { 'M' , char(181)};


XGL4015681ME_graph = componentPlotData(inductor, XGL4015681ME_DP);
XGL4015681ME.addGraph(XGL4015681ME_graph);

indDB.add(XGL4015681ME)


%% XGL3520101ME

XGL3520101ME = inductor('XGL3520-101ME','Power','Composite');
XGL3520101ME.manufacturer = 'Coilcraft';
XGL3520101ME.material = 'Composite';
XGL3520101ME.addParameter('L', 0.1e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3520101ME.addParameter('Rdc', 2.3e-3 ,'typ')
XGL3520101ME.addParameter('Rdc', 2.7e-3 ,'max')
XGL3520101ME.addParameter('SRF', 360000000 ,'typ' )
XGL3520101ME.addParameter('Isat', 7.5 ,'typ', {'10% drop in L'} )
XGL3520101ME.addParameter('Isat', 11.7 ,'typ', {'20% drop in L'} )
XGL3520101ME.addParameter('Isat', 15.9 ,'typ', {'30% drop in L'} )
XGL3520101ME.addParameter('Irms', 19.2 ,'typ', {'20°C rise in L'} )
XGL3520101ME.addParameter('Irms', 26.2 ,'typ', {'40°C rise in L'} )
XGL3520101ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3520101ME.addParameter('Fspice', 500e6 ,'max')
XGL3520101ME.addParameter('R1spice', 10 ,'typ')
XGL3520101ME.addParameter('R2spice', 2.70e-3 ,'typ')
XGL3520101ME.addParameter('Cspice', 2.4e-12 ,'typ')
XGL3520101ME.addParameter('K1spice', 2e-5 ,'typ')
XGL3520101ME.addParameter('K2spice', 0.014 ,'typ')
XGL3520101ME.addParameter('K3spice', 0.1 ,'typ')
XGL3520101ME.addParameter('K4spice', 5.60e-5 ,'typ')
XGL3520101ME.addParameter('K5spice', 3.70e-3 ,'typ')
XGL3520101ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3520101ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3520101ME.addParameter('Height', 2e-3 ,'typ')

XGL3520101ME_DP = digitizedPlot([]);

XGL3520101ME_DP.plotData = {[0 0.1002055978444066
4.228246237306459, 0.0963185149462
5.820504730274115, 0.09338256105370511
7.727707554567237, 0.08979366187110854
9.70480617124146, 0.08522005495354533
13.728899196971017, 0.07508813338345226
16.58075533454031, 0.06789768549696444
18.645430324496164, 0.06398231366931706
19.695232783534355, 0.06169631518007238

]};

XGL3520101ME_DP.logAxes = [0 0];
XGL3520101ME_DP.normAxes = [0 0];
XGL3520101ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3520101ME_DP.dataLabels = {'Y1'};
XGL3520101ME_DP.SIUnits = { '' , char(181)};


XGL3520101ME_graph = componentPlotData(inductor, XGL3520101ME_DP);
XGL3520101ME.addGraph(XGL3520101ME_graph);

XGL3520101ME_DP = digitizedPlot([]);


XGL3520101ME_DP.plotData = {[0.10, 0.10595757429370295
7.383026841154094, 0.10015485033876644
16.495644939388846, 0.09948566658444041
45.755326136647334, 0.09080633188403736
120.22258287845663, 0.07974285395467641
211.35690863302392, 0.07560473850567233

]};

XGL3520101ME_DP.logAxes = [1 1];
XGL3520101ME_DP.normAxes = [0 0];
XGL3520101ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3520101ME_DP.dataLabels = {'Y1'};
XGL3520101ME_DP.SIUnits = { 'M' , char(181)};


XGL3520101ME_graph = componentPlotData(inductor, XGL3520101ME_DP);
XGL3520101ME.addGraph(XGL3520101ME_graph);

indDB.add(XGL3520101ME)


%% XGL3520301ME

XGL3520301ME = inductor('XGL3520-301ME','Power','Composite');
XGL3520301ME.manufacturer = 'Coilcraft';
XGL3520301ME.material = 'Composite';
XGL3520301ME.addParameter('L', 0.3e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3520301ME.addParameter('Rdc', 4.1e-3 ,'typ')
XGL3520301ME.addParameter('Rdc', 4.8e-3 ,'max')
XGL3520301ME.addParameter('SRF', 165000000 ,'typ' )
XGL3520301ME.addParameter('Isat', 4.5 ,'typ', {'10% drop in L'} )
XGL3520301ME.addParameter('Isat', 7.2 ,'typ', {'20% drop in L'} )
XGL3520301ME.addParameter('Isat', 9.7 ,'typ', {'30% drop in L'} )
XGL3520301ME.addParameter('Irms', 14.1 ,'typ', {'20°C rise in L'} )
XGL3520301ME.addParameter('Irms', 19.1 ,'typ', {'40°C rise in L'} )
XGL3520301ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3520301ME.addParameter('Fspice', 300e6 ,'max')
XGL3520301ME.addParameter('R1spice', 7 ,'typ')
XGL3520301ME.addParameter('R2spice', 4.8e-3 ,'typ')
XGL3520301ME.addParameter('Cspice', 3.3e-12 ,'typ')
XGL3520301ME.addParameter('K1spice', 1e-5 ,'typ')
XGL3520301ME.addParameter('K2spice', 0.040 ,'typ')
XGL3520301ME.addParameter('K3spice', 0.3 ,'typ')
XGL3520301ME.addParameter('K4spice', 5.0e-4 ,'typ')
XGL3520301ME.addParameter('K5spice', 1.0e-4 ,'typ')
XGL3520301ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3520301ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3520301ME.addParameter('Height', 2e-3 ,'typ')

XGL3520301ME_DP = digitizedPlot([]);

XGL3520301ME_DP.plotData = {[0, 0.29590271227895265
1.3988000416004704, 0.2910245550363466
2.8686135208701664, 0.28041184465810515
4.328546171142901, 0.26652015650383903
5.749276833200042, 0.25262732075038635
7.356076358380881, 0.23300222705967155
8.717963872143107, 0.21746832446932507
11.29479169281639, 0.18885523394885292
12.38238578698407, 0.17823133447854156
12.90167154994029, 0.17250883113436583

]};

XGL3520301ME_DP.logAxes = [0 0];
XGL3520301ME_DP.normAxes = [0 0];
XGL3520301ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3520301ME_DP.dataLabels = {'Y1'};
XGL3520301ME_DP.SIUnits = { '' , char(181)};


XGL3520301ME_graph = componentPlotData(inductor, XGL3520301ME_DP);
XGL3520301ME.addGraph(XGL3520301ME_graph);

XGL3520301ME_DP = digitizedPlot([]);


XGL3520301ME_DP.plotData = {[0.10, 0.29560157140556265
35.24915393701304, 0.290962119370196
53.92501757033564, 0.2911175671803768
75.7710967709404, 0.29124198522144656
116.82203971712616, 0.30053186418262023

]};

XGL3520301ME_DP.logAxes = [1 1];
XGL3520301ME_DP.normAxes = [0 0];
XGL3520301ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3520301ME_DP.dataLabels = {'Y1'};
XGL3520301ME_DP.SIUnits = { 'M' , char(181)};


XGL3520301ME_graph = componentPlotData(inductor, XGL3520301ME_DP);
XGL3520301ME.addGraph(XGL3520301ME_graph);

indDB.add(XGL3520301ME)



%% XGL3520561ME

XGL3520561ME = inductor('XGL3520-561ME','Power','Composite');
XGL3520561ME.manufacturer = 'Coilcraft';
XGL3520561ME.material = 'Composite';
XGL3520561ME.addParameter('L', 0.56e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3520561ME.addParameter('Rdc', 7.9e-3 ,'typ')
XGL3520561ME.addParameter('Rdc', 9.1e-3 ,'max')
XGL3520561ME.addParameter('SRF', 105000000 ,'typ' )
XGL3520561ME.addParameter('Isat', 3.2 ,'typ', {'10% drop in L'} )
XGL3520561ME.addParameter('Isat', 5.0 ,'typ', {'20% drop in L'} )
XGL3520561ME.addParameter('Isat', 6.9 ,'typ', {'30% drop in L'} )
XGL3520561ME.addParameter('Irms', 9.9 ,'typ', {'20°C rise in L'} )
XGL3520561ME.addParameter('Irms', 13.6 ,'typ', {'40°C rise in L'} )
XGL3520561ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3520561ME.addParameter('Fspice', 300e6 ,'max')
XGL3520561ME.addParameter('R1spice', 3 ,'typ')
XGL3520561ME.addParameter('R2spice', 9.1e-3 ,'typ')
XGL3520561ME.addParameter('Cspice', 4.0e-12 ,'typ')
XGL3520561ME.addParameter('K1spice', 1e-6 ,'typ')
XGL3520561ME.addParameter('K2spice', 0.080 ,'typ')
XGL3520561ME.addParameter('K3spice', 0.56 ,'typ')
XGL3520561ME.addParameter('K4spice', 1.0e-6 ,'typ')
XGL3520561ME.addParameter('K5spice', 1.0e-6 ,'typ')
XGL3520561ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3520561ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3520561ME.addParameter('Height', 2e-3 ,'typ')

XGL3520561ME_DP = digitizedPlot([]);

XGL3520561ME_DP.plotData = {[0, 0.5995906141147013
1.7160543695450194, 0.5768268647828412
3.250859703184471, 0.5359692063753995
4.71860495983986, 0.4901742365268924
5.825694596407649, 0.4574715793929522
8.216026029652854, 0.38223411620790687
9.41534156744715, 0.3528295050563387

]};

XGL3520561ME_DP.logAxes = [0 0];
XGL3520561ME_DP.normAxes = [0 0];
XGL3520561ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3520561ME_DP.dataLabels = {'Y1'};
XGL3520561ME_DP.SIUnits = { '' , char(181)};


XGL3520561ME_graph = componentPlotData(inductor, XGL3520561ME_DP);
XGL3520561ME.addGraph(XGL3520561ME_graph);

XGL3520561ME_DP = digitizedPlot([]);


XGL3520561ME_DP.plotData = {[0.10, 0.6056929684131727
0.44761487502604297, 0.5975332927129521
1.123042424487935, 0.5936273447871598
3.4985879588931335, 0.5944753738503119
12.525987314549862, 0.5908532970171976
44.85373550771152, 0.6442086436818085
67.57790161476764, 0.7180444767745833
]};

XGL3520561ME_DP.logAxes = [1 1];
XGL3520561ME_DP.normAxes = [0 0];
XGL3520561ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3520561ME_DP.dataLabels = {'Y1'};
XGL3520561ME_DP.SIUnits = { 'M' , char(181)};


XGL3520561ME_graph = componentPlotData(inductor, XGL3520561ME_DP);
XGL3520561ME.addGraph(XGL3520561ME_graph);

indDB.add(XGL3520561ME)


%% XGL3520102ME

XGL3520102ME = inductor('XGL3520-102ME','Power','Composite');
XGL3520102ME.manufacturer = 'Coilcraft';
XGL3520102ME.material = 'Composite';
XGL3520102ME.addParameter('L', 1e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL3520102ME.addParameter('Rdc', 12.8e-3 ,'typ')
XGL3520102ME.addParameter('Rdc', 14.8e-3 ,'max')
XGL3520102ME.addParameter('SRF', 80000000 ,'typ' )
XGL3520102ME.addParameter('Isat', 2.4 ,'typ', {'10% drop in L'} )
XGL3520102ME.addParameter('Isat', 3.9 ,'typ', {'20% drop in L'} )
XGL3520102ME.addParameter('Isat', 5.4 ,'typ', {'30% drop in L'} )
XGL3520102ME.addParameter('Irms', 7.5 ,'typ', {'20°C rise in L'} )
XGL3520102ME.addParameter('Irms', 10.1 ,'typ', {'40°C rise in L'} )
XGL3520102ME.addParameter('Fspice', 0.10e6 ,'min')
XGL3520102ME.addParameter('Fspice', 300e6 ,'max')
XGL3520102ME.addParameter('R1spice', 4 ,'typ')
XGL3520102ME.addParameter('R2spice', 1.48e-2 ,'typ')
XGL3520102ME.addParameter('Cspice', 4e-12 ,'typ')
XGL3520102ME.addParameter('K1spice', 1e-6 ,'typ')
XGL3520102ME.addParameter('K2spice', 0.140 ,'typ')
XGL3520102ME.addParameter('K3spice', 1.0 ,'typ')
XGL3520102ME.addParameter('K4spice', 4.07e-3 ,'typ')
XGL3520102ME.addParameter('K5spice', 3.34e-6 ,'typ')
XGL3520102ME.addParameter('Length', 3.5e-3 ,'typ')
XGL3520102ME.addParameter('Width', 3.2e-3 ,'typ')
XGL3520102ME.addParameter('Height', 2e-3 ,'typ')

XGL3520102ME_DP = digitizedPlot([]);

XGL3520102ME_DP.plotData = {[0, 0.9709401709401708
0.7320202500925389, 0.9641025641025639
1.6046135724850579, 0.9321937321937324
2.4641323436504163, 0.8820512820512818
3.225950695830214, 0.8273504273504271
4.111488689801942, 0.7658119658119658
4.853775616827825, 0.7133903133903131
5.706788036105305, 0.6655270655270654
6.553306472985189, 0.6245014245014242
7.41930067231272, 0.5606837606837604
7.777403600696774, 0.5287749287749288

]};

XGL3520102ME_DP.logAxes = [0 0];
XGL3520102ME_DP.normAxes = [0 0];
XGL3520102ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL3520102ME_DP.dataLabels = {'Y1'};
XGL3520102ME_DP.SIUnits = { '' , char(181)};


XGL3520102ME_graph = componentPlotData(inductor, XGL3520102ME_DP);
XGL3520102ME.addGraph(XGL3520102ME_graph);

XGL3520102ME_DP = digitizedPlot([]);


XGL3520102ME_DP.plotData = {[0.10, 0.9923350795510473
0.26484139611563395, 0.9858962660341897
6.448428959975131, 0.9597821768028321
13.03022263925994, 0.960630722947018
22.912188840675693, 1.0225041404177933
44.55217080400996, 1.157788501935463
60.71161112064646, 1.3619144847925386
]};

XGL3520102ME_DP.logAxes = [1 1];
XGL3520102ME_DP.normAxes = [0 0];
XGL3520102ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL3520102ME_DP.dataLabels = {'Y1'};
XGL3520102ME_DP.SIUnits = { 'M' , char(181)};


XGL3520102ME_graph = componentPlotData(inductor, XGL3520102ME_DP);
XGL3520102ME.addGraph(XGL3520102ME_graph);

indDB.add(XGL3520102ME)




%% XGL4018121ME

XGL4018121ME = inductor('XGL4018-121ME','Power','Composite');
XGL4018121ME.manufacturer = 'Coilcraft';
XGL4018121ME.material = 'Composite';
XGL4018121ME.addParameter('L', 0.12e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4018121ME.addParameter('Rdc', 1.6e-3 ,'typ')
XGL4018121ME.addParameter('Rdc', 2.0e-3 ,'max')
XGL4018121ME.addParameter('SRF', 255000000 ,'typ' )
XGL4018121ME.addParameter('Isat', 9.9 ,'typ', {'10% drop in L'} )
XGL4018121ME.addParameter('Isat', 17.4 ,'typ', {'20% drop in L'} )
XGL4018121ME.addParameter('Isat', 24.5 ,'typ', {'30% drop in L'} )
XGL4018121ME.addParameter('Irms', 19.5 ,'typ', {'20°C rise in L'} )
XGL4018121ME.addParameter('Irms', 27.2 ,'typ', {'40°C rise in L'} )
XGL4018121ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4018121ME.addParameter('Fspice', 500e6 ,'max')
XGL4018121ME.addParameter('R1spice', 2 ,'typ')
XGL4018121ME.addParameter('R2spice', 2e-3 ,'typ')
XGL4018121ME.addParameter('Cspice', 3.4e-12 ,'typ')
XGL4018121ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4018121ME.addParameter('K2spice', 0.017 ,'typ')
XGL4018121ME.addParameter('K3spice', 0.12 ,'typ')
XGL4018121ME.addParameter('K4spice', 1e-6 ,'typ')
XGL4018121ME.addParameter('K5spice', 1e-6 ,'typ')
XGL4018121ME.addParameter('Length', 4e-3 ,'typ')
XGL4018121ME.addParameter('Width', 4e-3 ,'typ')
XGL4018121ME.addParameter('Height', 1.8e-3 ,'typ')

XGL4018121ME_DP = digitizedPlot([]);

XGL4018121ME_DP.plotData = {[0, 0.11318526867244816
1.6715976331360924, 0.11240449578911117
4.659763313609463, 0.10926086638907151
7.381656804733724, 0.10611573421829834
10.636094674556206, 0.1009982989470169
12.766272189349106, 0.09784982728572475
14.7189349112426, 0.0947003537772769
17.82544378698224, 0.0899771453617608
20.51775147928993, 0.0856466610312765
22.61834319526627, 0.08210296066706331
26.2869822485207, 0.07698786303914518
27.440828402366854, 0.07541412813207687
]};

XGL4018121ME_DP.logAxes = [0 0];
XGL4018121ME_DP.normAxes = [0 0];
XGL4018121ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4018121ME_DP.dataLabels = {'Y1'};
XGL4018121ME_DP.SIUnits = { '' , char(181)};


XGL4018121ME_graph = componentPlotData(inductor, XGL4018121ME_DP);
XGL4018121ME.addGraph(XGL4018121ME_graph);

XGL4018121ME_DP = digitizedPlot([]);


XGL4018121ME_DP.plotData = {[0.10, 0.11060712133333808
1.1487326845687726, 0.10746234930494815
3.7884987935549987, 0.10760921460992635
8.852792643314842, 0.10272365844666864
12.706168154883365, 0.10276615039835404
19.835596055962675, 0.10120559332804646
26.84306585074303, 0.10124064107794037
]};

XGL4018121ME_DP.logAxes = [1 1];
XGL4018121ME_DP.normAxes = [0 0];
XGL4018121ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4018121ME_DP.dataLabels = {'Y1'};
XGL4018121ME_DP.SIUnits = { 'M' , char(181)};


XGL4018121ME_graph = componentPlotData(inductor, XGL4018121ME_DP);
XGL4018121ME.addGraph(XGL4018121ME_graph);

indDB.add(XGL4018121ME)


%% XGL4018221ME

XGL4018221ME = inductor('XGL4018-221ME','Power','Composite');
XGL4018221ME.manufacturer = 'Coilcraft';
XGL4018221ME.material = 'Composite';
XGL4018221ME.addParameter('L', 0.22e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4018221ME.addParameter('Rdc', 2.6e-3 ,'typ')
XGL4018221ME.addParameter('Rdc', 3.2e-3 ,'max')
XGL4018221ME.addParameter('SRF', 155000000 ,'typ' )
XGL4018221ME.addParameter('Isat', 6.7 ,'typ', {'10% drop in L'} )
XGL4018221ME.addParameter('Isat', 11.1 ,'typ', {'20% drop in L'} )
XGL4018221ME.addParameter('Isat', 16.1 ,'typ', {'30% drop in L'} )
XGL4018221ME.addParameter('Irms', 15.6 ,'typ', {'20°C rise in L'} )
XGL4018221ME.addParameter('Irms', 21.7 ,'typ', {'40°C rise in L'} )
XGL4018221ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4018221ME.addParameter('Fspice', 400e6 ,'max')
XGL4018221ME.addParameter('R1spice', 6 ,'typ')
XGL4018221ME.addParameter('R2spice', 3.2e-3 ,'typ')
XGL4018221ME.addParameter('Cspice', 4.2e-12 ,'typ')
XGL4018221ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4018221ME.addParameter('K2spice', 0.030 ,'typ')
XGL4018221ME.addParameter('K3spice', 0.22 ,'typ')
XGL4018221ME.addParameter('K4spice', 1e-6 ,'typ')
XGL4018221ME.addParameter('K5spice', 1e-6 ,'typ')
XGL4018221ME.addParameter('Length', 4e-3 ,'typ')
XGL4018221ME.addParameter('Width', 4e-3 ,'typ')
XGL4018221ME.addParameter('Height', 1.8e-3 ,'typ')

XGL4018221ME_DP = digitizedPlot([]);

XGL4018221ME_DP.plotData = {[0, 0.2166666666666666
1.7645815722738802, 0.2144444444444444
3.518596787827558, 0.2092592592592592
7.259087066779376, 0.19444444444444436
9.372358410819952, 0.18333333333333326
10.703719357565516, 0.17666666666666658
12.795857988165684, 0.16629629629629625
14.761200338123416, 0.15666666666666657
17.656382079459007, 0.14629629629629623
18.353761622992398, 0.14333333333333323
]};

XGL4018221ME_DP.logAxes = [0 0];
XGL4018221ME_DP.normAxes = [0 0];
XGL4018221ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4018221ME_DP.dataLabels = {'Y1'};
XGL4018221ME_DP.SIUnits = { '' , char(181)};


XGL4018221ME_graph = componentPlotData(inductor, XGL4018221ME_DP);
XGL4018221ME.addGraph(XGL4018221ME_graph);

XGL4018221ME_DP = digitizedPlot([]);


XGL4018221ME_DP.plotData = {[0.10, 0.21487749364301886
4.915931689524404, 0.2091177826080268
12.813395565823393, 0.2160732709177844
16.07697659677429, 0.2230733710581258
20.00298874143695, 0.2376972421618196


]};

XGL4018221ME_DP.logAxes = [1 1];
XGL4018221ME_DP.normAxes = [0 0];
XGL4018221ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4018221ME_DP.dataLabels = {'Y1'};
XGL4018221ME_DP.SIUnits = { 'M' , char(181)};


XGL4018221ME_graph = componentPlotData(inductor, XGL4018221ME_DP);
XGL4018221ME.addGraph(XGL4018221ME_graph);

indDB.add(XGL4018221ME)



%% XGL4018471ME

XGL4018471ME = inductor('XGL4018-471ME','Power','Composite');
XGL4018471ME.manufacturer = 'Coilcraft';
XGL4018471ME.material = 'Composite';
XGL4018471ME.addParameter('L', 0.47e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4018471ME.addParameter('Rdc', 5.1e-3 ,'typ')
XGL4018471ME.addParameter('Rdc', 6.1e-3 ,'max')
XGL4018471ME.addParameter('SRF', 95000000 ,'typ' )
XGL4018471ME.addParameter('Isat', 4.9 ,'typ', {'10% drop in L'} )
XGL4018471ME.addParameter('Isat', 8.1 ,'typ', {'20% drop in L'} )
XGL4018471ME.addParameter('Isat', 11.4 ,'typ', {'30% drop in L'} )
XGL4018471ME.addParameter('Irms', 11.0 ,'typ', {'20°C rise in L'} )
XGL4018471ME.addParameter('Irms', 15.2 ,'typ', {'40°C rise in L'} )
XGL4018471ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4018471ME.addParameter('Fspice', 300e6 ,'max')
XGL4018471ME.addParameter('R1spice', 6 ,'typ')
XGL4018471ME.addParameter('R2spice', 6.1e-3 ,'typ')
XGL4018471ME.addParameter('Cspice', 5.7e-12 ,'typ')
XGL4018471ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4018471ME.addParameter('K2spice', 0.062 ,'typ')
XGL4018471ME.addParameter('K3spice', 0.47 ,'typ')
XGL4018471ME.addParameter('K4spice', 1e-6 ,'typ')
XGL4018471ME.addParameter('K5spice', 1e-6 ,'typ')
XGL4018471ME.addParameter('Length', 4e-3 ,'typ')
XGL4018471ME.addParameter('Width', 4e-3 ,'typ')
XGL4018471ME.addParameter('Height', 1.8e-3 ,'typ')

XGL4018471ME_DP = digitizedPlot([]);

XGL4018471ME_DP.plotData = {[0, 0.43649815043156603
1.8174133558748942, 0.43057953144266353
3.3896872358410826, 0.41726263871763264
5.705832628909552, 0.38766954377311963
7.599323753169907, 0.3580764488286067
9.74640743871513, 0.3270036991368683
12.704987320371936, 0.29297163995067826
11.808960270498732, 0.30184956843403227

]};

XGL4018471ME_DP.logAxes = [0 0];
XGL4018471ME_DP.normAxes = [0 0];
XGL4018471ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4018471ME_DP.dataLabels = {'Y1'};
XGL4018471ME_DP.SIUnits = { '' , char(181)};


XGL4018471ME_graph = componentPlotData(inductor, XGL4018471ME_DP);
XGL4018471ME.addGraph(XGL4018471ME_graph);

XGL4018471ME_DP = digitizedPlot([]);


XGL4018471ME_DP.plotData = {[0.10, 0.4308606699868846
0.20686805475259984, 0.43121719734682756
0.6116406629525076, 0.4249794694890934
1.5286447508535093, 0.4254252284171073
3.9510754990380224, 0.43267538827539664
6.822469547382085, 0.4539777353401821
12.813395565823381, 0.5237817985078824
14.172936913996356, 0.558043864345431

]};

XGL4018471ME_DP.logAxes = [1 1];
XGL4018471ME_DP.normAxes = [0 0];
XGL4018471ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4018471ME_DP.dataLabels = {'Y1'};
XGL4018471ME_DP.SIUnits = { 'M' , char(181)};


XGL4018471ME_graph = componentPlotData(inductor, XGL4018471ME_DP);
XGL4018471ME.addGraph(XGL4018471ME_graph);

indDB.add(XGL4018471ME)



%% XGL4018102ME

XGL4018102ME = inductor('XGL4018-102ME','Power','Composite');
XGL4018102ME.manufacturer = 'Coilcraft';
XGL4018102ME.material = 'Composite';
XGL4018102ME.addParameter('L', 1e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL4018102ME.addParameter('Rdc', 10.7e-3 ,'typ')
XGL4018102ME.addParameter('Rdc', 12.9e-3 ,'max')
XGL4018102ME.addParameter('SRF', 60000000 ,'typ' )
XGL4018102ME.addParameter('Isat', 3.4 ,'typ', {'10% drop in L'} )
XGL4018102ME.addParameter('Isat', 5.6 ,'typ', {'20% drop in L'} )
XGL4018102ME.addParameter('Isat', 7.8 ,'typ', {'30% drop in L'} )
XGL4018102ME.addParameter('Irms', 7.3 ,'typ', {'20°C rise in L'} )
XGL4018102ME.addParameter('Irms', 10.4 ,'typ', {'40°C rise in L'} )
XGL4018102ME.addParameter('Fspice', 0.10e6 ,'min')
XGL4018102ME.addParameter('Fspice', 250e6 ,'max')
XGL4018102ME.addParameter('R1spice', 6 ,'typ')
XGL4018102ME.addParameter('R2spice', 1.29e-2 ,'typ')
XGL4018102ME.addParameter('Cspice', 6.6e-12 ,'typ')
XGL4018102ME.addParameter('K1spice', 1e-6 ,'typ')
XGL4018102ME.addParameter('K2spice', 0.150 ,'typ')
XGL4018102ME.addParameter('K3spice', 1.0 ,'typ')
XGL4018102ME.addParameter('K4spice', 1e-6 ,'typ')
XGL4018102ME.addParameter('K5spice', 1e-6 ,'typ')
XGL4018102ME.addParameter('Length', 4e-3 ,'typ')
XGL4018102ME.addParameter('Width', 4e-3 ,'typ')
XGL4018102ME.addParameter('Height', 1.8e-3 ,'typ')

XGL4018102ME_DP = digitizedPlot([]);

XGL4018102ME_DP.plotData = {[0, 0.9588090199409075
0.303600020229227, 0.9616857410544326
0.8363585932762858, 0.9585687718096346
1.62276455318175, 0.9405801929805626
2.43449875574617, 0.910746880762788
3.4237809153933467, 0.869023788631696
4.184754358277327, 0.8332868791048258
5.647504688436027, 0.7618305781439896
6.349292003147001, 0.7290704931186747
7.14409204923627, 0.6933235732530006
7.854346442184822, 0.6635202920516349
8.412422830454581, 0.6455992830094837
]};

XGL4018102ME_DP.logAxes = [0 0];
XGL4018102ME_DP.normAxes = [0 0];
XGL4018102ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL4018102ME_DP.dataLabels = {'Y1'};
XGL4018102ME_DP.SIUnits = { '' , char(181)};


XGL4018102ME_graph = componentPlotData(inductor, XGL4018102ME_DP);
XGL4018102ME.addGraph(XGL4018102ME_graph);

XGL4018102ME_DP = digitizedPlot([]);


XGL4018102ME_DP.plotData = {[0.10, 0.9804340766909574
0.19505001381205572, 0.9811793034263102
0.6652619749230598, 0.9825580489861621
1.5158524649961744, 0.9834845948432784
4.915931689524404, 1.0326501131336618
6.822469547382096, 1.1539440937624412
9.232694536044212, 1.4633204155934525
10.212312304989291, 1.6608299876492254
]};

XGL4018102ME_DP.logAxes = [1 1];
XGL4018102ME_DP.normAxes = [0 0];
XGL4018102ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL4018102ME_DP.dataLabels = {'Y1'};
XGL4018102ME_DP.SIUnits = { 'M' , char(181)};


XGL4018102ME_graph = componentPlotData(inductor, XGL4018102ME_DP);
XGL4018102ME.addGraph(XGL4018102ME_graph);

indDB.add(XGL4018102ME)


%% XGL5020161ME

XGL5020161ME = inductor('XGL5020-161ME','Power','Composite');
XGL5020161ME.manufacturer = 'Coilcraft';
XGL5020161ME.material = 'Composite';
XGL5020161ME.addParameter('L', 0.16e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL5020161ME.addParameter('Rdc', 1.8e-3 ,'typ')
XGL5020161ME.addParameter('Rdc', 2.1e-3 ,'max')
XGL5020161ME.addParameter('SRF', 150000000 ,'typ' )
XGL5020161ME.addParameter('Isat', 11.5 ,'typ', {'10% drop in L'} )
XGL5020161ME.addParameter('Isat', 19.0 ,'typ', {'20% drop in L'} )
XGL5020161ME.addParameter('Isat', 27.0 ,'typ', {'30% drop in L'} )
XGL5020161ME.addParameter('Irms', 21.2 ,'typ', {'20°C rise in L'} )
XGL5020161ME.addParameter('Irms', 29.4 ,'typ', {'40°C rise in L'} )
XGL5020161ME.addParameter('Fspice', 0.10e6 ,'min')
XGL5020161ME.addParameter('Fspice', 250e6 ,'max')
XGL5020161ME.addParameter('R1spice', 5 ,'typ')
XGL5020161ME.addParameter('R2spice', 1.8e-3 ,'typ')
XGL5020161ME.addParameter('Cspice', 7.0e-12 ,'typ')
XGL5020161ME.addParameter('K1spice', 2e-4 ,'typ')
XGL5020161ME.addParameter('K2spice', 0.027 ,'typ')
XGL5020161ME.addParameter('K3spice', 0.16 ,'typ')
XGL5020161ME.addParameter('K4spice', 2e-3 ,'typ')
XGL5020161ME.addParameter('K5spice', 1e-6 ,'typ')
XGL5020161ME.addParameter('Length', 5.28e-3 ,'typ')
XGL5020161ME.addParameter('Width', 5.28e-3 ,'typ')
XGL5020161ME.addParameter('Height', 2.1e-3 ,'typ')

XGL5020161ME_DP = digitizedPlot([]);

XGL5020161ME_DP.plotData = {[0, 0.15586673369410567
11.640238313964757, 0.1415074854249239
15.753799838120164, 0.13433839383357085
20.594713805807753, 0.12461364310673607
23.784639124720663, 0.11796064613403884
28.5416640632838, 0.10874663449540004
35.81707364532749, 0.09441452091388201

]};

XGL5020161ME_DP.logAxes = [0 0];
XGL5020161ME_DP.normAxes = [0 0];
XGL5020161ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL5020161ME_DP.dataLabels = {'Y1'};
XGL5020161ME_DP.SIUnits = { '' , char(181)};


XGL5020161ME_graph = componentPlotData(inductor, XGL5020161ME_DP);
XGL5020161ME.addGraph(XGL5020161ME_graph);

XGL5020161ME_DP = digitizedPlot([]);


XGL5020161ME_DP.plotData = {[0.10, 0.16574449605296718
1.2983544788836978, 0.16105070836904536
8.409958879273816, 0.1613029107302274
34.544082222952326, 0.16149388795520056
114.31673879673103, 0.17461928986574818

]};

XGL5020161ME_DP.logAxes = [1 1];
XGL5020161ME_DP.normAxes = [0 0];
XGL5020161ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL5020161ME_DP.dataLabels = {'Y1'};
XGL5020161ME_DP.SIUnits = { 'M' , char(181)};


XGL5020161ME_graph = componentPlotData(inductor, XGL5020161ME_DP);
XGL5020161ME.addGraph(XGL5020161ME_graph);

indDB.add(XGL5020161ME)


%% XGL5020471ME

XGL5020471ME = inductor('XGL5020-471ME','Power','Composite');
XGL5020471ME.manufacturer = 'Coilcraft';
XGL5020471ME.material = 'Composite';
XGL5020471ME.addParameter('L', 0.47e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL5020471ME.addParameter('Rdc', 3.7e-3 ,'typ')
XGL5020471ME.addParameter('Rdc', 4.3e-3 ,'max')
XGL5020471ME.addParameter('SRF', 75000000 ,'typ' )
XGL5020471ME.addParameter('Isat', 6.4 ,'typ', {'10% drop in L'} )
XGL5020471ME.addParameter('Isat', 10.7 ,'typ', {'20% drop in L'} )
XGL5020471ME.addParameter('Isat', 15.7 ,'typ', {'30% drop in L'} )
XGL5020471ME.addParameter('Irms', 16.0 ,'typ', {'20°C rise in L'} )
XGL5020471ME.addParameter('Irms', 22.1 ,'typ', {'40°C rise in L'} )
XGL5020471ME.addParameter('Fspice', 0.10e6 ,'min')
XGL5020471ME.addParameter('Fspice', 250e6 ,'max')
XGL5020471ME.addParameter('R1spice', 4 ,'typ')
XGL5020471ME.addParameter('R2spice', 3.7e-3 ,'typ')
XGL5020471ME.addParameter('Cspice', 9.0e-12 ,'typ')
XGL5020471ME.addParameter('K1spice', 1e-4 ,'typ')
XGL5020471ME.addParameter('K2spice', 0.085 ,'typ')
XGL5020471ME.addParameter('K3spice', 0.47 ,'typ')
XGL5020471ME.addParameter('K4spice', 2e-3 ,'typ')
XGL5020471ME.addParameter('K5spice', 1e-6 ,'typ')
XGL5020471ME.addParameter('Length', 5.28e-3 ,'typ')
XGL5020471ME.addParameter('Width', 5.28e-3 ,'typ')
XGL5020471ME.addParameter('Height', 2.1e-3 ,'typ')

XGL5020471ME_DP = digitizedPlot([]);

XGL5020471ME_DP.plotData = {[0, 0.4716479017400206
1.5529171647436333, 0.47328556806550665
3.5464780543619012, 0.46018423746161724
5.854762725620362, 0.4405322415557831
9.404561293091936, 0.40450358239508705
12.587098580395615, 0.36847492323439107
14.87782301790428, 0.3422722620266123
17.6057453647964, 0.3144319344933473
19.96653173390213, 0.29805527123848535

]};

XGL5020471ME_DP.logAxes = [0 0];
XGL5020471ME_DP.normAxes = [0 0];
XGL5020471ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL5020471ME_DP.dataLabels = {'Y1'};
XGL5020471ME_DP.SIUnits = { '' , char(181)};


XGL5020471ME_graph = componentPlotData(inductor, XGL5020471ME_DP);
XGL5020471ME.addGraph(XGL5020471ME_graph);

XGL5020471ME_DP = digitizedPlot([]);


XGL5020471ME_DP.plotData = {[0.10, 0.4659592637498647
1.1138958502593403, 0.45974635395918245
6.526326174062693, 0.47120687017772267
11.64510447369847, 0.4787650029534631
23.511646527746613, 0.5056250055170035
35.67666942590808, 0.5548571985555075
54.985416242788254, 0.6835823259604876

]};

XGL5020471ME_DP.logAxes = [1 1];
XGL5020471ME_DP.normAxes = [0 0];
XGL5020471ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL5020471ME_DP.dataLabels = {'Y1'};
XGL5020471ME_DP.SIUnits = { 'M' , char(181)};


XGL5020471ME_graph = componentPlotData(inductor, XGL5020471ME_DP);
XGL5020471ME.addGraph(XGL5020471ME_graph);

indDB.add(XGL5020471ME)


%% XGL5020102ME

XGL5020102ME = inductor('XGL5020-102ME','Power','Composite');
XGL5020102ME.manufacturer = 'Coilcraft';
XGL5020102ME.material = 'Composite';
XGL5020102ME.addParameter('L', 1.0e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL5020102ME.addParameter('Rdc', 7.5e-3 ,'typ')
XGL5020102ME.addParameter('Rdc', 8.7e-3 ,'max')
XGL5020102ME.addParameter('SRF', 50000000 ,'typ' )
XGL5020102ME.addParameter('Isat', 4.9 ,'typ', {'10% drop in L'} )
XGL5020102ME.addParameter('Isat', 8.1 ,'typ', {'20% drop in L'} )
XGL5020102ME.addParameter('Isat', 11.4 ,'typ', {'30% drop in L'} )
XGL5020102ME.addParameter('Irms', 11.0 ,'typ', {'20°C rise in L'} )
XGL5020102ME.addParameter('Irms', 15.0 ,'typ', {'40°C rise in L'} )
XGL5020102ME.addParameter('Fspice', 0.10e6 ,'min')
XGL5020102ME.addParameter('Fspice', 150e6 ,'max')
XGL5020102ME.addParameter('R1spice', 4 ,'typ')
XGL5020102ME.addParameter('R2spice', 7.5e-3 ,'typ')
XGL5020102ME.addParameter('Cspice', 10.5e-12 ,'typ')
XGL5020102ME.addParameter('K1spice', 1.5e-4 ,'typ')
XGL5020102ME.addParameter('K2spice', 0.180 ,'typ')
XGL5020102ME.addParameter('K3spice', 1.0 ,'typ')
XGL5020102ME.addParameter('K4spice', 6.5e-3 ,'typ')
XGL5020102ME.addParameter('K5spice', 6e-7 ,'typ')
XGL5020102ME.addParameter('Length', 5.28e-3 ,'typ')
XGL5020102ME.addParameter('Width', 5.28e-3 ,'typ')
XGL5020102ME.addParameter('Height', 2.1e-3 ,'typ')

XGL5020102ME_DP = digitizedPlot([]);

XGL5020102ME_DP.plotData = {[0, 0.955599102585595
2.3516637478108593, 0.9285568249506885
3.801751313485114, 0.891798502644804
5.264448336252192, 0.8477073487396607
7.130647985989494, 0.7962568260207373
8.94640980735552, 0.7374777518984451
10.736952714535906, 0.6860323654145903
12.666199649737308, 0.6296895787930558
14.998949211908934, 0.5806513744779053

]};

XGL5020102ME_DP.logAxes = [0 0];
XGL5020102ME_DP.normAxes = [0 0];
XGL5020102ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL5020102ME_DP.dataLabels = {'Y1'};
XGL5020102ME_DP.SIUnits = { '' , char(181)};


XGL5020102ME_graph = componentPlotData(inductor, XGL5020102ME_DP);
XGL5020102ME.addGraph(XGL5020102ME_graph);

XGL5020102ME_DP = digitizedPlot([]);


XGL5020102ME_DP.plotData = {[0.10, 0.9621857616534344
0.33951439624743923, 0.9631505667471425
1.5901604216115177, 0.9718648955116697
9.461493310166398, 0.9733175875942015
15.994669401496443, 1.0120359910939867
24.64952344725096, 1.1722084917730633
37.123955284214546, 1.5720161205461665

]};

XGL5020102ME_DP.logAxes = [1 1];
XGL5020102ME_DP.normAxes = [0 0];
XGL5020102ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL5020102ME_DP.dataLabels = {'Y1'};
XGL5020102ME_DP.SIUnits = { 'M' , char(181)};


XGL5020102ME_graph = componentPlotData(inductor, XGL5020102ME_DP);
XGL5020102ME.addGraph(XGL5020102ME_graph);

indDB.add(XGL5020102ME)


%% XGL5020222ME

XGL5020222ME = inductor('XGL5020-222ME','Power','Composite');
XGL5020222ME.manufacturer = 'Coilcraft';
XGL5020222ME.material = 'Composite';
XGL5020222ME.addParameter('L', 2.2e-6 ,'typ', {'1 MHz, 0.1 Vrms, 0 Adc'} )
XGL5020222ME.addParameter('Rdc', 16.3e-3 ,'typ')
XGL5020222ME.addParameter('Rdc', 18.8e-3 ,'max')
XGL5020222ME.addParameter('SRF', 30000000 ,'typ' )
XGL5020222ME.addParameter('Isat', 3.3 ,'typ', {'10% drop in L'} )
XGL5020222ME.addParameter('Isat', 5.4 ,'typ', {'20% drop in L'} )
XGL5020222ME.addParameter('Isat', 7.6 ,'typ', {'30% drop in L'} )
XGL5020222ME.addParameter('Irms', 7.8 ,'typ', {'20°C rise in L'} )
XGL5020222ME.addParameter('Irms', 10.8 ,'typ', {'40°C rise in L'} )
XGL5020222ME.addParameter('Fspice', 0.10e6 ,'min')
XGL5020222ME.addParameter('Fspice', 150e6 ,'max')
XGL5020222ME.addParameter('R1spice', 5 ,'typ')
XGL5020222ME.addParameter('R2spice', 1.63e-2 ,'typ')
XGL5020222ME.addParameter('Cspice', 10.0e-12 ,'typ')
XGL5020222ME.addParameter('K1spice', 2.0e-4 ,'typ')
XGL5020222ME.addParameter('K2spice', 0.50 ,'typ')
XGL5020222ME.addParameter('K3spice', 2.2 ,'typ')
XGL5020222ME.addParameter('K4spice', 7.0e-3 ,'typ')
XGL5020222ME.addParameter('K5spice', 1e-6 ,'typ')
XGL5020222ME.addParameter('Length', 5.28e-3 ,'typ')
XGL5020222ME.addParameter('Width', 5.28e-3 ,'typ')
XGL5020222ME.addParameter('Height', 2.1e-3 ,'typ')

XGL5020222ME_DP = digitizedPlot([]);

XGL5020222ME_DP.plotData = {[0, 2.320704252425736
1.2516625831291552, 2.271838079608898
2.3185159257962895, 2.173750695731507
4.4018200910045495, 1.9529728535607107
5.4434721736086775, 1.8302888505080994
6.627931396569825, 1.7014938861697184
7.879593979698982, 1.572716135806791
9.307665383269159, 1.4562786540966408
10.046902345117251, 1.40728767995777

]};

XGL5020222ME_DP.logAxes = [0 0];
XGL5020222ME_DP.normAxes = [0 0];
XGL5020222ME_DP.axisLabels = {'I, Current', 'L, Inductance'};
XGL5020222ME_DP.dataLabels = {'Y1'};
XGL5020222ME_DP.SIUnits = { '' , char(181)};


XGL5020222ME_graph = componentPlotData(inductor, XGL5020222ME_DP);
XGL5020222ME.addGraph(XGL5020222ME_graph);

XGL5020222ME_DP = digitizedPlot([]);


XGL5020222ME_DP.plotData = {[0.10, 2.336217656185142
1.4508726359560018, 2.323442966144176
4.655005303685238, 2.3437221757508047
11.667312988237201, 2.633247004440903
17.30108462094412, 3.269168049441129
24.68865624617203, 4.8464177943619084
27.513545798214643, 6.061954663211667
]};

XGL5020222ME_DP.logAxes = [1 1];
XGL5020222ME_DP.normAxes = [0 0];
XGL5020222ME_DP.axisLabels = {'F, Frequency', 'L, Inductance'};
XGL5020222ME_DP.dataLabels = {'Y1'};
XGL5020222ME_DP.SIUnits = { 'M' , char(181)};


XGL5020222ME_graph = componentPlotData(inductor, XGL5020222ME_DP);
XGL5020222ME.addGraph(XGL5020222ME_graph);

indDB.add(XGL5020222ME)



% indDB(end).datasheet