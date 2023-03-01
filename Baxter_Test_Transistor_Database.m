%% Test Cases adding and subtracting devices to transistor database

clear

transDB = transistorDB();

%% EPC2055
EPC2055 = transistor('EPC2055','en-nMOS','GaN');
EPC2055.manufacturer = 'EPC';
EPC2055.material = 'GaN';
EPC2055.addParameter('Vds', 40 ,'max', {'VGS = 0 V, ID = 0.5 mA (Continuous)'} )
EPC2055.addParameter('Ids', 29 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2055.addParameter('Idspulse', 161 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2055.addParameter('Ron', 0.003 ,'typ', {'VGS = 5 V, ID = 15 A'} )
EPC2055.addParameter('Ron', 0.0036 ,'max', {'VGS = 5 V, ID = 15 A'} )
EPC2055.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 7 mA'} )
EPC2055.addParameter('Vth', 1.1 ,'typ', {'VDS = VGS, ID = 7 mA'} )
EPC2055.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 7 mA'} )
EPC2055.addParameter('Vsd', 1.9 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2055.addParameter('Coss', 574e-12 ,'typ', {'VDS = 20 V, VGS = 0 V'} )
EPC2055.addParameter('Coss', 612e-12 ,'max', {'VDS = 20 V, VGS = 0 V'} )
EPC2055.addParameter('Rg', 0.4 ,'typ')
EPC2055.addParameter('Qg', 6.6e-9 ,'typ', {'VDS = 20 V, VGS = 5 V, ID = 15 A'} )
EPC2055.addParameter('Qg', 8.5e-9 ,'max', {'VDS = 20 V, VGS = 5 V, ID = 15 A'} )
EPC2055.addParameter('Length', 2.5e-3 ,'typ')
EPC2055.addParameter('Width', 1.5e-3 ,'typ')
EPC2055.addParameter('Height', 0.638e-3 ,'typ')

transDB.add(EPC2055)


%% EPC2069
EPC2069 = transistor('EPC2069','en-nMOS','GaN');
EPC2069.manufacturer = 'EPC';
EPC2069.material = 'GaN';
EPC2069.addParameter('Vds', 40 ,'max', {'VGS = 0 V, ID = 0.7 mA (Continuous)' } )
EPC2069.addParameter('Ids', 80 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2069.addParameter('Idspulse', 422 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2069.addParameter('Ron', 0.0016 ,'typ', {'VGS = 5 V, ID = 30 A'} )
EPC2069.addParameter('Ron', 0.00225 ,'max', {'VGS = 5 V, ID = 30 A'} )
EPC2069.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 14 mA'} )
EPC2069.addParameter('Vth', 1.4 ,'typ', {'VDS = VGS, ID = 14 mA'} )
EPC2069.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 14 mA'} )
EPC2069.addParameter('Vsd', 1.6 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2069.addParameter('Coss', 1044e-12 ,'typ', {'VDS = 20 V, VGS = 0 V'} )
EPC2069.addParameter('Coss', 1566e-12 ,'max', {'VDS = 20 V, VGS = 0 V'} )
EPC2069.addParameter('Rg', 0.4 ,'typ')
EPC2069.addParameter('Qg', 12.5e-9 ,'typ', {'VDS = 20 V, VGS = 5 V, ID = 30 A'} )
EPC2069.addParameter('Qg', 16.2e-9 ,'max', {'VDS = 20 V, VGS = 5 V, ID = 30 A'} )
EPC2069.addParameter('Length', 3.25e-3 ,'typ')
EPC2069.addParameter('Width', 3.25e-3 ,'typ')
EPC2069.addParameter('Height', 0.638e-3 ,'typ')


transDB.add(EPC2069)


%% EPC2066
EPC2066 = transistor('EPC2066','en-nMOS','GaN');
EPC2066.manufacturer = 'EPC';
EPC2066.material = 'GaN';
EPC2066.addParameter('Vds', 40 ,'max', {'VGS = 0 V, ID = 1.2 mA (Continuous)'} )
EPC2066.addParameter('Ids', 90 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2066.addParameter('Idspulse', 639 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2066.addParameter('Ron', 0.0008 ,'typ', {'VGS = 5 V, ID = 50 A'} )
EPC2066.addParameter('Ron', 0.0011 ,'max', {'VGS = 5 V, ID = 50 A'} )
EPC2066.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 28 mA'} )
EPC2066.addParameter('Vth', 1.2 ,'typ', {'VDS = VGS, ID = 28 mA'} )
EPC2066.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 28 mA'} )
EPC2066.addParameter('Vsd', 1.5 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2066.addParameter('Coss', 1670e-12 ,'typ', {'VDS = 20 V, VGS = 0 V'} )
EPC2066.addParameter('Coss', 1919e-12 ,'max', {'VDS = 20 V, VGS = 0 V'} )
EPC2066.addParameter('Rg', 0.4 ,'typ')
EPC2066.addParameter('Qg', 25e-9 ,'typ', {'VDS = 20 V, VGS = 5 V, ID = 50 A'} )
EPC2066.addParameter('Qg', 33e-9 ,'max', {'VDS = 20 V, VGS = 5 V, ID = 50 A'} )
EPC2066.addParameter('Length', 6.05e-3 ,'typ')
EPC2066.addParameter('Width', 2.3e-3 ,'typ')
EPC2066.addParameter('Height', 0.618e-3 ,'typ')

transDB.add(EPC2066)


%% EPC2067
EPC2067 = transistor('EPC2067','en-nMOS','GaN');
EPC2067.manufacturer = 'EPC';
EPC2067.material = 'GaN';
EPC2067.addParameter('Vds', 40 ,'max', {'VGS = 0 V, ID = 1.1 mA (Continuous)'} )
EPC2067.addParameter('Ids', 69 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2067.addParameter('Idspulse', 409 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2067.addParameter('Ron', 1.3e-3 ,'typ', {'VGS = 5 V, ID = 37 A'} )
EPC2067.addParameter('Ron', 1.55e-3 ,'max', {'VGS = 5 V, ID = 37 A'} )
EPC2067.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 18 mA'} )
EPC2067.addParameter('Vth', 1 ,'typ', {'VDS = VGS, ID = 18 mA'} )
EPC2067.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 18 mA'} )
EPC2067.addParameter('Vsd', 1.2 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2067.addParameter('Coss', 1071e-12 ,'typ', {'VDS = 20 V, VGS = 0 V'} )
EPC2067.addParameter('Coss', 1607e-12 ,'max', {'VDS = 20 V, VGS = 0 V'} )
EPC2067.addParameter('Rg', 0.4 ,'typ')
EPC2067.addParameter('Qg', 17.1e-9 ,'typ', {'VDS = 20 V, VGS = 5 V, ID = 37 A'} )
EPC2067.addParameter('Qg', 22.3e-9 ,'max', {'VDS = 20 V, VGS = 5 V, ID = 37 A'} )
EPC2067.addParameter('Length', 6.05e-3 ,'typ')
EPC2067.addParameter('Width', 2.3e-3 ,'typ')
EPC2067.addParameter('Height', 0.618e-3 ,'typ')

transDB.add(EPC2067)



%% EPC2065
EPC2065 = transistor('EPC2065','en-nMOS','GaN');
EPC2065.manufacturer = 'EPC';
EPC2065.material = 'GaN';
EPC2065.addParameter('Vds', 80 ,'max', {'VGS = 0 V, ID = 0.4 mA (Continuous)'} )
EPC2065.addParameter('Ids', 60 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2065.addParameter('Idspulse', 215 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2065.addParameter('Ron', 2.7e-3 ,'typ', {'VGS = 5 V, ID = 25 A'} )
EPC2065.addParameter('Ron', 3.6e-3 ,'max', {'VGS = 5 V, ID = 25 A'} )
EPC2065.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 7 mA'} )
EPC2065.addParameter('Vth', 1.2 ,'typ', {'VDS = VGS, ID = 7 mA'} )
EPC2065.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 7 mA'} )
EPC2065.addParameter('Vsd', 1.4 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2065.addParameter('Coss', 534e-12 ,'typ', {'VDS = 40 V, VGS = 0 V'} )
EPC2065.addParameter('Coss', 801e-12 ,'max', {'VDS = 40 V, VGS = 0 V'} )
EPC2065.addParameter('Rg', 0.5 ,'typ')
EPC2065.addParameter('Qg', 9.4e-9 ,'typ', {'VDS = 40 V, VGS = 5 V, ID = 25 A'} )
EPC2065.addParameter('Qg', 12.2e-9 ,'max', {'VDS = 40 V, VGS = 5 V, ID = 25 A'} )
EPC2065.addParameter('Length', 3.5e-3 ,'typ')
EPC2065.addParameter('Width', 1.95e-3 ,'typ')
EPC2065.addParameter('Height', 0.638e-3 ,'typ')

transDB.add(EPC2065)


%% EPC2051
EPC2051 = transistor('EPC2051','en-nMOS','GaN');
EPC2051.manufacturer = 'EPC';
EPC2051.material = 'GaN';
EPC2051.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = 300 μA (Continuous)'} )
EPC2051.addParameter('Ids', 1.7 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2051.addParameter('Idspulse', 37 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2051.addParameter('Ron', 20e-3 ,'typ', {'VGS = 5 V, ID = 3 A'} )
EPC2051.addParameter('Ron', 25e-3 ,'max', {'VGS = 5 V, ID = 3 A'} )
EPC2051.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 1.5 mA'} )
EPC2051.addParameter('Vth', 1.4 ,'typ', {'VDS = VGS, ID = 1.5 mA'} )
EPC2051.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 1.5 mA'} )
EPC2051.addParameter('Vsd', 1.9 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2051.addParameter('Coss', 86e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2051.addParameter('Coss', 129e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2051.addParameter('Rg', 0.8 ,'typ')
EPC2051.addParameter('Qg', 1.8e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 3 A'} )
EPC2051.addParameter('Qg', 2.3e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 3 A'} )
EPC2051.addParameter('Length', 1.3e-3 ,'typ')
EPC2051.addParameter('Width', 0.85e-3 ,'typ')
EPC2051.addParameter('Height', 0.685e-3 ,'typ')

transDB.add(EPC2051)


%% EPC2052
EPC2052 = transistor('EPC2052','en-nMOS','GaN');
EPC2052.manufacturer = 'EPC';
EPC2052.material = 'GaN';
EPC2052.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = 0.2 mA (Continuous)'} )
EPC2052.addParameter('Ids', 8.2 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2052.addParameter('Idspulse', 74 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2052.addParameter('Ron', 10e-3 ,'typ', {'VGS = 5 V, ID = 11 A'} )
EPC2052.addParameter('Ron', 13.5e-3 ,'max', {'VGS = 5 V, ID = 11 A'} )
EPC2052.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 3 mA'} )
EPC2052.addParameter('Vth', 1.4 ,'typ', {'VDS = VGS, ID = 3 mA'} )
EPC2052.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 3 mA'} )
EPC2052.addParameter('Vsd', 2 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2052.addParameter('Coss', 195e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2052.addParameter('Coss', 293e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2052.addParameter('Rg', 0.7 ,'typ')
EPC2052.addParameter('Qg', 3.5e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 11 A'} )
EPC2052.addParameter('Qg', 4.5e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 11 A'} )
EPC2052.addParameter('Length', 1.5e-3 ,'typ')
EPC2052.addParameter('Width', 1.5e-3 ,'typ')
EPC2052.addParameter('Height', 0.885e-3 ,'typ')

transDB.add(EPC2052)



%% EPC2044
EPC2044 = transistor('EPC2044','en-nMOS','GaN');
EPC2044.manufacturer = 'EPC';
EPC2044.material = 'GaN';
EPC2044.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = 0.2 mA (Continuous)'} )
EPC2044.addParameter('Ids', 9.4 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2044.addParameter('Idspulse', 89 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2044.addParameter('Ron', 7e-3 ,'typ', {'VGS = 5 V, ID = 10 A'} )
EPC2044.addParameter('Ron', 10.5e-3 ,'max', {'VGS = 5 V, ID = 10 A'} )
EPC2044.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 3 mA'} )
EPC2044.addParameter('Vth', 1.4 ,'typ', {'VDS = VGS, ID = 3 mA'} )
EPC2044.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 3 mA'} )
EPC2044.addParameter('Vsd', 2 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2044.addParameter('Coss', 196e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2044.addParameter('Coss', 294e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2044.addParameter('Rg', 0.5 ,'typ')
EPC2044.addParameter('Qg', 4.3e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 10 A'} )
EPC2044.addParameter('Qg', 5.5e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 10 A'} )
EPC2044.addParameter('Length', 2.15e-3 ,'typ')
EPC2044.addParameter('Width', 1.25e-3 ,'typ')
EPC2044.addParameter('Height', 0.518e-3 ,'typ')

transDB.add(EPC2044)


%% EPC2204
EPC2204 = transistor('EPC2204','en-nMOS','GaN');
EPC2204.manufacturer = 'EPC';
EPC2204.material = 'GaN';
EPC2204.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = 0.25 mA (Continuous)'} )
EPC2204.addParameter('Ids', 29 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2204.addParameter('Idspulse', 125 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2204.addParameter('Ron', 4.4e-3 ,'typ', {'VGS = 5 V, ID = 16 A'} )
EPC2204.addParameter('Ron', 6e-3 ,'max', {'VGS = 5 V, ID = 16 A'} )
EPC2204.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 4 mA'} )
EPC2204.addParameter('Vth', 1.1 ,'typ', {'VDS = VGS, ID = 4 mA'} )
EPC2204.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 4 mA'} )
EPC2204.addParameter('Vsd', 1.6 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2204.addParameter('Coss', 304e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2204.addParameter('Coss', 456e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2204.addParameter('Rg', 0.4 ,'typ')
EPC2204.addParameter('Qg', 5.7e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 16 A'} )
EPC2204.addParameter('Qg', 7.4e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 16 A'} )
EPC2204.addParameter('Length', 2.5e-3 ,'typ')
EPC2204.addParameter('Width', 1.5e-3 ,'typ')
EPC2204.addParameter('Height', 0.638e-3 ,'typ')

transDB.add(EPC2204)


%% EPC2088
EPC2088 = transistor('EPC2088','en-nMOS','GaN');
EPC2088.manufacturer = 'EPC';
EPC2088.material = 'GaN';
EPC2088.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = 0.1 mA (Continuous)'} )
EPC2088.addParameter('Ids', 60 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2088.addParameter('Idspulse', 231 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2088.addParameter('Ron', 2.4e-3 ,'typ', {'VGS = 5 V, ID = 25 A'} )
EPC2088.addParameter('Ron', 3.2e-3 ,'max', {'VGS = 5 V, ID = 25 A'} )
EPC2088.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 7 mA'} )
EPC2088.addParameter('Vth', 1.3 ,'typ', {'VDS = VGS, ID = 7 mA'} )
EPC2088.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 7 mA'} )
EPC2088.addParameter('Vsd', 1.5 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2088.addParameter('Coss', 557e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2088.addParameter('Coss', 659e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2088.addParameter('Rg', 0.4 ,'typ')
EPC2088.addParameter('Qg', 12.5e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 25 A'} )
EPC2088.addParameter('Qg', 17.8e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 25 A'} )
EPC2088.addParameter('Length', 3.5e-3 ,'typ')
EPC2088.addParameter('Width', 1.95e-3 ,'typ')
EPC2088.addParameter('Height', 0.638e-3 ,'typ')


EPC2088_DP = digitizedPlot([]);


EPC2088_DP.plotData = {[0, 1889.6801046574656
7.80952380952381, 1645.6672543793616
10.85714285714284, 1367.8614663685528
13.269841269841283, 1020.5158638351236
19.619047619047638, 875.0458085098675
25.46031746031747, 810.4489030540457
29.142857142857167, 673.6696729953237
34.73015873015875, 623.9263661915011
39.93650793650792, 595.9464721366903
45.142857142857146, 578.0712030014258
55.17460317460318, 560.9408550074761
68.8888888888889, 544.4727854776243
84.50793650793652, 536.783627448352
99.49206349206351, 537.4044293999186
    ]};

EPC2088_DP.logAxes = [0 1];
EPC2088_DP.normAxes = [0 0];
EPC2088_DP.axisLabels = {'Vds, Drain-Source Voltage', 'Coss, Output Capacitance'};
EPC2088_DP.dataLabels = {'Y1'};
EPC2088_DP.SIUnits = { '' , 'p'};


EPC2088_graph = componentPlotData(transistor, EPC2088_DP);
EPC2088.addGraph(EPC2088_graph);


transDB.add(EPC2088)


%% EPC2071
EPC2071 = transistor('EPC2071','en-nMOS','GaN');
EPC2071.manufacturer = 'EPC';
EPC2071.material = 'GaN';
EPC2071.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = 0.15 mA (Continuous)'} )
EPC2071.addParameter('Ids', 64 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2071.addParameter('Idspulse', 350 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2071.addParameter('Ron', 1.7e-3 ,'typ', {'VGS = 5 V, ID = 30 A'} )
EPC2071.addParameter('Ron', 2.2e-3 ,'max', {'VGS = 5 V, ID = 30 A'} )
EPC2071.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 13 mA'} )
EPC2071.addParameter('Vth', 1.3 ,'typ', {'VDS = VGS, ID = 13 mA'} )
EPC2071.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 13 mA'} )
EPC2071.addParameter('Vsd', 1.5 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2071.addParameter('Coss', 878e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2071.addParameter('Coss', 976e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2071.addParameter('Rg', 0.3 ,'typ')
EPC2071.addParameter('Qg', 18e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 30 A'} )
EPC2071.addParameter('Qg', 26e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 30 A'} )
EPC2071.addParameter('Length', 4.45e-3 ,'typ')
EPC2071.addParameter('Width', 2.3e-3 ,'typ')
EPC2071.addParameter('Height', 0.618e-3 ,'typ')

EPC2071_DP = digitizedPlot([]);

EPC2071_DP.plotData = {[0, 3023.8173567033323
4.799578059071743, 2764.8187540101085
7.964135021097055, 2528.004122200026
9.546413502109717, 1982.541204051158
11.128691983122371, 1515.4986560002747
15.03164556962026, 1385.6918628031274
20.094936708860764, 1283.31518555765
30.643459915611817, 1100.6941712522091
35.07383966244726, 1032.497215000986
45.094936708860764, 920.2141114901077
54.27215189873417, 874.312458022074
64.71518987341773, 830.7004475455904
74.52531645569618, 820.1416901494699
88.02742616033754, 809.7171416105634
99.52531645569621, 789.2638692505287
    ]};

EPC2071_DP.logAxes = [0 1];
EPC2071_DP.normAxes = [0 0];
EPC2071_DP.axisLabels = {'Vds, Drain-Source Voltage', 'Coss, Output Capacitance'};
EPC2071_DP.dataLabels = {'Y1'};
EPC2071_DP.SIUnits = { '' , 'p'};


EPC2071_graph = componentPlotData(transistor, EPC2071_DP);
EPC2071.addGraph(EPC2071_graph);

transDB.add(EPC2071)


%% EPC2302
EPC2302 = transistor('EPC2302','en-nMOS','GaN');
EPC2302.manufacturer = 'EPC';
EPC2302.material = 'GaN';
EPC2302.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = 0.15 mA (Continuous)'} )
EPC2302.addParameter('Ids', 101 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2302.addParameter('Idspulse', 408 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2302.addParameter('Ron', 1.4e-3 ,'typ', {'VGS = 5 V, ID = 50 A'} )
EPC2302.addParameter('Ron', 1.8e-3 ,'max', {'VGS = 5 V, ID = 50 A'} )
EPC2302.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 14 mA'} )
EPC2302.addParameter('Vth', 1.3 ,'typ', {'VDS = VGS, ID = 14 mA'} )
EPC2302.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 14 mA'} )
EPC2302.addParameter('Vsd', 1.5 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2302.addParameter('Coss', 1000e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2302.addParameter('Coss', 1200e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2302.addParameter('Rg', 0.5 ,'typ')
EPC2302.addParameter('Qg', 23e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 50 A'} )
EPC2302.addParameter('Qg', 29e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 50 A'} )
EPC2302.addParameter('Length', 5e-3 ,'typ')
EPC2302.addParameter('Width', 3e-3 ,'typ')
EPC2302.addParameter('Height', 0.65e-3 ,'typ')

EPC2302_DP = digitizedPlot([]);

EPC2302_DP.plotData = {[0, 3511.6630244618227
4.485825485859081, 3465.3439149935966
7.480174851526029, 3289.113756496106
10.152687513845448, 2636.489436711749
12.18296084483645, 1929.6173427317497
15.390459096585362, 1608.2013582834913
18.384808462252316, 1526.4162347233355
24.908190872412206, 1357.3043983641735
33.03547336374832, 1130.8567857422427
46.51072480810752, 1005.1147503719706
55.708506821858784, 991.5400825196496
67.4732082643499, 990.7821492429877
84.05074211513285, 989.7151356136741
99.45180582166667, 988.7248752379451

    ]};

EPC2302_DP.logAxes = [0 1];
EPC2302_DP.normAxes = [0 0];
EPC2302_DP.axisLabels = {'Vds, Drain-Source Voltage', 'Coss, Output Capacitance'};
EPC2302_DP.dataLabels = {'Y1'};
EPC2302_DP.SIUnits = { '' , 'p'};


EPC2302_graph = componentPlotData(transistor, EPC2302_DP);
EPC2302.addGraph(EPC2302_graph);


transDB.add(EPC2302)


%% EPC2059
EPC2059 = transistor('EPC2059','en-nMOS','GaN');
EPC2059.manufacturer = 'EPC';
EPC2059.material = 'GaN';
EPC2059.addParameter('Vds', 170 ,'max', {'VGS = 0 V, ID = 0.15 mA (Continuous)'} )
EPC2059.addParameter('Ids', 24 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2059.addParameter('Idspulse', 102 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2059.addParameter('Ron', 6.8e-3 ,'typ', {'VGS = 5 V, ID = 10 A'} )
EPC2059.addParameter('Ron', 9e-3 ,'max', {'VGS = 5 V, ID = 10 A'} )
EPC2059.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 3 mA'} )
EPC2059.addParameter('Vth', 1 ,'typ', {'VDS = VGS, ID = 3 mA'} )
EPC2059.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 3 mA'} )
EPC2059.addParameter('Vsd', 1.6 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2059.addParameter('Coss', 267e-12 ,'typ', {'VDS = 85 V, VGS = 0 V'} )
EPC2059.addParameter('Coss', 401e-12 ,'max', {'VDS = 85 V, VGS = 0 V'} )
EPC2059.addParameter('Rg', 0.5 ,'typ')
EPC2059.addParameter('Qg', 5.7e-9 ,'typ', {'VDS = 85 V, VGS = 5 V, ID = 10 A'} )
EPC2059.addParameter('Qg', 7.4e-9 ,'max', {'VDS = 85 V, VGS = 5 V, ID = 10 A'} )
EPC2059.addParameter('Length', 2.8e-3 ,'typ')
EPC2059.addParameter('Width', 1.4e-3 ,'typ')
EPC2059.addParameter('Height', 0.638e-3 ,'typ')

transDB.add(EPC2059)


%% EPC2054
EPC2054 = transistor('EPC2054','en-nMOS','GaN');
EPC2054.manufacturer = 'EPC';
EPC2054.material = 'GaN';
EPC2054.addParameter('Vds', 200 ,'max', {'VGS = 0 V, ID = 0.12 mA (Continuous)'} )
EPC2054.addParameter('Ids', 3 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2054.addParameter('Idspulse', 32 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2054.addParameter('Ron', 32e-3 ,'typ', {'VGS = 5 V, ID = 1 A'} )
EPC2054.addParameter('Ron', 43e-3 ,'max', {'VGS = 5 V, ID = 1 A'} )
EPC2054.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 1 mA'} )
EPC2054.addParameter('Vth', 1.2 ,'typ', {'VDS = VGS, ID = 1 mA'} )
EPC2054.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 1 mA'} )
EPC2054.addParameter('Vsd', 1.5 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2054.addParameter('Coss', 89e-12 ,'typ', {'VDS = 100 V, VGS = 0 V'} )
EPC2054.addParameter('Coss', 134e-12 ,'max', {'VDS = 100 V, VGS = 0 V'} )
EPC2054.addParameter('Rg', 0.8 ,'typ')
EPC2054.addParameter('Qg', 2.9e-9 ,'typ', {'VDS = 100 V, VGS = 5 V, ID = 1 A'} )
EPC2054.addParameter('Qg', 4.3e-9 ,'max', {'VDS = 100 V, VGS = 5 V, ID = 1 A'} )
EPC2054.addParameter('Length', 1.3e-3 ,'typ')
EPC2054.addParameter('Width', 1.3e-3 ,'typ')
EPC2054.addParameter('Height', 0.638e-3 ,'typ')

transDB.add(EPC2054)


%% EPC2207
EPC2207 = transistor('EPC2207','en-nMOS','GaN');
EPC2207.manufacturer = 'EPC';
EPC2207.material = 'GaN';
EPC2207.addParameter('Vds', 200 ,'max', {'VGS = 0 V, ID = 0.2 mA (Continuous)'} )
EPC2207.addParameter('Ids', 14 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2207.addParameter('Idspulse', 54 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2207.addParameter('Ron', 15e-3 ,'typ', {'VGS = 5 V, ID = 14 A'} )
EPC2207.addParameter('Ron', 22e-3 ,'max', {'VGS = 5 V, ID = 14 A'} )
EPC2207.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 2 mA'} )
EPC2207.addParameter('Vth', 1.1 ,'typ', {'VDS = VGS, ID = 2 mA'} )
EPC2207.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 2 mA'} )
EPC2207.addParameter('Vsd', 1.7 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2207.addParameter('Coss', 130e-12 ,'typ', {'VDS = 100 V, VGS = 0 V'} )
EPC2207.addParameter('Coss', 195e-12 ,'max', {'VDS = 100 V, VGS = 0 V'} )
EPC2207.addParameter('Rg', 0.3 ,'typ')
EPC2207.addParameter('Qg', 4.5e-9 ,'typ', {'VDS = 100 V, VGS = 5 V, ID = 14 A'} )
EPC2207.addParameter('Qg', 5.9e-9 ,'max', {'VDS = 100 V, VGS = 5 V, ID = 14 A'} )
EPC2207.addParameter('Length', 2.8e-3 ,'typ')
EPC2207.addParameter('Width', 0.925e-3 ,'typ')
EPC2207.addParameter('Height', 0.638e-3 ,'typ')


EPC2207_DP = digitizedPlot([]);

EPC2207_DP.plotData = {[0, 434.7877000443586
5.89942076153889, 410.29945517634394
13.543484744627476, 376.1200021449432
18.480961525380188, 328.46792343184757
26.60398821639549, 260.3946207408678
44.43997007006316, 216.7223393990962
62.9128107454752, 182.13279291281492
77.24502257826524, 162.20553434076626
96.99199949767942, 132.41745496767834
108.13891258051515, 126.19203819072058
124.5406853536913, 119.1153772551344
137.91664617268316, 116.86818097316868
159.73201195104465, 115.8029384361192
182.98049844856192, 114.75138718350468
198.90414571483893, 113.68921989869752
    ]};

EPC2207_DP.logAxes = [0 1];
EPC2207_DP.normAxes = [0 0];
EPC2207_DP.axisLabels = {'Vds, Drain-Source Voltage', 'Coss, Output Capacitance'};
EPC2207_DP.dataLabels = {'Y1'};
EPC2207_DP.SIUnits = { '' , 'p'};


EPC2207_graph = componentPlotData(transistor, EPC2207_DP);
EPC2207.addGraph(EPC2207_graph);


transDB.add(EPC2207)


%% EPC2215
EPC2215 = transistor('EPC2215','en-nMOS','GaN');
EPC2215.manufacturer = 'EPC';
EPC2215.material = 'GaN';
EPC2215.addParameter('Vds', 200 ,'max', {'VGS = 0 V, ID = 0.6 mA (Continuous)'} )
EPC2215.addParameter('Ids', 32 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2215.addParameter('Idspulse', 162 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2215.addParameter('Ron', 6e-3 ,'typ', {'VGS = 5 V, ID = 20 A'} )
EPC2215.addParameter('Ron', 8e-3 ,'max', {'VGS = 5 V, ID = 20 A'} )
EPC2215.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 6 mA'} )
EPC2215.addParameter('Vth', 1.1 ,'typ', {'VDS = VGS, ID = 6 mA'} )
EPC2215.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 6 mA'} )
EPC2215.addParameter('Vsd', 1.7 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2215.addParameter('Coss', 390e-12 ,'typ', {'VDS = 100 V, VGS = 0 V'} )
EPC2215.addParameter('Coss', 585e-12 ,'max', {'VDS = 100 V, VGS = 0 V'} )
EPC2215.addParameter('Rg', 0.4 ,'typ')
EPC2215.addParameter('Qg', 13.6e-9 ,'typ', {'VDS = 100 V, VGS = 5 V, ID = 20 A'} )
EPC2215.addParameter('Qg', 17.7e-9 ,'max', {'VDS = 100 V, VGS = 5 V, ID = 20 A'} )
EPC2215.addParameter('Length', 4.6e-3 ,'typ')
EPC2215.addParameter('Width', 1.6e-3 ,'typ')
EPC2215.addParameter('Height', 0.638e-3 ,'typ')

transDB.add(EPC2215)



%% EPC2050
EPC2050 = transistor('EPC2050','en-nMOS','GaN');
EPC2050.manufacturer = 'EPC';
EPC2050.material = 'GaN';
EPC2050.addParameter('Vds', 350 ,'max', {'VGS = 0 V, ID = 0.08 mA (Continuous)'} )
EPC2050.addParameter('Ids', 6.3 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2050.addParameter('Idspulse', 26 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2050.addParameter('Ron', 55e-3 ,'typ', {'VGS = 5 V, ID = 6 A'} )
EPC2050.addParameter('Ron', 80e-3 ,'max', {'VGS = 5 V, ID = 6 A'} )
EPC2050.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 1 mA'} )
EPC2050.addParameter('Vth', 1.3 ,'typ', {'VDS = VGS, ID = 1 mA'} )
EPC2050.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 1 mA'} )
EPC2050.addParameter('Vsd', 1.5 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2050.addParameter('Coss', 81e-12 ,'typ', {'VDS = 280 V, VGS = 0 V'} )
EPC2050.addParameter('Coss', 122e-12 ,'max', {'VDS = 280 V, VGS = 0 V'} )
EPC2050.addParameter('Rg', 0.5 ,'typ')
EPC2050.addParameter('Qg', 2.9e-9 ,'typ', {'VDS = 280 V, VGS = 5 V, ID = 6 A'} )
EPC2050.addParameter('Qg', 4e-9 ,'max', {'VDS = 280 V, VGS = 5 V, ID = 6 A'} )
EPC2050.addParameter('Length', 1.95e-3 ,'typ')
EPC2050.addParameter('Width', 1.95e-3 ,'typ')
EPC2050.addParameter('Height', 0.638e-3 ,'typ')

EPC2050_DP = digitizedPlot([]);

EPC2050_DP.plotData = {[0, 402.39636625246595
7.566309242722563, 385.56427533465546
10.351856421807248, 324.8603375812928
10.96277264673862, 245.90678845876414
13.435649007042121, 192.22514839542342
21.808158498190807, 184.187937586689
34.212570819760046, 176.49766132347892
40.72095636364188, 153.57378354419112
48.469566758491275, 132.20634252350297
62.11492520165944, 128.05291497978007
79.17081182615931, 120.11460585940976
99.32852280168302, 113.88694543199583
115.76483795367285, 109.13863182282357
137.78342715862135, 103.48299607526283
163.21470835679602, 101.32944106940626
179.96189115098815, 99.20752119962266
202.29098736282324, 95.07981179727312
239.50710941005732, 91.14460044669862
270.82936849149854, 84.60205791658454
296.2566827011997, 73.635080037181
320.7553609727674, 67.61425120543247
348.97759887976366, 65.50467218894435
    ]};

EPC2050_DP.logAxes = [0 1];
EPC2050_DP.normAxes = [0 0];
EPC2050_DP.axisLabels = {'Vds, Drain-Source Voltage', 'Coss, Output Capacitance'};
EPC2050_DP.dataLabels = {'Y1'};
EPC2050_DP.SIUnits = { '' , 'p'};


EPC2050_graph = componentPlotData(transistor, EPC2050_DP);
EPC2050.addGraph(EPC2050_graph);



transDB.add(EPC2050)



%% EPC2306
EPC2306 = transistor('EPC2306','en-nMOS','GaN');
EPC2306.manufacturer = 'EPC';
EPC2306.material = 'GaN';
EPC2306.addParameter('Vds', 100 ,'max', {'VGS = 0 V, ID = [TBD] (Continuous)'} )
EPC2306.addParameter('Ids', 48 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2306.addParameter('Idspulse', 197 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2306.addParameter('Ron', 3e-3 ,'typ', {'VGS = 5 V, ID = 25 A'} )
EPC2306.addParameter('Ron', 3.8e-3 ,'max', {'VGS = 5 V, ID = 25 A'} )
EPC2306.addParameter('Vth', 0.8 ,'min', {'VDS = VGS, ID = 7 mA'} )
EPC2306.addParameter('Vth', 1.3 ,'typ', {'VDS = VGS, ID = 7 mA'} )
EPC2306.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 7 mA'} )
EPC2306.addParameter('Vsd', 1.6 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2306.addParameter('Coss', 482e-12 ,'typ', {'VDS = 50 V, VGS = 0 V'} )
EPC2306.addParameter('Coss', 559e-12 ,'max', {'VDS = 50 V, VGS = 0 V'} )
EPC2306.addParameter('Rg', 0.4 ,'typ')
EPC2306.addParameter('Qg', 11.6e-9 ,'typ', {'VDS = 50 V, VGS = 5 V, ID = 25 A'} )
EPC2306.addParameter('Qg', 16.3e-9 ,'max', {'VDS = 50 V, VGS = 5 V, ID = 25 A'} )
EPC2306.addParameter('Length', 5.0e-3 ,'typ')
EPC2306.addParameter('Width', 3.0e-3 ,'typ')
EPC2306.addParameter('Height', 0.65e-3 ,'typ')

EPC2306_DP = digitizedPlot([]);

EPC2306_DP.plotData = {[0, 1591.7229838937205
4.405874499332441, 1523.6728989923452
6.364040943480193, 1474.739467821668
8.945260347129496, 1267.2620214505598
11.170449488206488, 1020.6227322556088
13.751668891855806, 839.9182364708204
16.599910992434353, 778.465626391075
20.605251446372936, 729.2599822850665
23.275478415665315, 698.1965155421425
25.50066755674232, 654.1904411004236
28.704939919893185, 568.2257294614107
31.464174454828647, 526.656566606708
34.57943925233643, 504.198912085617
40.36493101913662, 487.8044725043856
44.28126390743212, 482.35532444090387
51.134846461949245, 476.8156229384583
55.85224744103247, 476.5725264853143
63.68491321762347, 476.16916834064966
70.36048064085443, 475.8256666732336
77.03604806408542, 470.3700079323163
84.24566088117487, 470.0035548147243
99.64396973742765, 474.32176463924367

    ]};

EPC2306_DP.logAxes = [0 1];
EPC2306_DP.normAxes = [0 0];
EPC2306_DP.axisLabels = {'Vds, Drain-Source Voltage', 'Coss, Output Capacitance'};
EPC2306_DP.dataLabels = {'Y1'};
EPC2306_DP.SIUnits = { '' , 'p'};


EPC2306_graph = componentPlotData(transistor, EPC2306_DP);
EPC2306.addGraph(EPC2306_graph);



transDB.add(EPC2306)


%% EPC2308
EPC2308 = transistor('EPC2308','en-nMOS','GaN');
EPC2308.manufacturer = 'EPC';
EPC2308.material = 'GaN';
EPC2308.addParameter('Vds', 180 ,'max', {'VGS = 0 V, ID = [TBD] (Continuous)'} )
EPC2308.addParameter('Ids', 48 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2308.addParameter('Idspulse', 157 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2308.addParameter('Ron', 4.6e-3 ,'typ', {'VGS = 5 V, ID = 15 A'} )
EPC2308.addParameter('Ron', 6e-3 ,'max', {'VGS = 5 V, ID = 15 A'} )
EPC2308.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 5 mA'} )
EPC2308.addParameter('Vth', 1.2 ,'typ', {'VDS = VGS, ID = 5 mA'} )
EPC2308.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 5 mA'} )
EPC2308.addParameter('Vsd', 1.5 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2308.addParameter('Coss', 405e-12 ,'typ', {'VDS = 75 V, VGS = 0 V'} )
EPC2308.addParameter('Coss', 592e-12 ,'max', {'VDS = 75 V, VGS = 0 V'} )
EPC2308.addParameter('Rg', 0.4 ,'typ')
EPC2308.addParameter('Qg', 10.6e-9 ,'typ', {'VDS = 75 V, VGS = 5 V, ID = 15 A'} )
EPC2308.addParameter('Qg', 13.8e-9 ,'max', {'VDS = 75 V, VGS = 5 V, ID = 15 A'} )
EPC2308.addParameter('Length', 5.0e-3 ,'typ')
EPC2308.addParameter('Width', 3.0e-3 ,'typ')
EPC2308.addParameter('Height', 0.65e-3 ,'typ')



transDB.add(EPC2308)


%% EPC2305
EPC2305 = transistor('EPC2305','en-nMOS','GaN');
EPC2305.manufacturer = 'EPC';
EPC2305.material = 'GaN';
EPC2305.addParameter('Vds', 180 ,'max', {'VGS = 0 V, ID = [TBD] (Continuous)'} )
EPC2305.addParameter('Ids', 80 ,'max', {'(Continuous) TA = 25°C)'} )
EPC2305.addParameter('Idspulse', 329 ,'max', {'(Pulsed) 25°C Tpulse = 300μs)'} )
EPC2305.addParameter('Ron', 2.2e-3 ,'typ', {'VGS = 5 V, ID = 30 A'} )
EPC2305.addParameter('Vth', 0.7 ,'min', {'VDS = VGS, ID = 11 mA'} )
EPC2305.addParameter('Vth', 1.1 ,'typ', {'VDS = VGS, ID = 11 mA'} )
EPC2305.addParameter('Vth', 2.5 ,'max', {'VDS = VGS, ID = 11 mA'} )
EPC2305.addParameter('Vsd', 1.4 ,'typ', {'IS = 0.5 A, VGS = 0 V '} )
EPC2305.addParameter('Coss', 920e-12 ,'typ', {'VDS = 75 V, VGS = 0 V'} )
EPC2305.addParameter('Rg', 0.5 ,'typ')
EPC2305.addParameter('Qg', 21e-9 ,'typ', {'VDS = 75 V, VGS = 5 V, ID = 30 A'} )
EPC2305.addParameter('Length', 5.0e-3 ,'typ')
EPC2305.addParameter('Width', 3.0e-3 ,'typ')
EPC2305.addParameter('Height', 0.65e-3 ,'typ')


EPC2305_DP = digitizedPlot([]);

EPC2305_DP.plotData = {[0, 2932.488966344632
5.844155844155807, 2796.2252837656283
9.268004722550192, 2590.410066059327
10.68476977567887, 2287.165181240119
12.337662337662321, 1870.3998017823067
13.754427390791038, 1589.3302008252851
16.942148760330596, 1430.595819964196
20.247933884297463, 1338.050076580898
24.970484061393094, 1239.6662150584023
30.04722550177093, 1159.602980500648
35.950413223140465, 1095.2143974753365
43.27036599763868, 1064.6706257085516
54.25029515938602, 1025.3464427917693
64.04958677685948, 977.9820726828559
76.32821723730814, 915.2406390127885
89.78748524203066, 832.3106774948415
105.13577331759147, 756.986361150932
118.12278630460449, 729.1201592700078
134.76977567886658, 729.8960969282458
147.8748524203069, 730.5075223352148
    ]};

EPC2305_DP.logAxes = [0 1];
EPC2305_DP.normAxes = [0 0];
EPC2305_DP.axisLabels = {'Vds, Drain-Source Voltage', 'Coss, Output Capacitance'};
EPC2305_DP.dataLabels = {'Y1'};
EPC2305_DP.SIUnits = { '' , 'p'};


EPC2305_graph = componentPlotData(transistor, EPC2305_DP);
EPC2305.addGraph(EPC2305_graph);


transDB.add(EPC2305)

%transDB(end).datasheet


