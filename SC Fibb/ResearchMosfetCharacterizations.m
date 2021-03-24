% Different Mosfet characterization
clear
clc
close all
%% 8HVnLDMOS nbl #1
%fingers =2

L1=800*10^-9;
W1=[60 125 250 500 999].*(400*10^-6);

Ron1=[57.7954 27.7375 13.8648 6.9288 3.4649].*10^-3; %at fingers =4 Ron(W=400)=Ron(W=200)
aw1=0.001385;
b1=-1;
aa1=1.102*10^-9;
Area1=L1*W1;
for i=1:10000
    w(i)=i*10^-4+24*10^-3;
    ar1(i)=w(i)*L1;
    armm1(i)=ar1(i)*10^6;
Roncalcw1(i)=aw1*w(i)^b1;
Roncalca1(i)=aa1*ar1(i)^b1;
end
%% #1 Coss
vds1_m60=[56.12E-3
66.41E-3
79.07E-3
83.04E-3
91.04E-3
107.3E-3
140.5E-3
186.5E-3
259.0E-3
387.1E-3
448.1E-3
576.3E-3
855.1E-3
1.488
2.650
4.379
6.798];
vds1_m125=[50.09E-3
65.02E-3
85.22E-3
116.2E-3
169.0E-3
193.5E-3
244.3E-3
352.6E-3
410.0E-3
530.6E-3
793.1E-3
1.389
2.563
4.236
6.593];
vds1_m250=[56.33E-3
80.80E-3
91.95E-3
114.7E-3
161.8E-3
186.2E-3
236.7E-3
344.5E-3
584.6E-3
1.145
2.159
3.599
5.690];
vds1_m500=[55.59E-3
77.44E-3
88.57E-3
111.2E-3
158.2E-3
259.0E-3
486.6E-3
1.027
2.002
3.342
5.311];
vds1_m999=[54.00E-3
75.85E-3
121.1E-3
218.5E-3
439.9E-3
971.2E-3
1.946
3.247
5.163];

% armm1_20=Area1(1).*ones(length(vds),1);
% area1_m60=L1*400*60.*ones(length(vds),1);
% area1_m125=L1*400*125.*ones(length(vds),1);
% area1_m250=L1*400*250.*ones(length(vds),1);
% area1_m500=L1*400*500.*ones(length(vds),1);
% area1_m999=L1*400*999.*ones(length(vds),1);

c1_m60=8.876*10^-11;
d1_m60=-0.2564;

coss1_m60=[153.0E-12
151.5E-12
149.6E-12
148.2E-12
147.3E-12
145.4E-12
141.6E-12
136.0E-12
128.7E-12
118.9E-12
111.7E-12
106.3E-12
97.75E-12
86.10E-12
73.75E-12
63.38E-12
55.31E-12];
coss1_m125=[319.5E-12
315.6E-12
309.6E-12
301.2E-12
288.4E-12
277.6E-12
268.2E-12
251.5E-12
237.5E-12
225.9E-12
207.7E-12
182.8E-12
155.7E-12
133.3E-12
116.3E-12];
coss1_m250=[636.1E-12
622.6E-12
610.7E-12
599.7E-12
578.4E-12
558.2E-12
539.3E-12
505.6E-12
454.0E-12
389.1E-12
327.3E-12
279.7E-12
243.4E-12];
coss1_m500=[1.269E-9
1.247E-9
1.225E-9
1.203E-9
1.160E-9
1.081E-9
957.7E-12
806.2E-12
670.6E-12
571.4E-12
496.9E-12];

coss1_m999=[2.540E-9
2.494E-9
2.406E-9
2.239E-9
1.969E-9
1.641E-9
1.354E-9
1.151E-9
1.001E-9];
c1_m60=1.525E-10;
d1_m60=-.1275;
e1_m60=-5.901E-11;
for i=1:5000
    vds=i*.001;
coss1_m60eq(i)=c1_m60*vds^d1_m60+e1_m60; 
end
c1_m250=6.247E-10;
d1_m250=-.1302;
e1_m250=-2.329E-10;
for i=1:5000
    vds=i*.001;
coss1_m250eq(i)=c1_m250*vds^d1_m250+e1_m250; 
end
c1_m125=3.137E-10;
d1_m125=-.1283;
e1_m125=-1.184E-10;
for i=1:5000
    vds=i*.001;
coss1_m125eq(i)=c1_m125*vds^d1_m125+e1_m125; 
end
c1_m500=1.273E-09;
d1_m500=-.1259;
e1_m500=-4.891E-10;
for i=1:5000
    vds=i*.001;
coss1_m500eq(i)=c1_m500*vds^d1_m500+e1_m500; 
end
c1_m999=2.409E-09;
d1_m999=-.1321;
e1_m999=-8.386E-10;
for i=1:5000
    vds=i*.001;
coss1_m999eq(i)=c1_m999*vds^d1_m999+e1_m999; 
end
coss1_vds_4=[coss1_m60eq(4000) coss1_m125eq(4000) coss1_m250eq(4000) coss1_m500eq(4000) coss1_m999eq(4000)];


%% #1 Qg
vgs1=[250.0E-3
500.0E-3
750.0E-3
1.000
1.250
1.500
1.750
2.000
2.250
2.500
2.750
3.000
3.250
3.500
3.750
4.000
4.250
4.500
4.750
5.000];
Qg1_m60=[5.978E-12
11.84E-12
158.5E-12
194.2E-12
222.2E-12
248.3E-12
274.1E-12
299.9E-12
325.7E-12
351.4E-12
377.2E-12
403.0E-12
428.8E-12
454.7E-12
480.5E-12
506.3E-12
532.2E-12
558.1E-12
583.9E-12
609.7E-12];
Qg1_m125=[12.48E-12
24.71E-12
334.8E-12
408.0E-12
464.3E-12
518.1E-12
571.8E-12
625.4E-12
679.0E-12
732.7E-12
786.4E-12
840.2E-12
894.0E-12
947.8E-12
1.002E-9
1.055E-9
1.109E-9
1.163E-9
1.217E-9
1.271E-9];
Qg1_m250=[24.97E-12
49.54E-12
677.7E-12
820.1E-12
930.0E-12
1.037E-9
1.144E-9
1.251E-9
1.359E-9
1.466E-9
1.573E-9
1.681E-9
1.788E-9
1.896E-9
2.004E-9
2.111E-9
2.219E-9
2.327E-9
2.435E-9
2.543E-9];
Qg1_m500=[49.96E-12
99.42E-12
1.372E-9
1.645E-9
1.861E-9
2.075E-9
2.289E-9
2.503E-9
2.718E-9
2.932E-9
3.147E-9
3.362E-9
3.577E-9
3.793E-9
4.008E-9
4.223E-9
4.439E-9
4.654E-9
4.870E-9
5.086E-9];
Qg1_m999=[99.84E-12
200.7E-12
2.777E-9
3.292E-9
3.720E-9
4.147E-9
4.574E-9
5.002E-9
5.431E-9
5.859E-9
6.289E-9
6.718E-9
7.148E-9
7.578E-9
8.009E-9
8.439E-9
8.870E-9
9.300E-9
9.731E-9
10.16E-9];
Qg1_vgs_5=[Qg1_m60(20,1) Qg1_m125(20,1) Qg1_m250(20,1) Qg1_m500(20,1) Qg1_m999(20,1)];
%% 8HVnLDMOS iso #2
%fingers =2
L2=900*10^-9;
W2=[60 125 250 500 999].*(400*10^-6);
Ron2=[60.9994 29.2749 14.6331 7.31258 3.65664]*10^-3;
a2=0.001462;
b2=-1;
aa2=1.308*10^-9;
Area2=L2*W2;
for i=1:10000
    w(i)=i*10^-4+24*10^-3;
    ar2(i)=w(i)*L2;
    armm2(i)=ar2(i)*10^6;
Roncalcw2(i)=a2*w(i)^b2;
Roncalca2(i)=aa2*ar2(i)^b2;
end
%% #2 Coss 
vds2_m60=[61.39E-3
105.2E-3
153.0E-3
177.5E-3
190.0E-3
215.1E-3
266.7E-3
343.1E-3
448.3E-3
557.6E-3
679.4E-3
809.6E-3
959.6E-3
1.134
1.332
1.528
1.702
1.746
1.834
1.972
2.013
2.097
2.140
2.225
2.399
2.460
2.584
2.646
2.772
2.836
2.965
3.228
3.742
4.446
5.549];
vds2_m125=[84.63E-3
175.2E-3
374.1E-3
834.1E-3
1.639
2.736
4.269
6.332];
vds2_m250=[86.57E-3
180.1E-3
386.1E-3
863.5E-3
1.706
2.826
4.393
6.486];
vds2_m500=[86.87E-3
180.5E-3
386.5E-3
863.9E-3
1.702
2.821
4.387
6.484];
vds2_m999=[87.11E-3
180.8E-3
387.1E-3
865.2E-3
1.709
2.831
4.400
6.500];
% area2_m60=L2*400*60.*ones(length(vds),1);
% area2_m125=L2*400*125.*ones(length(vds),1);
% area2_m250=L2*400*250.*ones(length(vds),1);
% area2_m500=L2*400*500.*ones(length(vds),1);
% area2_m999=L2*400*999.*ones(length(vds),1);
coss2_m60=[134.5E-12
128.3E-12
122.3E-12
117.9E-12
116.0E-12
113.6E-12
109.3E-12
103.2E-12
96.42E-12
90.32E-12
85.04E-12
80.38E-12
76.08E-12
71.97E-12
68.07E-12
64.69E-12
62.00E-12
60.77E-12
59.87E-12
58.49E-12
57.62E-12
56.91E-12
56.32E-12
55.66E-12
54.37E-12
53.49E-12
52.67E-12
52.01E-12
51.25E-12
50.64E-12
49.93E-12
48.61E-12
46.30E-12
43.45E-12
40.09E-12];
coss2_m125=[272.6E-12
253.1E-12
216.0E-12
175.7E-12
139.1E-12
114.4E-12
94.98E-12
81.53E-12];
coss2_m250=[547.7E-12
500.4E-12
431.7E-12
345.3E-12
276.7E-12
223.9E-12
189.4E-12
160.0E-12];
coss2_m500=[1.089E-9
1.005E-9
858.3E-12
694.5E-12
550.3E-12
451.1E-12
376.6E-12
322.3E-12];
coss2_m999=[2.181E-9
2.002E-9
1.719E-9
1.383E-9
1.102E-9
897.0E-12
754.1E-12
640.9E-12];
c2_m60=9.756E-11;
d2_m60=-.2125;
e2_m60=-2.532E-11;
for i=1:5000
    vds=i*.001;
coss2_m60eq(i)=c2_m60*vds^d2_m60+e2_m60; 
end
c2_m250=3.89E-10;
d2_m250=-.2177;
e2_m250=-7.875E-11;
for i=1:5000
    vds=i*.001;
coss2_m250eq(i)=c2_m250*vds^d2_m250+e2_m250; 
end
c2_m125=1.997E-10;
d2_m125=-.2109;
e2_m125=-4.428E-11;
for i=1:5000
    vds=i*.001;
coss2_m125eq(i)=c2_m125*vds^d2_m125+e2_m125; 
end
c2_m500=7.889E-10;
d2_m500=-.2144;
e2_m500=-1.68E-10;
for i=1:5000
    vds=i*.001;
coss2_m500eq(i)=c2_m500*vds^d2_m500+e2_m500; 
end
c2_m999=1.564E-09;
d2_m999=-.2164;
e2_m999=-3.236E-10;
for i=1:5000
    vds=i*.001;
coss2_m999eq(i)=c2_m999*vds^d2_m999+e2_m999; 
end
coss2_vds_4=[coss2_m60eq(4000) coss2_m125eq(4000) coss2_m250eq(4000) coss2_m500eq(4000) coss2_m999eq(4000)];

%% #2 Qg
vgs2=[250.0E-3
500.0E-3
750.0E-3
1.000
1.250
1.500
1.750
2.000
2.250
2.500
2.750
3.000
3.250
3.500
3.750
4.000
4.250
4.500
4.750
5.000];
Qg2_m60=[9.230E-12
17.23E-12
26.08E-12
201.3E-12
236.6E-12
266.0E-12
294.7E-12
323.3E-12
351.9E-12
380.5E-12
409.1E-12
437.8E-12
466.5E-12
495.2E-12
524.0E-12
552.6E-12
581.4E-12
610.2E-12
639.1E-12
667.8E-12];
Qg2_m125=[19.25E-12
35.94E-12
54.59E-12
424.0E-12
495.1E-12
555.4E-12
614.8E-12
674.2E-12
733.7E-12
793.3E-12
852.9E-12
912.7E-12
972.5E-12
1.032E-9
1.092E-9
1.152E-9
1.212E-9
1.272E-9
1.332E-9
1.392E-9];
Qg2_m250=[38.52E-12
71.92E-12
115.6E-12
857.1E-12
992.7E-12
1.112E-9
1.230E-9
1.349E-9
1.468E-9
1.587E-9
1.706E-9
1.826E-9
1.945E-9
2.065E-9
2.185E-9
2.305E-9
2.425E-9
2.544E-9
2.665E-9
2.784E-9];
Qg2_m500=[77.10E-12
143.9E-12
287.6E-12
1.731E-9
1.988E-9
2.225E-9
2.462E-9
2.698E-9
2.937E-9
3.175E-9
3.413E-9
3.652E-9
3.891E-9
4.131E-9
4.370E-9
4.610E-9
4.850E-9
5.090E-9
5.329E-9
5.570E-9];
Qg2_m999=[154.0E-12
287.6E-12
627.4E-12
3.483E-9
3.975E-9
4.446E-9
4.919E-9
5.393E-9
5.868E-9
6.344E-9
6.820E-9
7.298E-9
7.775E-9
8.254E-9
8.732E-9
9.211E-9
9.690E-9
10.17E-9
10.65E-9
11.13E-9];
Qg2_vgs_5=[Qg2_m60(20,1) Qg2_m125(20,1) Qg2_m250(20,1) Qg2_m500(20,1) Qg2_m999(20,1)];
%% 12HVnLDMOS iso hp mac #3
%fingers =2
L3=900*10^-9;
W3=[60 125 250 500 999].*(400*10^-6);
Ron3=[60.4805 29.0297 14.5124 7.25405 3.62885]*10^-3;
a3=.001451;
b3=-1;
aa3=1.302*10^-9;
Area3=L3*W3;
for i=1:10000
    w(i)=i*10^-4+24*10^-3;
    ar3(i)=w(i)*L3;
    armm3(i)=ar3(i)*10^6;
Roncalcw3(i)=a3*w(i)^b3;
Roncalca3(i)=aa3*ar3(i)^b3;
end
%% #3 Coss 
vds3_m60=[56.66E-3
105.8E-3
209.1E-3
435.4E-3
956.2E-3
1.837
3.001
4.584
6.675];
vds3_m125=[54.29E-3
105.7E-3
214.1E-3
452.4E-3
1.003
1.904
3.093
4.706
6.835];
vds3_m250=[51.83E-3
103.2E-3
211.4E-3
449.4E-3
999.5E-3
1.904
3.090
4.704
6.832];
vds3_m500=[67.83E-3
137.5E-3
286.7E-3
621.4E-3
1.412
2.425
3.795
5.640];
vds3_m999=[67.11E-3
136.8E-3
286.1E-3
621.1E-3
1.412
2.425
3.796
5.641];
% area3_m60=L3*400*60.*ones(length(vds),1);
% area3_m125=L3*400*125.*ones(length(vds),1);
% area3_m250=L3*400*250.*ones(length(vds),1);
% area3_m500=L3*400*500.*ones(length(vds),1);
% area3_m999=L3*400*999.*ones(length(vds),1);
coss3_m60=[132.6E-12
128.3E-12
119.5E-12
105.2E-12
87.42E-12
72.32E-12
61.07E-12
53.16E-12
46.75E-12];
coss3_m125=[277.1E-12
266.9E-12
248.3E-12
217.0E-12
180.3E-12
148.7E-12
126.2E-12
109.6E-12
96.81E-12];
coss3_m250=[554.3E-12
535.5E-12
496.8E-12
435.3E-12
360.4E-12
298.0E-12
252.0E-12
219.6E-12
193.4E-12];
coss3_m500=[1.098E-9
1.044E-9
949.1E-12
809.2E-12
659.0E-12
547.5E-12
471.1E-12
412.9E-12];
coss3_m999=[2.191E-9
2.090E-9
1.894E-9
1.619E-9
1.315E-9
1.096E-9
939.8E-12
826.2E-12];
c3_m60=1.173E-10;
d3_m60=-.1487;
e3_m60=-3.608E-11;
for i=1:5000
    vds=i*.001;
coss3_m60eq(i)=c3_m60*vds^d3_m60+e3_m60; 
end
c3_m250=4.889E-10;
d3_m250=-.1464;
e3_m250=-1.513E-10;
for i=1:5000
    vds=i*.001;
coss3_m250eq(i)=c3_m250*vds^d3_m250+e3_m250; 
end
c3_m125=2.418E-10;
d3_m125=-.1492;
e3_m125=-7.288E-11;
for i=1:5000
    vds=i*.001;
coss3_m125eq(i)=c3_m125*vds^d3_m125+e3_m125; 
end
c3_m500=9.374E-10;
d3_m500=-.1578;
e3_m500=-2.566E-10;
for i=1:5000
    vds=i*.001;
coss3_m500eq(i)=c3_m500*vds^d3_m500+e3_m500; 
end
c3_m999=1.888E-09;
d3_m999=-.1562;
e3_m999=-5.27E-10;
for i=1:5000
    vds=i*.001;
coss3_m999eq(i)=c3_m999*vds^d3_m999+e3_m999; 
end
coss3_vds_4=[coss3_m60eq(4000) coss3_m125eq(4000) coss3_m250eq(4000) coss3_m500eq(4000) coss3_m999eq(4000)];

%% #3 Qg
vgs3=[250.0E-3
500.0E-3
750.0E-3
1.000
1.250
1.500
1.750
2.000
2.250
2.500
2.750
3.000
3.250
3.500
3.750
4.000
4.250
4.500
4.750
5.000];
Qg3_m60=[9.952E-12
18.54E-12
26.11E-12
32.83E-12
229.9E-12
256.2E-12
279.6E-12
302.6E-12
325.6E-12
348.6E-12
371.5E-12
394.5E-12
417.5E-12
440.5E-12
463.5E-12
486.5E-12
509.6E-12
532.6E-12
555.7E-12
578.7E-12];
Qg3_m125=[20.73E-12
38.63E-12
54.41E-12
109.5E-12
482.6E-12
535.1E-12
583.2E-12
630.9E-12
678.7E-12
726.5E-12
774.2E-12
822.1E-12
870.0E-12
917.9E-12
965.8E-12
1.014E-9
1.062E-9
1.110E-9
1.158E-9
1.206E-9];
Qg3_m250=[41.47E-12
77.26E-12
108.8E-12
272.7E-12
971.0E-12
1.072E-9
1.167E-9
1.262E-9
1.358E-9
1.453E-9
1.549E-9
1.644E-9
1.740E-9
1.836E-9
1.932E-9
2.028E-9
2.124E-9
2.220E-9
2.316E-9
2.412E-9];
Qg3_m500=[82.93E-12
154.5E-12
217.7E-12
592.6E-12
1.950E-9
2.145E-9
2.335E-9
2.525E-9
2.716E-9
2.907E-9
3.098E-9
3.289E-9
3.481E-9
3.672E-9
3.864E-9
4.056E-9
4.247E-9
4.439E-9
4.631E-9
4.823E-9];
Qg3_m999=[165.7E-12
308.7E-12
435.1E-12
1.206E-9
3.907E-9
4.287E-9
4.665E-9
5.045E-9
5.426E-9
5.808E-9
6.190E-9
6.572E-9
6.954E-9
7.337E-9
7.720E-9
8.103E-9
8.487E-9
8.870E-9
9.254E-9
9.637E-9];
Qg3_vgs_5=[Qg3_m60(20,1) Qg3_m125(20,1) Qg3_m250(20,1) Qg3_m500(20,1) Qg3_m999(20,1)];
%% 12HVnLDMOS iso mac #4
%fingers =2
% Ron
L4=1*10^-6;
W4=[60 125 250 500 999].*(400*10^-6);
Ron4=[69.112 33.1725 16.5851 8.29153 4.14903]*10^-3;
a4=.001658;
b4=-1;
aa4=1.656*10^-9;
Area4=L4*W4;
for i=1:10000
    w(i)=i*10^-4+24*10^-3;
    ar4(i)=w(i)*L4;
    armm4(i)=ar4(i)*10^6;
Roncalcw4(i)=a4*w(i)^b4;
Roncalca4(i)=aa4*ar4(i)^b4;
end

%% #4 Coss 
vds4_m60=[59.18E-3
86.72E-3
143.4E-3
170.8E-3
184.7E-3
212.8E-3
270.3E-3
317.9E-3
375.0E-3
449.1E-3
579.9E-3
711.2E-3
857.4E-3
1.011
1.175
1.344
1.516
1.572
1.684
1.863
1.928
1.961
2.027
2.161
2.212
2.238
2.289
2.393
2.446
2.552
2.768
2.866
3.066
3.477
4.347
5.968];
vds4_m125=[78.51E-3
155.7E-3
322.0E-3
699.0E-3
1.607
2.809
4.421
6.568];
vds4_m250=[85.72E-3
177.4E-3
377.2E-3
836.6E-3
1.855
3.135
4.853];
vds4_m500=[85.99E-3
177.7E-3
377.5E-3
837.0E-3
1.849
3.126
4.845];
vds4_m999=[86.22E-3
178.1E-3
378.1E-3
838.2E-3
1.856
3.136
4.858];
% area4_m60=L4*400*60.*ones(length(vds),1);
% area4_m125=L4*400*125.*ones(length(vds),1);
% area4_m250=L4*400*250.*ones(length(vds),1);
% area4_m500=L4*400*500.*ones(length(vds),1);
% area4_m999=L4*400*999.*ones(length(vds),1);
coss4_m60=[135.7E-12
132.5E-12
127.5E-12
123.3E-12
121.5E-12
119.3E-12
115.1E-12
110.8E-12
106.8E-12
102.4E-12
96.21E-12
90.19E-12
84.87E-12
80.15E-12
75.94E-12
72.21E-12
68.91E-12
67.22E-12
65.75E-12
63.51E-12
62.07E-12
61.50E-12
60.84E-12
59.56E-12
58.66E-12
58.29E-12
57.85E-12
56.99E-12
56.30E-12
55.49E-12
53.97E-12
52.82E-12
51.60E-12
49.35E-12
45.47E-12
40.12E-12];
coss4_m125=[278.5E-12
264.6E-12
237.5E-12
197.9E-12
151.6E-12
119.3E-12
96.84E-12
81.25E-12];
coss4_m250=[555.3E-12
521.2E-12
461.0E-12
375.6E-12
287.1E-12
226.8E-12
185.9E-12];
coss4_m500=[1.109E-9
1.042E-9
921.1E-12
751.4E-12
574.5E-12
454.3E-12
372.0E-12];
coss4_m999=[2.215E-9
2.082E-9
1.839E-9
1.501E-9
1.146E-9
906.7E-12
742.2E-12];
c4_m60=1.112E-10;
d4_m60=-.1902;
e4_m60=-3.486E-11;
for i=1:5000
    vds=i*.001;
coss4_m60eq(i)=c4_m60*vds^d4_m60+e4_m60; 
end
c4_m250=4.557E-10;
d4_m250=-.1865;
e4_m250=-1.262E-10;
for i=1:5000
    vds=i*.001;
coss4_m250eq(i)=c4_m250*vds^d4_m250+e4_m250; 
end
c4_m125=2.271E-10;
d4_m125=-.189;
e4_m125=-6.4E-11;
for i=1:5000
    vds=i*.001;
coss4_m125eq(i)=c4_m125*vds^d4_m125+e4_m125; 
end
c4_m500=9.172E-10;
d4_m500=-.186;
e4_m500=-2.542E-10;
for i=1:5000
    vds=i*.001;
coss4_m500eq(i)=c4_m500*vds^d4_m500+e4_m500; 
end
c4_m999=1.833E-09;
d4_m999=-.1861;
e4_m999=-5.085E-10;
for i=1:5000
    vds=i*.001;
coss4_m999eq(i)=c4_m999*vds^d4_m999+e4_m999; 
end
coss4_vds_4=[coss4_m60eq(4000) coss4_m125eq(4000) coss4_m250eq(4000) coss4_m500eq(4000) coss4_m999eq(4000)];

%% #4 Qg
vgs4=[250.0E-3
500.0E-3
750.0E-3
1.000
1.250
1.500
1.750
2.000
2.250
2.500
2.750
3.000
3.250
3.500
3.750
4.000
4.250
4.500
4.750
5.000];
Qg4_m60=[7.570E-12
14.86E-12
22.66E-12
184.5E-12
210.1E-12
232.2E-12
253.6E-12
274.9E-12
296.2E-12
317.6E-12
338.9E-12
360.2E-12
381.6E-12
402.9E-12
424.3E-12
445.7E-12
467.1E-12
488.5E-12
509.9E-12
531.3E-12];
Qg4_m125=[15.79E-12
31.00E-12
47.35E-12
387.3E-12
439.5E-12
484.6E-12
528.9E-12
573.3E-12
617.6E-12
662.0E-12
706.4E-12
750.9E-12
795.4E-12
839.9E-12
884.4E-12
929.0E-12
973.6E-12
1.018E-9
1.063E-9
1.107E-9];
Qg4_m250=[31.60E-12
62.05E-12
95.47E-12
780.3E-12
880.9E-12
970.0E-12
1.058E-9
1.147E-9
1.236E-9
1.324E-9
1.413E-9
1.502E-9
1.591E-9
1.680E-9
1.769E-9
1.858E-9
1.948E-9
2.037E-9
2.126E-9
2.215E-9];
Qg4_m500=[63.22E-12
124.1E-12
219.3E-12
1.572E-9
1.764E-9
1.941E-9
2.117E-9
2.294E-9
2.472E-9
2.649E-9
2.827E-9
3.005E-9
3.183E-9
3.361E-9
3.539E-9
3.717E-9
3.896E-9
4.074E-9
4.253E-9
4.431E-9];
Qg4_m999=[126.3E-12
248.1E-12
482.9E-12
3.159E-9
3.527E-9
3.878E-9
4.231E-9
4.585E-9
4.939E-9
5.294E-9
5.649E-9
6.004E-9
6.360E-9
6.715E-9
7.071E-9
7.428E-9
7.784E-9
8.140E-9
8.497E-9
8.854E-9];
Qg4_vgs_5=[Qg4_m60(20,1) Qg4_m125(20,1) Qg4_m250(20,1) Qg4_m500(20,1) Qg4_m999(20,1)];
%% 12HVnLDMOS nbl hp mac #5
%fingers =2
L5=0.9*10^-6;
W5=[60 125 250 500 999].*(400*10^-6);
Ron5=[60.4806 29.0297 14.5124 7.25405 3.62885]*10^-3;
a5=.001451;
b5=-1;
aa5=1.302*10^-9;
Area5=L5*W5;
for i=1:10000
    w(i)=i*10^-4+24*10^-3;
    ar5(i)=w(i)*L5;
    armm5(i)=ar5(i)*10^6;
Roncalcw5(i)=a5*w(i)^b5;
Roncalca5(i)=aa5*ar5(i)^b5;
end
%% #5 Coss 
vds5_m60=[56.66E-3
105.8E-3
209.1E-3
435.4E-3
956.2E-3
1.837
3.001
4.584
6.675];
vds5_m125=[54.29E-3
105.7E-3
214.1E-3
452.4E-3
1.003
1.904
3.093
4.706
6.835];
vds5_m250=[51.83E-3
103.2E-3
211.4E-3
449.4E-3
999.5E-3
1.904
3.090
4.704
6.832];
vds5_m500=[67.83E-3
137.5E-3
286.7E-3
621.4E-3
1.412
2.425
3.795
5.640];
vds5_m999=[67.11E-3
136.8E-3
286.1E-3
621.1E-3
1.412
2.425
3.796
5.641];
% area5_m60=L5*400*60.*ones(length(vds),1);
% area5_m125=L5*400*125.*ones(length(vds),1);
% area5_m250=L5*400*250.*ones(length(vds),1);
% area5_m500=L5*400*500.*ones(length(vds),1);
% area5_m999=L5*400*999.*ones(length(vds),1);
coss5_m60=[132.6E-12
128.3E-12
119.5E-12
105.2E-12
87.42E-12
72.32E-12
61.07E-12
53.16E-12
46.75E-12];
coss5_m125=[277.1E-12
266.9E-12
248.3E-12
217.0E-12
180.3E-12
148.7E-12
126.2E-12
109.6E-12
96.81E-12];
coss5_m250=[554.3E-12
535.5E-12
496.8E-12
435.3E-12
360.4E-12
298.0E-12
252.0E-12
219.6E-12
193.4E-12];
coss5_m500=[1.098E-9
1.044E-9
949.1E-12
809.2E-12
659.0E-12
547.5E-12
471.1E-12
412.9E-12];
coss5_m999=[2.191E-9
2.090E-9
1.894E-9
1.619E-9
1.315E-9
1.096E-9
939.8E-12
826.2E-12];
c5_m60=1.173E-10;
d5_m60=-.1487;
e5_m60=-3.608E-11;
for i=1:5000
    vds=i*.001;
coss5_m60eq(i)=c5_m60*vds^d5_m60+e5_m60; 
end
c5_m250=4.889E-10;
d5_m250=-.1464;
e5_m250=-1.513E-10;
for i=1:5000
    vds=i*.001;
coss5_m250eq(i)=c5_m250*vds^d5_m250+e5_m250; 
end
c5_m125=2.418E-10;
d5_m125=-.1492;
e5_m125=-7.288E-11;
for i=1:5000
    vds=i*.001;
coss5_m125eq(i)=c5_m125*vds^d5_m125+e5_m125; 
end
c5_m500=9.374E-10;
d5_m500=-.1578;
e5_m500=-2.566E-10;
for i=1:5000
    vds=i*.001;
coss5_m500eq(i)=c5_m500*vds^d5_m500+e5_m500; 
end
c5_m999=1.888E-09;
d5_m999=-.1562;
e5_m999=-5.27E-10;
for i=1:5000
    vds=i*.001;
coss5_m999eq(i)=c5_m999*vds^d5_m999+e5_m999; 
end
coss5_vds_4=[coss5_m60eq(4000) coss5_m125eq(4000) coss5_m250eq(4000) coss5_m500eq(4000) coss5_m999eq(4000)];

%% #5 Qg
vgs5=[250.0E-3
500.0E-3
750.0E-3
1.000
1.250
1.500
1.750
2.000
2.250
2.500
2.750
3.000
3.250
3.500
3.750
4.000
4.250
4.500
4.750
5.000];
Qg5_m60=[9.952E-12
18.54E-12
26.11E-12
32.83E-12
229.9E-12
256.2E-12
279.6E-12
302.6E-12
325.6E-12
348.6E-12
371.5E-12
394.5E-12
417.5E-12
440.5E-12
463.5E-12
486.5E-12
509.6E-12
532.6E-12
555.7E-12
578.7E-12];
Qg5_m125=[20.73E-12
38.63E-12
54.41E-12
109.5E-12
482.6E-12
535.1E-12
583.2E-12
630.9E-12
678.7E-12
726.5E-12
774.2E-12
822.1E-12
870.0E-12
917.9E-12
965.8E-12
1.014E-9
1.062E-9
1.110E-9
1.158E-9
1.206E-9];
Qg5_m250=[41.47E-12
77.26E-12
108.8E-12
272.7E-12
971.0E-12
1.072E-9
1.167E-9
1.262E-9
1.358E-9
1.453E-9
1.549E-9
1.644E-9
1.740E-9
1.836E-9
1.932E-9
2.028E-9
2.124E-9
2.220E-9
2.316E-9
2.412E-9];
Qg5_m500=[82.93E-12
154.5E-12
217.7E-12
592.6E-12
1.950E-9
2.145E-9
2.335E-9
2.525E-9
2.716E-9
2.907E-9
3.098E-9
3.289E-9
3.481E-9
3.672E-9
3.864E-9
4.056E-9
4.247E-9
4.439E-9
4.631E-9
4.823E-9];
Qg5_m999=[165.7E-12
308.7E-12
435.1E-12
1.206E-9
3.907E-9
4.287E-9
4.665E-9
5.045E-9
5.426E-9
5.808E-9
6.190E-9
6.572E-9
6.954E-9
7.337E-9
7.720E-9
8.103E-9
8.487E-9
8.870E-9
9.254E-9
9.637E-9];
Qg5_vgs_5=[Qg5_m60(20,1) Qg5_m125(20,1) Qg5_m250(20,1) Qg5_m500(20,1) Qg5_m999(20,1)];
%% 12HVnLDMOS nbl mac #6
%fingers =2
L6=1*10^-6;
W6=[60 125 250 500 999].*(400*10^-6);
Ron6=[66.3225 31.8317 15.913 7.95387 3.9787]*10^-3;
a6=.00159;
b6=-1;
aa6=1.585*10^-9;
Area6=L6*W6;
for i=1:10000
    w(i)=i*10^-4+24*10^-3;
    ar6(i)=w(i)*L6;
    armm6(i)=ar6(i)*10^6;
Roncalcw6(i)=a6*w(i)^b6;
Roncalca6(i)=aa6*ar6(i)^b6;
end
%% #6 Coss 
vds6_m60=[84.77E-3
159.7E-3
272.9E-3
524.1E-3
1.113
2.078
3.308
4.970];
vds6_m125=[89.05E-3
162.4E-3
321.0E-3
683.5E-3
1.562
2.619
4.031
5.944];
vds6_m250=[97.37E-3
198.1E-3
420.9E-3
945.3E-3
1.891
3.061
4.639
6.753];
vds6_m500=[50.17E-3
102.9E-3
214.8E-3
464.8E-3
1.058
2.062
3.286
4.943];
vds6_m999=[50.46E-3
103.3E-3
215.4E-3
465.7E-3
1.060
2.059
3.285
4.940];
% area6_m60=L6*400*60.*ones(length(vds),1);
% area6_m125=L6*400*125.*ones(length(vds),1);
% area6_m250=L6*400*250.*ones(length(vds),1);
% area6_m500=L6*400*500.*ones(length(vds),1);
% area6_m999=L6*400*999.*ones(length(vds),1);
coss6_m60=[148.0E-12
141.1E-12
132.4E-12
119.4E-12
101.8E-12
85.25E-12
73.28E-12
64.41E-12];
coss6_m125=[305.8E-12
293.1E-12
271.0E-12
237.1E-12
195.8E-12
164.4E-12
143.4E-12
126.8E-12];
coss6_m250=[608.8E-12
576.4E-12
520.7E-12
442.7E-12
367.2E-12
312.9E-12
274.4E-12
243.5E-12];
coss6_m500=[1.252E-9
1.214E-9
1.144E-9
1.024E-9
863.2E-12
715.3E-12
611.9E-12
537.9E-12];
coss6_m999=[2.500E-9
2.424E-9
2.284E-9
2.045E-9
1.724E-9
1.429E-9
1.223E-9
1.075E-9];
c6_m60=1.377E-10;
d6_m60=-.1432;
e6_m60=-3.893E-11;
for i=1:5000
    vds=i*.001;
coss6_m60eq(i)=c6_m60*vds^d6_m60+e6_m60; 
end
c6_m250=5.403E-10;
d6_m250=-.1585;
e6_m250=-1.327E-10;
for i=1:5000
    vds=i*.001;
coss6_m250eq(i)=c6_m250*vds^d6_m250+e6_m250; 
end
c6_m125=2.812E-10;
d6_m125=-.1481;
e6_m125=-7.683E-11;
for i=1:5000
    vds=i*.001;
coss6_m125eq(i)=c6_m125*vds^d6_m125+e6_m125; 
end
c6_m500=1.258E-09;
d6_m500=-.1176;
e6_m500=-4.374E-10;
for i=1:5000
    vds=i*.001;
coss6_m500eq(i)=c6_m500*vds^d6_m500+e6_m500; 
end
c6_m999=2.512E-09;
d6_m999=-.1177;
e6_m999=-8.719E-10;
for i=1:5000
    vds=i*.001;
coss6_m999eq(i)=c6_m999*vds^d6_m999+e6_m999; 
end
coss6_vds_4=[coss6_m60eq(4000) coss6_m125eq(4000) coss6_m250eq(4000) coss6_m500eq(4000) coss6_m999eq(4000)];

%% #6 Qg
vgs6=[250.0E-3
500.0E-3
750.0E-3
1.000
1.250
1.500
1.750
2.000
2.250
2.500
2.750
3.000
3.250
3.500
3.750
4.000
4.250
4.500
4.750
5.000];
Qg6_m60=[6.941E-12
12.90E-12
173.2E-12
199.8E-12
222.1E-12
243.5E-12
264.8E-12
286.1E-12
307.4E-12
328.7E-12
350.0E-12
371.3E-12
392.7E-12
414.0E-12
435.4E-12
456.7E-12
478.1E-12
499.5E-12
520.9E-12
542.3E-12];
Qg6_m125=[14.48E-12
27.28E-12
364.1E-12
418.2E-12
463.6E-12
508.0E-12
552.3E-12
596.5E-12
640.9E-12
685.2E-12
729.6E-12
774.1E-12
818.5E-12
863.0E-12
907.5E-12
952.0E-12
996.6E-12
1.041E-9
1.086E-9
1.130E-9];
Qg6_m250=[28.98E-12
103.4E-12
734.4E-12
838.8E-12
928.1E-12
1.017E-9
1.105E-9
1.194E-9
1.282E-9
1.371E-9
1.460E-9
1.549E-9
1.638E-9
1.726E-9
1.815E-9
1.905E-9
1.994E-9
2.083E-9
2.172E-9
2.261E-9];
Qg6_m500=[57.99E-12
299.8E-12
1.481E-9
1.680E-9
1.857E-9
2.034E-9
2.210E-9
2.388E-9
2.565E-9
2.742E-9
2.920E-9
3.098E-9
3.275E-9
3.453E-9
3.632E-9
3.810E-9
3.988E-9
4.166E-9
4.345E-9
4.523E-9];
Qg6_m999=[115.9E-12
682.2E-12
2.982E-9
3.360E-9
3.712E-9
4.064E-9
4.417E-9
4.771E-9
5.124E-9
5.479E-9
5.834E-9
6.189E-9
6.544E-9
6.900E-9
7.256E-9
7.611E-9
7.968E-9
8.325E-9
8.681E-9
9.037E-9];
Qg6_vgs_5=[Qg6_m60(20,1) Qg6_m125(20,1) Qg6_m250(20,1) Qg6_m500(20,1) Qg6_m999(20,1)];
%% 12HVnLDMOS nbl mr mac #7
%fingers =2
L7=.9*10^-6;
W7=[60 125 250 500 999].*(400*10^-6);
Ron7=[106.539 51.1386 25.568 12.787 6.39517]*10^-3;
a7=.002557;
b7=-1;
aa7=2.3*10^-9;
Area7=L7*W7;
for i=1:10000
    w(i)=i*10^-4+24*10^-3;
    ar7(i)=w(i)*L7;
    armm7(i)=ar7(i)*10^6;
Roncalcw7(i)=a7*w(i)^b7;
Roncalca7(i)=aa7*ar7(i)^b7;
end
%% #7 Coss 
vds7_m60=[53.45E-3
95.84E-3
163.6E-3
296.7E-3
611.7E-3
1.401
2.316
3.570
5.347];
vds7_m125=[73.99E-3
130.5E-3
255.3E-3
550.9E-3
1.295
2.196
3.406
5.122];
vds7_m250=[62.66E-3
118.6E-3
242.1E-3
535.2E-3
1.275
2.182
3.388
5.098];
vds7_m500=[57.05E-3
112.7E-3
235.5E-3
527.4E-3
1.265
2.174
3.378
5.084];
vds7_m999=[54.31E-3
109.8E-3
232.6E-3
524.2E-3
1.261
2.171
3.373
5.078];
% area7_m60=L7*400*60.*ones(length(vds),1);
% area7_m125=L7*400*125.*ones(length(vds),1);
% area7_m250=L7*400*250.*ones(length(vds),1);
% area7_m500=L7*400*500.*ones(length(vds),1);
% area7_m999=L7*400*999.*ones(length(vds),1);
coss7_m60=[162.1E-12
159.9E-12
149.2E-12
134.3E-12
113.5E-12
90.59E-12
74.69E-12
64.50E-12
56.37E-12];
coss7_m125=[336.2E-12
318.5E-12
288.5E-12
243.6E-12
193.6E-12
158.8E-12
136.6E-12
119.2E-12];
coss7_m250=[680.1E-12
644.1E-12
582.8E-12
491.3E-12
389.5E-12
318.5E-12
273.7E-12
238.8E-12];
coss7_m500=[1.367E-9
1.295E-9
1.172E-9
986.7E-12
781.2E-12
638.1E-12
547.9E-12
478.1E-12];
coss7_m999=[2.739E-9
2.593E-9
2.347E-9
1.975E-9
1.562E-9
1.276E-9
1.095E-9
955.6E-12];
c7_m60=1.501E-10;
d7_m60=-.1482;
e7_m60=-5.477E-11;
for i=1:5000
    vds=i*.001;
coss7_m60eq(i)=c7_m60*vds^d7_m60+e7_m60; 
end
c7_m250=5.36E-10;
d7_m250=-.1737;
e7_m250=-1.409E-10;
for i=1:5000
    vds=i*.001;
coss7_m250eq(i)=c7_m250*vds^d7_m250+e7_m250; 
end
c7_m125=2.629E-10;
d7_m125=-.1819;
e7_m125=-6.527E-11;
for i=1:5000
    vds=i*.001;
coss7_m125eq(i)=c7_m125*vds^d7_m125+e7_m125; 
end
c7_m500=1.085E-09;
d7_m500=-.1689;
e7_m500=-2.95E-10;
for i=1:5000
    vds=i*.001;
coss7_m500eq(i)=c7_m500*vds^d7_m500+e7_m500; 
end
c7_m999=2.176E-09;
d7_m999=-.1669;
e7_m999=-5.98E-10;
for i=1:5000
    vds=i*.001;
coss7_m999eq(i)=c7_m999*vds^d7_m999+e7_m999; 
end
coss7_vds_4=[coss7_m60eq(4000) coss7_m125eq(4000) coss7_m250eq(4000) coss7_m500eq(4000) coss7_m999eq(4000)];

%% #7 Qg
vgs7=[250.0E-3
500.0E-3
750.0E-3
1.000
1.250
1.500
1.750
2.000
2.250
2.500
2.750
3.000
3.250
3.500
3.750
4.000
4.250
4.500
4.750
5.000];
Qg7_m60=[6.991E-12
12.72E-12
175.9E-12
214.3E-12
244.8E-12
274.0E-12
303.0E-12
332.0E-12
361.0E-12
390.0E-12
419.0E-12
448.1E-12
477.2E-12
506.2E-12
535.3E-12
564.4E-12
593.5E-12
622.6E-12
651.7E-12
680.8E-12];
Qg7_m125=[14.56E-12
26.54E-12
371.2E-12
449.3E-12
511.2E-12
571.6E-12
631.9E-12
692.2E-12
752.6E-12
813.0E-12
873.5E-12
934.0E-12
994.5E-12
1.055E-9
1.116E-9
1.176E-9
1.237E-9
1.297E-9
1.358E-9
1.419E-9];
Qg7_m250=[29.13E-12
53.22E-12
751.4E-12
901.9E-12
1.024E-9
1.144E-9
1.264E-9
1.385E-9
1.506E-9
1.627E-9
1.747E-9
1.868E-9
1.989E-9
2.110E-9
2.231E-9
2.353E-9
2.474E-9
2.595E-9
2.716E-9
2.838E-9];
Qg7_m500=[58.26E-12
107.3E-12
1.522E-9
1.808E-9
2.049E-9
2.289E-9
2.530E-9
2.771E-9
3.012E-9
3.254E-9
3.495E-9
3.737E-9
3.979E-9
4.221E-9
4.463E-9
4.706E-9
4.948E-9
5.191E-9
5.433E-9
5.676E-9];
Qg7_m999=[116.4E-12
250.8E-12
3.077E-9
3.616E-9
4.095E-9
4.574E-9
5.055E-9
5.536E-9
6.018E-9
6.501E-9
6.984E-9
7.467E-9
7.951E-9
8.434E-9
8.918E-9
9.403E-9
9.887E-9
10.37E-9
10.86E-9
11.34E-9];
Qg7_vgs_5=[Qg7_m60(20,1) Qg7_m125(20,1) Qg7_m250(20,1) Qg7_m500(20,1) Qg7_m999(20,1)];
%% Plotting
plot(armm1,Roncalca1)
hold on
plot(armm2,Roncalca2)
hold on
plot(armm3,Roncalca3)
hold on
plot(armm4,Roncalca4)
hold on
plot(armm5,Roncalca5)
hold on
plot(armm6,Roncalca6)
hold on
plot(armm7,Roncalca7)
title('Ron vs. Area')
xlabel('Area (mm^2)')
ylabel('Ron (ohm)')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
ylim([0 .01])
% figure
% plot(Area1,Ron1)
% hold on
% plot(Area2,Ron2)
% hold on
% plot(Area3,Ron3)
% hold on
% plot(Area4,Ron4)
% hold on
% plot(Area5,Ron5)
% hold on
% plot(Area6,Ron6)
% hold on
% plot(Area7,Ron7)
% title('Ron vs. Area')
% xlabel('Area (m)')
% ylabel('Ron (ohm)')
% legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
% ymin([0 .1])
%% Plot Coss
figure
plot(vds1_m999,coss1_m999)
hold on
plot(vds2_m999,coss2_m999)
hold on
plot(vds3_m999,coss3_m999)
hold on
plot(vds4_m999,coss4_m999)
hold on
plot(vds5_m999,coss5_m999)
hold on
plot(vds6_m999,coss6_m999)
hold on
plot(vds7_m999,coss7_m999)
hold on
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds at W=0.4m')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(vds1_m500,coss1_m500)
hold on
plot(vds2_m500,coss2_m500)
hold on
plot(vds3_m500,coss3_m500)
hold on
plot(vds4_m500,coss4_m500)
hold on
plot(vds5_m500,coss5_m500)
hold on
plot(vds6_m500,coss6_m500)
hold on
plot(vds7_m500,coss7_m500)
hold on
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds at W=0.2m')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(vds1_m250,coss1_m250)
hold on
plot(vds2_m250,coss2_m250)
hold on
plot(vds3_m250,coss3_m250)
hold on
plot(vds4_m250,coss4_m250)
hold on
plot(vds5_m250,coss5_m250)
hold on
plot(vds6_m250,coss6_m250)
hold on
plot(vds7_m250,coss7_m250)
hold on
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds at W=0.1m')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(vds1_m125,coss1_m125)
hold on
plot(vds2_m125,coss2_m125)
hold on
plot(vds3_m125,coss3_m125)
hold on
plot(vds4_m125,coss4_m125)
hold on
plot(vds5_m125,coss5_m125)
hold on
plot(vds6_m125,coss6_m125)
hold on
plot(vds7_m125,coss7_m125)
hold on
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds at W=0.05m')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(vds1_m60,coss1_m60)
hold on
plot(vds2_m60,coss2_m60)
hold on
plot(vds3_m60,coss3_m60)
hold on
plot(vds4_m60,coss4_m60)
hold on
plot(vds5_m60,coss5_m60)
hold on
plot(vds6_m60,coss6_m60)
hold on
plot(vds7_m60,coss7_m60)
hold on
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds at W=0.024m')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(vds1_m60,coss1_m60)
hold on
plot(vds1_m125,coss1_m125)
hold on
plot(vds1_m250,coss1_m250)
hold on
plot(vds1_m500,coss1_m500)
hold on
plot(vds1_m999,coss1_m999)
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds for 8Vnbl')
legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
figure
plot(vds2_m60,coss2_m60)
hold on
plot(vds2_m125,coss2_m125)
hold on
plot(vds2_m250,coss2_m250)
hold on
plot(vds2_m500,coss2_m500)
hold on
plot(vds2_m999,coss2_m999)
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds for 8Viso')
legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
figure
plot(vds3_m60,coss3_m60)
hold on
plot(vds3_m125,coss3_m125)
hold on
plot(vds3_m250,coss3_m250)
hold on
plot(vds3_m500,coss3_m500)
hold on
plot(vds3_m999,coss3_m999)
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds for 12Visohpmac')
legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
figure
plot(vds4_m60,coss4_m60)
hold on
plot(vds4_m125,coss4_m125)
hold on
plot(vds4_m250,coss4_m250)
hold on
plot(vds4_m500,coss4_m500)
hold on
plot(vds4_m999,coss4_m999)
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds for 12Visomac')
legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
figure
plot(vds5_m60,coss5_m60)
hold on
plot(vds5_m125,coss5_m125)
hold on
plot(vds5_m250,coss5_m250)
hold on
plot(vds5_m500,coss5_m500)
hold on
plot(vds5_m999,coss5_m999)
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds for 12Vnblisomac')
legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
figure
plot(vds6_m60,coss6_m60)
hold on
plot(vds6_m125,coss6_m125)
hold on
plot(vds6_m250,coss6_m250)
hold on
plot(vds6_m500,coss6_m500)
hold on
plot(vds6_m999,coss6_m999)
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds for 12Vnblmac')
legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
figure
plot(vds7_m60,coss7_m60)
hold on
plot(vds7_m125,coss7_m125)
hold on
plot(vds7_m250,coss7_m250)
hold on
plot(vds7_m500,coss7_m500)
hold on
plot(vds7_m999,coss7_m999)
xlabel('Vds (V)')
ylabel('Coss (F)')
title('Coss vs. Vds for 12Vnblmrmac')
legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
figure
plot(W1,coss1_vds_4)
hold on
plot(W2,coss2_vds_4)
hold on
plot(W3,coss3_vds_4)
hold on
plot(W4,coss4_vds_4)
hold on
plot(W5,coss5_vds_4)
hold on
plot(W6,coss6_vds_4)
hold on
plot(W7,coss7_vds_4)
xlabel('Width (m)')
ylabel('Coss (F)')
title('Coss vs. Width at Vds=4V')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
%% Plot Qg
% plot(Qg1_m60,vgs1)
% hold on
% plot(Qg1_m125,vgs1)
% hold on
% plot(Qg1_m250,vgs1)
% hold on
% plot(Qg1_m500,vgs1)
% hold on
% plot(Qg1_m999,vgs1)
% xlabel('Qg (C)')
% ylabel('Vgs (V)')
% title('Vgs vs. Qg for 8Vnbl')
% legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
% figure
% plot(Qg2_m60,vgs2)
% hold on
% plot(Qg2_m125,vgs2)
% hold on
% plot(Qg2_m250,vgs2)
% hold on
% plot(Qg2_m500,vgs2)
% hold on
% plot(Qg2_m999,vgs2)
% xlabel('Qg (C)')
% ylabel('Vgs (V)')
% title('Vgs vs. Qg for 8Viso')
% legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
% figure
% plot(Qg3_m60,vgs3)
% hold on
% plot(Qg3_m125,vgs3)
% hold on
% plot(Qg3_m250,vgs3)
% hold on
% plot(Qg3_m500,vgs3)
% hold on
% plot(Qg3_m999,vgs3)
% xlabel('Qg (C)')
% ylabel('Vgs (V)')
% title('Vgs vs. Qg for 12Visohpmac')
% legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
% figure
% plot(Qg4_m60,vgs4)
% hold on
% plot(Qg4_m125,vgs4)
% hold on
% plot(Qg4_m250,vgs4)
% hold on
% plot(Qg4_m500,vgs4)
% hold on
% plot(Qg4_m999,vgs4)
% xlabel('Qg (C)')
% ylabel('Vgs (V)')
% title('Vgs vs. Qg for 12Visomac')
% legend('w=24mm','w=50mm','w=100mm','w=200mm','w=400mm')
% figure
% plot(W5,Qg5)
% xlabel('Width (m)')
% ylabel('Qg (C)')
% title('Qg vs. Width for 12Vnblisomac')
% figure
% plot(W6,Qg6)
% xlabel('Width (m)')
% ylabel('Qg (C)')
% title('Qg vs. Width for 12Vnblmac')
% figure
% plot(W7,Qg7)
% xlabel('Width (m)')
% ylabel('Qg (C)')
% title('Qg vs. Width for 12Vnblmrmac')
figure
plot(Qg1_m60,vgs1)
hold on
plot(Qg2_m60,vgs2)
hold on
plot(Qg3_m60,vgs3)
hold on
plot(Qg4_m60,vgs4)
hold on
plot(Qg5_m60,vgs5)
hold on
plot(Qg6_m60,vgs6)
hold on
plot(Qg7_m60,vgs7)
hold on
xlabel('Qg (C)')
ylabel('Vgs (V)')
title('Vgs vs. Qg for W=24mm')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(Qg1_m125,vgs1)
hold on
plot(Qg2_m125,vgs2)
hold on
plot(Qg3_m125,vgs3)
hold on
plot(Qg4_m125,vgs4)
hold on
plot(Qg5_m125,vgs5)
hold on
plot(Qg6_m125,vgs6)
hold on
plot(Qg7_m125,vgs7)
hold on
xlabel('Qg (C)')
ylabel('Vgs (V)')
title('Vgs vs. Qg for W=50mm')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(Qg1_m250,vgs1)
hold on
plot(Qg2_m250,vgs2)
hold on
plot(Qg3_m250,vgs3)
hold on
plot(Qg4_m250,vgs4)
hold on
plot(Qg5_m250,vgs5)
hold on
plot(Qg6_m250,vgs6)
hold on
plot(Qg7_m250,vgs7)
hold on
xlabel('Qg (C)')
ylabel('Vgs (V)')
title('Vgs vs. Qg for W=100mm')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(Qg1_m500,vgs1)
hold on
plot(Qg2_m500,vgs2)
hold on
plot(Qg3_m500,vgs3)
hold on
plot(Qg4_m500,vgs4)
hold on
plot(Qg5_m500,vgs5)
hold on
plot(Qg6_m500,vgs6)
hold on
plot(Qg7_m500,vgs7)
hold on
xlabel('Qg (C)')
ylabel('Vgs (V)')
title('Vgs vs. Qg for W=200mm')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(Qg1_m999,vgs1)
hold on
plot(Qg2_m999,vgs2)
hold on
plot(Qg3_m999,vgs3)
hold on
plot(Qg4_m999,vgs4)
hold on
plot(Qg5_m999,vgs5)
hold on
plot(Qg6_m999,vgs6)
hold on
plot(Qg7_m999,vgs7)
hold on
xlabel('Qg (C)')
ylabel('Vgs (V)')
title('Vgs vs. Qg for W=400mm')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
figure
plot(Qg1_vgs_5,W1)
hold on 
plot(Qg2_vgs_5,W2)
hold on 
plot(Qg3_vgs_5,W2)
hold on 
plot(Qg4_vgs_5,W3)
hold on 
plot(Qg5_vgs_5,W5)
hold on 
plot(Qg6_vgs_5,W6)
hold on 
plot(Qg7_vgs_5,W7)
hold on 
xlabel('Qg (C)')
ylabel('Width (m)')
title('Width vs. Qg for Vgs=5V')
legend('8HVnLDMOS nbl','8HVnLDMOS iso','12HVnLDMOS iso hp mac','12HVnLDMOS iso mac','12HVnLDMOS nbl hp mac','12HVnLDMOS nbl mac','12HVnLDMOS nbl mr mac')
