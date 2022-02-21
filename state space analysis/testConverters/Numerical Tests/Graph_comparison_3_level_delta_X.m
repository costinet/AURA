% Script to compare changes in Numberical Components vs 

x = [4.0000    9.0000    6.0000    6.0000    4.0000    8.0000    8.0000    2.4372    0.2737    0.7965];

[stick1,graph1]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_X_compare(x,0.00153);

x = [4.0000    9.0000    6.0000    6.0000    4.0000    8.0000    8.0000    2.4372    0.2737    0.7965];

[stick2,graph2]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_X_compare(x,0.0016);

    Vgrange = 10:.25:20;
    PoutRange = 20:1:60;
    
    [f,c] = contourf(Vgrange,PoutRange,graph1-graph2,'ShowText','on');
xlabel('Vg');
ylabel('P_{out}');
% ylim([0 80])
colorbar
