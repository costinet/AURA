x = [2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.4372    0.2737    0.7965];

Coss_adj = [0 0 0 0 0 0]; 
Ron_adj = [0 0 0 0 0 0];

[stick1,graph1]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);

Coss_adj = [-0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9]; 
Ron_adj = [-4e-3 -4e-3 -4e-3 -4e-3 -4e-3 -4e-3];

[stick2,graph2]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);

Coss_adj = [0 0 0 0 0 0]; 
Ron_adj = [-4e-3 -4e-3 -4e-3 -4e-3 -4e-3 -4e-3];

[stick3,graph3]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);


Coss_adj = [-0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9]; 
Ron_adj = [0 0 0 0 0 0];

[stick9,graph9]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);




A = 1;
B = -(0.4e-9/4e-3);
C = (0.4e-9/4e-3)*16e-3-1530e-12;
x0 = 2.6e-3;
y0 = 980e-12;

B*(B*x0-A*y0)-A*C/(A^2+B^2)





Coss_adj = [0 0 0 0 0 0].*3; 
Ron_adj = [-4e-3 -4e-3 -4e-3 -4e-3 -4e-3 -4e-3].*3.35;

[stick4,graph4]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);



Coss_adj = [-0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9 -0.4e-9].*2.75; 
Ron_adj = [-4e-3 -4e-3 -4e-3 -4e-3 -4e-3 -4e-3].*2.75;


[stick5,graph5]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_CR_compare(x,Coss_adj,Ron_adj);

x = [7.0000    7.0000    7.0000    7.0000    7.0000    7.0000    7.0000    2.4372    0.2737    0.7965];