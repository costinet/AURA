function  [ron,Coss]=Select_IC_FET(FET_Number)

W_selection=linspace(0.005,0.4,100);

a =  0.001447;
b = -1;
ron = a*W_selection(FET_Number)^b;

p1 = 2.487e-9;
p2 = -6.393e-13;
Coss = p1*W_selection(FET_Number)+p2;

%L3=900*10^-9;

% tot_area = tot_area+W_selection*L3;

end