



% for loop to cycle though width percentages

% Initial_guess = [0.0536  0.0341  0.1439  0.1299  0.0051  0.0417  0.0354  0.1115  1.6428e6];

Initial_guess2 = [0.1217   0.0669   0.2847   0.2795   0.0051   0.0881   0.0572   0.2082   1.5709e6];

to_change = Initial_guess2([1:4,6:8]);



lb = 0.5;
ub = 2.75;

lb2 = 0.25;
ub2 = 1.35;


L = 900*10^-9;
ratio = 2.75;

changed = ratio*to_change;

x = [changed(1:4) 0.005 changed(5:end) Initial_guess2(end)];

x(1:8)

adjust = [ones(1,length(x)-1)*20, 1e-6];

area = 2*(sum(x(1:8).*L));


rrange = linspace(lb2,ub2,30);

fval2 = zeros(1,30);
areas2 = zeros(1,30);

for i = 1:length(rrange)
    
    ratio = rrange(i);
    changed = ratio*to_change;
    
    x = [changed(1:4) 0.005 changed(5:end) Initial_guess2(end)];
    
    areas2(i) =  2*(sum(x(1:8).*L));
    
    fval2(i) = AURA_Eff_Sweep_IC(x.*adjust);
    
end

