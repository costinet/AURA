x = [8.0000    8.0000    8.0000    8.0000    8.0000    8.0000    8.0000    2.4372    0.2737    0.7965];


for i = 1:1:10
    
    for j = 1:1:4
        
        Order = [6
            1
            2
            5
            11
            8
            7
            10
            3
            9
            4];
        
        [stick1,graph1]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_X_compare(x,0.00153);
        
        if j ==1
            adjust = 1;
        end
        
        
        if j ==2
            adjust = 2;
        end
        
        if j == 3
            adjust = [3 4];
        end
        if j == 4
            adjust = [5 6];
        end
        
        index=find(x(adjust(1))==Order);
        if index < 10
            x(adjust) = Order(index+1);
        else
            continue
        end
        
        [stick2,graph2]=AURA_Eff_Sweep_Disc_Buck_3L_Boost_delta_X_compare(x,0.00153);
        
        if stick2>stick1
            
            if index > 1
                x(adjust) = Order(index-1);
            end
        end
    end
    
end