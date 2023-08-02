%PWM_COMPEL_2023 Gives the gate values for the SC Fib Converter
%controller

% The indexes for the rise time do not need to subtract 1 because the
% rise is actually starts during the d ON time so it cancels with the
% -1 for index starting at zero

% The fall time indexes need to subtract 1 because they occur duing
% the d ON period and must account for indexing starting at zero

%% Variable Declaration
fpga_fs = 133e6;
fpga_multi = 1;
recomended_fs = 1e6;
%max_mod_value = 2*floor(fpga_fs*fpga_multi/recomended_fs/2);
max_mod_value = 132;
Ts = 1/(fpga_fs*fpga_multi)*max_mod_value;

d = 0.8;
dt1 = 0.001;
dt2 = 0.001;

dt1_length = ceil(dt1*max_mod_value);
dt2_length = ceil(dt2*max_mod_value);
dt1_length = 0;
dt2_length = 0;

mod_for_d = max_mod_value - (dt1_length*2 + dt1_length*2);
d_length = floor(d*mod_for_d/2);
d_prime_length = (mod_for_d-2*d_length)/2;

%% Added for the use of approx 24 ns delay of half bridge gate drivers

gate_driver_delay = 16e-9;


clk_ts = 1/(fpga_fs*fpga_multi);

offset = round(gate_driver_delay/clk_ts);

difference_between = abs(offset*clk_ts-gate_driver_delay);





        g1.rise = 0;
        g1.fall = (-1) + d_length + dt1_length + d_prime_length;

        g2.rise = d_length + dt1_length + d_prime_length + dt2_length;
        g2.fall = (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length + d_prime_length;

        g3.rise = offset + 0;
        g3.fall = offset + (-1) + d_length;

        g4.rise = 999;
        g4.fall = 0;

        g5.rise = d_length + dt1_length + d_prime_length + dt2_length;
        g5.fall = (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length + d_prime_length;

        g6.rise = offset + 0;
        g6.fall = offset + (-1) + d_length + dt1_length + d_prime_length;

        g7.rise = offset + d_length + dt1_length + d_prime_length + dt2_length;
        g7.fall = offset + (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length;

        g8.rise = offset + 0;
        g8.fall = offset + (-1) + d_length;

        g9.rise = offset + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length ;
        g9.fall = offset + (-1) + d_length + dt1_length + d_prime_length;

        g10.rise = offset + d_length + dt1_length;
        g10.fall = offset + (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length + d_prime_length;

        g11.rise = d_length + dt1_length + d_prime_length + dt2_length;
        g11.fall = (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length + d_prime_length;

        g12.rise = 0;
        g12.fall = (-1) + d_length + dt1_length + d_prime_length;

        g13.rise = offset + d_length + dt1_length;
        g13.fall = offset + (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length + d_prime_length;

        g14.rise = 0;
        g14.fall = (-1) + d_length + dt1_length + d_prime_length;

        g15.rise = d_length + dt1_length + d_prime_length + dt2_length;
        g15.fall = (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length + d_prime_length;

        g16.rise = offset + d_length + dt1_length + d_prime_length + dt2_length;
        g16.fall = offset + (-1) + d_length + dt1_length + d_prime_length + dt2_length + d_length + dt1_length + d_prime_length;



% if ideal_input_Voltage_index ==4
%     g6.fall = g6.fall-max_mod_value;
% end
% 
% if ideal_input_Voltage_index ==5 % && dt2_length ==0
%     g16.fall = g16.fall-max_mod_value;
% end
%%{
%if ideal_input_Voltage_index ==5 && dt2_length >0
    g16.fall = g16.fall-max_mod_value;
    g13.fall = g13.fall-max_mod_value;
    g10.fall = g10.fall-max_mod_value;
%end
%}

fprintf('----------\n')
fprintf('Gate 1 fall: %d\n',g1.fall)
fprintf('Gate 1 rise: %d\n',g1.rise)
fprintf('----------\n')
fprintf('Gate 2 fall: %d\n',g2.fall)
fprintf('Gate 2 rise: %d\n',g2.rise)
fprintf('----------\n')
fprintf('Gate 3 fall: %d\n',g3.fall)
fprintf('Gate 3 rise: %d\n',g3.rise)
fprintf('----------\n')
fprintf('Gate 4 fall: %d\n',g4.fall)
fprintf('Gate 4 rise: %d\n',g4.rise)
fprintf('----------\n')
fprintf('Gate 5 fall: %d\n',g5.fall)
fprintf('Gate 5 rise: %d\n',g5.rise)
fprintf('----------\n')
fprintf('Gate 6 fall: %d\n',g6.fall)
fprintf('Gate 6 rise: %d\n',g6.rise)
fprintf('----------\n')
fprintf('Gate 7 fall: %d\n',g7.fall)
fprintf('Gate 7 rise: %d\n',g7.rise)
fprintf('----------\n')
fprintf('Gate 8 fall: %d\n',g8.fall)
fprintf('Gate 8 rise: %d\n',g8.rise)
fprintf('----------\n')
fprintf('Gate 9 fall: %d\n',g9.fall)
fprintf('Gate 9 rise: %d\n',g9.rise)
fprintf('----------\n')
fprintf('Gate 10 fall: %d\n',g10.fall)
fprintf('Gate 10 rise: %d\n',g10.rise)
fprintf('----------\n')
fprintf('Gate 11 fall: %d\n',g11.fall)
fprintf('Gate 11 rise: %d\n',g11.rise)
fprintf('----------\n')
fprintf('Gate 12 fall: %d\n',g12.fall)
fprintf('Gate 12 rise: %d\n',g12.rise)
fprintf('----------\n')
fprintf('Gate 13 fall: %d\n',g13.fall)
fprintf('Gate 13 rise: %d\n',g13.rise)
fprintf('----------\n')
fprintf('Gate 14 fall: %d\n',g14.fall)
fprintf('Gate 14 rise: %d\n',g14.rise)
fprintf('----------\n')
fprintf('Gate 15 fall: %d\n',g15.fall)
fprintf('Gate 15 rise: %d\n',g15.rise)
fprintf('----------\n')
fprintf('Gate 16 fall: %d\n',g16.fall)
fprintf('Gate 16 rise: %d\n',g16.rise)
fprintf('----------\n')



Numerical_Stack = [
  g1.rise g1.fall
  g2.rise g2.fall
  g3.rise g3.fall
  g4.rise g4.fall
  g5.rise g5.fall
  g6.rise g6.fall
  g7.rise g7.fall
  g8.rise g8.fall
  g9.rise g9.fall
  g10.rise g10.fall
  g11.rise g11.fall
  g12.rise g12.fall
  g13.rise g13.fall
  g14.rise g14.fall
  g15.rise g15.fall
  g16.rise g16.fall];

Numerical_Stack = Numerical_Stack +1; % MATLAB base 1 instead of zero

Waveform = zeros(16,max_mod_value*2);

for i = 1:max_mod_value*2

    if i > max_mod_value
        clock = i - max_mod_value;
    else
        clock = i;
    end



    for j = 1:16

        if Numerical_Stack(j,1) == clock
            Waveform(j,i) = 1;
        end

        if  Numerical_Stack(j,2) == clock
            Waveform(j,i) = 0;
        end

        if Numerical_Stack(j,1) ~= clock && Numerical_Stack(j,2) ~= clock && i~=1
             Waveform(j,i) = Waveform(j,i-1);
        end
    end
end

for j = 1:16
figure(129)
    ax = subplot(16,1,j);
    plot([0:max_mod_value-1],Waveform(j,max_mod_value+1:end),'Color','#7E2F8E','LineWidth',3)
    ax.YLim = [-0.1,1.1];
    ax.XLim = [0,max_mod_value-1];
    if(j<16)
        set(gca, 'Xticklabel', []);
    else
        xlabel('Count')
    end
end
