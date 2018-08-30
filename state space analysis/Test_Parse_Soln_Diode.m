% This script tests the combination of state space analysis and circuit
% parsing 

% Currently trying to get Buck converter to work 
% Then work on making it general for other tests

%% User Input

clear


Vg = 5;
L1 = 230e-9; %L
C1 = 4040e-9; %Cout
fs = 2e6;
Ts = 1/fs;
V = 1.8;
Io = 1; % was 1
M1_C = 3.4874e-10; % CHS
D1_C = 3.4874e-10; % LHS

R1 =  .01; % Rl
dt = Ts/1000;%5e-10;
Vdr = 5;
M1_R = .05; % ronHS
D1_R = .05; % ronLS

ts = [Ts*.5-dt dt Ts*.5-dt dt];


TestparseWaveform = false;



% Select .net file
filename = 'Buck2.net';
% Current options for filename:
% Boost.net
% Buck.net
% Buck2.net
% Dickson.net
% Flyback.net
% Forward.net

% Select test case file
testcase = 'TEST_ABCD_Buck_Diode';
% Current options for testcase:
% TEST_ABCD_Boost
% TEST_ABCD_Buck
% TEST_ABCD_Buck2
% TEST_ABCD_Flyback
% TEST_ABCD_Forward
% TEST_ABCD_Dickson

% Set Voltage and Current Nodes to add
Voltage = {'V1'
    'C1'
    'M1'
    'D1'};
Current = {'V1'
    'C1'
    'M1'
    'D1'};
% Change Voltage and Current based on desired output measurements (C and D
% matricies). Voltage and Current should be of type Cell 
% Example:
% Voltage = {'V1'
%     'M1'
%     'L3'};
% 
% Current = {'C1'
%     'D2'
%     'R3'};

%%% Note All Swithch, Inductor, and Capaictor elements will be included in
%%% the C and D matrix whether they are placed in the above cell array or
%%% not.


%% Run functions
parse = NetListParse();
parse.initialize(filename,Voltage,Current);
parse.ABCD();

if TestparseWaveform
testfun = str2func(testcase);
testfun(parse);
end

%% DC code

u = [Vg Io]';

top = SMPStopology();
top.Parse = parse;

conv = SMPSconverter();
conv.topology = top;
conv.ts = ts;
conv.u = u;

simulator = SMPSim();

%% To update and test
    % taken from TEST function

A = parse.Asym;
B = parse.Bsym;
C = parse.Csym;
D = parse.Dsym;

SortedTree = parse.SortedTree;
SortedCoTree = parse.SortedCoTree;

if isempty(A)
    for k = 1:1:size(parse.HtempAB,3)
        HtempAB(:,:,k) = eval(parse.HtempAB(:,:,k));
        HtempCD(:,:,k) = eval(parse.HtempCD(:,:,k));
        dependsAB(:,:,k) = eval(parse.dependsAB(:,:,k));
        savedCD(:,:,k) = eval(parse.savedCD(:,:,k));
        for j = 1:1:size(parse.DependentNames(:,k),1)
            DependentNames(j,k) = eval(parse.DependentNames{j,k});
        end
        for j = 1:1:size(parse.OutputNames(:,k),1)
            OutputNames(j,k) = eval(parse.OutputNames{j,k});
        end
    end
    for k = 1:1:size(parse.HtempAB,3)
        [A,B,C,D] = parse.loopfixAB_large(HtempAB(:,:,k),dependsAB(:,:,k),OutputNames(:,k),DependentNames(:,k));
        [C,D] = parse.loopfixCD_large(B,C,D,HtempCD(:,:,k),savedCD(:,:,k),DependentNames(:,k),SortedTree(:,:,k),SortedCoTree(:,:,k));
        parse.Anum(:,:,k)=A;
        parse.Bnum(:,:,k)=B;
        parse.Cnum(:,:,k)=C;
        parse.Dnum(:,:,k)=D;
    end
    
elseif ~isempty(parse.Anum)
    fprintf('Confirm Plecs used');
else
    for k = 1:1:size(A,3)
        parse.Anum(:,:,k) = eval(A(:,:,k));
        parse.Bnum(:,:,k) = eval(B(:,:,k));
        parse.Cnum(:,:,k) = eval(C(:,:,k));
        parse.Dnum(:,:,k) = eval(D(:,:,k));
    end
end

simulator.loadTestConverter2(conv);
Xss = simulator.SS_Soln();

%% Reconstruction of Dependent variables
% Dont need anymore due to fix of SS_Soln.m
% Dependent variables are calculated with independent variables
%{
Xss(end+1,:) = zeros(size(parse.DependentNames,1),size(Xss,2));


% For now diode is 8th row in C and D will char match dependent variables
% to find where they are in output and calculate
for i = 2:1:size(Xss,2)
    Xss(size(OutputNames,1)+1:end,i) = simulator.Cw(8,:,i-1)*Xss(:,i)+simulator.Dw(8,:,i-1)*u;
    Xss(:,1) = Xss(:,end);
end
%}

%% adjustDiodeCond


% User input to define output voltage/Cap
% identify output cap/component

Vopos = 2;
Xs = Xss;
ron = M1_R;

Order = [2 1 3 1];

parse.find_diode(Order);

check = parse.VfwdIrev(Xs,u,Order);

% Have two conditions to check physical validity of diodes
% Check to see if make sure they have no negative current (should be blocking)
% Check to see if the voltage exceeds the forward voltage of the diode.
% These do not need to know what the state of anything is wheter fet or
% power diode


simulator.Baxter_adjustDiodeCond();




Voerr = mean(Xs(Vopos,:)) - V; % Error on output
LSdiode_DT1 = Xs(1,3) < -10*ron*2;
HSdiode_DT1 = Xs(1,3) > Vg + 10*ron*2; 

LSdiode_DT2 = Xs(1,5) < -10*ron*2;  
HSdiode_DT2 = Xs(1,5) > Vg + 10*ron*2; 

hardSwNecessary_DT1 = 0;
hardSwNecessary_DT2 = 0;
hardSw_DT1 = Xs(1,3) > Vg*.01 && ts(2) < tsmax(2) && ~hardSwNecessary_DT1;
hardSw_DT2 = Xs(1,5) < Vg*.99 && ts(4) < tsmax(4) && ~hardSwNecessary_DT2;

modelError = [(abs(Voerr)>maxVoErr) 0 0;
    LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
    LSdiode_DT2 HSdiode_DT2 hardSw_DT2];
 
nattempts = 0;

while (nattempts < 100) && sum(sum(modelError))
    
	if(debug)
        [ ys, t ] = SS_WF_Reconstruct( Xs, As, Bs, ts, u );
        disp(modelError);
        disp(ts);
        plot(t,ys);
        hold on;
        ylims = ylim;
        for i = 1:length(ts)
            plot(sum(ts(1:i))*ones(1,2), ylims, ':r');
        end
%         pause
	end

    introduced_Voerr = 0;
    
    if(LSdiode_DT1 || HSdiode_DT1 || hardSw_DT1)
        [tsnew, dxsdtd, hardSwNecessary_DT1] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 3, 1, Vg, 0);
        introduced_Voerr = sum(dxsdtd(3,2:end).*(ts-tsnew));
        ts = tsnew;
    end
    
    if(LSdiode_DT2 || HSdiode_DT2 || hardSw_DT2)
        [tsnew, dxsdtd, hardSwNecessary_DT2] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 5, 1, Vg, 0);
        introduced_Voerr = introduced_Voerr + sum(dxsdtd(3,2:end).*(ts-tsnew));
        ts = tsnew;
    end
    

    %% compensate Vo for error (original) plus change from dead times
    delta_DTs = max(min(ts)/10, sum(ts)/10000);
    dXs = StateSensitivity( As, Bs, ts, u, 'ts', 1, delta_DTs, 3);
    dxsdt = (Xs - dXs)/delta_DTs;
    dt = (Voerr+introduced_Voerr)/mean(dxsdt(3,:));

    if(ts(3) - dt <0)
        dt = ts(3);
    elseif(ts(1) + dt < 0)
        dt = -ts(1);
    end
    dt = dt*.5;

    ts(1) = ts(1) + dt;
    ts(3) = ts(3) - dt;
    
    %% Recompute and reevaluate
    [ Xs] = SS_Soln( As, Bs, ts, u);
    
    Voerr = mean(Xs(3,:)) - V;
    LSdiode_DT1 = Xs(1,3) < -10*ron*2;
    HSdiode_DT1 = Xs(1,3) > Vg + 10*ron*2; 

    LSdiode_DT2 = Xs(1,5) < -10*ron*2;  
    HSdiode_DT2 = Xs(1,5) > Vg + 10*ron*2; 
    
    
    [~, ~, ~, mx1] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 3, 1, Vg, 0);
    [~, ~, ~, mx2] = adjustDiodeConduction( Xs, As, Bs, u, ts, tsmax, 5, 1, Vg, 0);
    hardSw_DT1 = (Xs(1,3) > Vg*.01 && ts(2) < tsmax(2) && ~hardSwNecessary_DT1) || mx1;
    hardSw_DT2 = (Xs(1,5) < Vg*.99 && ts(4) < tsmax(4) && ~hardSwNecessary_DT2) || mx2;
    
    modelError = [(abs(Voerr)>maxVoErr) 0 0;
        LSdiode_DT1 HSdiode_DT1 hardSw_DT1;
        LSdiode_DT2 HSdiode_DT2 hardSw_DT2];
    
    nattempts = nattempts + 1;

end



return

