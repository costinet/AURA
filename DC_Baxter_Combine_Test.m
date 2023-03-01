clear all;

niter = 0;
debug = 0;
debug2 = 0;

if ~debug
    w = warning ('off','all');
else
    w = warning ('on','all');
end

% sdir = mfilename('fullpath');
% sdir = sdir(1:find(sdir=='\',1,'last')-1);
% addpath(sdir);

%% Load test circuit
% modelfile = 'AsyncBoost'; PLECsModel = 'Boost_Async';
% modelfile = 'MRBuck'; PLECsModel = 'MRBuck';
% modelfile = 'DSC4to1'; PLECsModel = 'HDSC';
% modelfile = 'DAB'; PLECsModel = 'DAB_oneCap';
% modelfile = 'DABfull'; PLECsModel = 'DAB_8Cap';
% modelfile = 'DSC4to1Diodes'; PLECsModel = 'HDSC_withDiodes';
%modelfile = 'SC_FIB_AURA_L.net';
%modelfile = 'Buck_Boost_Vout.net';
%modelfile = 'SC_FIB_AURA.net';



%% Load for PLECS
find_system(modelfile,'SearchDepth',1, 'IncludeCommented', 'on')
open_system(modelfile,'loadonly');
circuitPath = [modelfile '/' PLECsModel];
set_param(circuitPath,'Commented','on');
simout = sim(modelfile,eps);

for i = 1:length(simout.properties)
   assignin('base',simout.properties{i},eval(['simout.' simout.properties{i}]));
end

set_param(circuitPath,'Commented','off');


%% Load for SC FIB:

%{
swvec = [
        0     1     1     0     1     0     1     0     0     1     0     1     1     0     0     1
        0     1     1     0     0     0     1     0     0     1     0     1     1     0     0     0
        0     1     1     1     0     1     1     0     0     1     0     1     1     0     1     0
        0     0     0     1     0     1     1     0     0     1     0     0     1     0     0     0
        1     0     0     1     0     1     1     0     0     1     0     0     1     1     0     1
        0     0     0     1     0     1     1     0     0     1     0     0     0     1     0     1
        0     1     1     1     0     1     1     0     0     1     1     0     0     1     0     1
        0     1     1     0     0     0     1     0     0     1     0     0     0     0     0     1
    ];

Ts = 1e-6;
Numerical_Components = [
    {'C1'   }    {[9.4000e-06]}
    {'C2'   }    {[9.4000e-06]}
    {'C3'   }    {[9.4000e-06]}
    {'C4'   }    {[9.4000e-06]}
    {'C5'   }    {[2.0000e-06]}
    {'C6'   }    {[2.0000e-06]}
    {'L1'   }    {[1.2800e-06]}
    {'L2'   }    {[1.2800e-06]}
    {'L3'   }    {[1.2800e-10]}
    {'L4'   }    {[1.2800e-10]}
    {'R1'   }    {[    0.0050]}
    {'R2'   }    {[    0.0050]}
    {'R3'   }    {[    0.0370]}
    {'R4'   }    {[    0.0370]}
    {'R5'   }    {[    0.0030]}
    {'R6'   }    {[    0.0030]}
    {'R7'   }    {[    0.0030]}
    {'R8'   }    {[    0.0030]}
    {'R9'   }    {[    0.0020]}
    {'R10'  }    {[    0.0020]}
    {'R11'  }    {[    0.0020]}
    {'R12'  }    {[    0.0020]}
    {'M1_C' }    {[4.0800e-10]}
    {'M2_C' }    {[4.0800e-10]}
    {'M3_C' }    {[4.0800e-10]}
    {'M4_C' }    {[4.0800e-10]}
    {'M5_C' }    {[4.0800e-10]}
    {'M6_C' }    {[4.0800e-10]}
    {'M7_C' }    {[4.0800e-10]}
    {'M8_C' }    {[4.0800e-10]}
    {'M9_C' }    {[4.0800e-10]}
    {'M10_C'}    {[4.0800e-10]}
    {'M11_C'}    {[4.0800e-10]}
    {'M12_C'}    {[4.0800e-10]}
    {'M13_C'}    {[4.0800e-10]}
    {'M14_C'}    {[4.0800e-10]}
    {'M15_C'}    {[4.0800e-10]}
    {'M16_C'}    {[4.0800e-10]}
    {'M1_R' }    {[    0.0036]}
    {'M2_R' }    {[    0.0036]}
    {'M3_R' }    {[    0.0036]}
    {'M4_R' }    {[    0.0036]}
    {'M5_R' }    {[    0.0036]}
    {'M6_R' }    {[    0.0036]}
    {'M7_R' }    {[    0.0036]}
    {'M8_R' }    {[    0.0036]}
    {'M9_R' }    {[    0.0036]}
    {'M10_R'}    {[    0.0036]}
    {'M11_R'}    {[    0.0036]}
    {'M12_R'}    {[    0.0036]}
    {'M13_R'}    {[    0.0036]}
    {'M14_R'}    {[    0.0036]}
    {'M15_R'}    {[    0.0036]}
    {'M16_R'}    {[    0.0036]}];

ts_Baxter = [0.25 0.005 0.25 0.005 0.25 0.005 0.25 0.005];
ts = (ts_Baxter./sum(ts_Baxter)).*Ts;

% The inital guess of time intervals % The inital guess of time intervals
% Assigned later dynamically


% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:


% List out all char variables in the
Switch_Resistors = {
    'M1_R'
    'M2_R'
    'M3_R'
    'M4_R'
    'M5_R'
    'M6_R'
    'M7_R'
    'M8_R'
    'M9_R'
    'M10_R'
    'M11_R'
    'M12_R'
    'M13_R'
    'M14_R'
    'M15_R'
    'M16_R'
    };
% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly

Switch_Names = {
    'M1'
    'M2'
    'M3'
    'M4'
    'M5'
    'M6'
    'M7'
    'M8'
    'M9'
    'M10'
    'M11'
    'M12'
    'M13'
    'M14'
    'M15'
    'M16'
    };



ON = 1;
OFF = 0;

ron = [0.0036 0.0036 0.0036 0.0036 0.0036 0.0036 0.0036 0.0036];

% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,16).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [ron(1) ron(2) ron(3) ron(2) ron(1) ron(3) ron(4) ron(5) ron(5) ron(4) ron(6) ron(7) ron(8) ron(7) ron(6) ron(8)];
SW = [SW_OFF;SW_ON;SW_ON];

Diode_Forward_Voltage = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'.*1.5;
u = [14.3 4 4 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]';
Diode_Forward_Voltage = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]'.*0;
u = [14.3 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
Order = [1 2 3 4 5 6 7 8];

%}


%% Load for 3 Level Buck

%{
ON = 1;
OFF = 0;
swvec = [
    ON  OFF ON  ON  OFF OFF   % Q1 ON
    OFF OFF ON  ON  OFF OFF  % Dead Q1 Q2
    OFF ON  ON  ON  OFF OFF  % Q2 ON
    OFF OFF ON  ON  OFF OFF];  % Dead


Ts = 0.5e-6;
Numerical_Components = [
    {'C1'  }    {[2.3500e-05]}
    {'C2'  }    {[2.0000e-06]}
    {'C3'  }    {[2.0000e-06]}
    {'L1'  }    {[1.2978e-07]}
    {'M1_C'}    {[4.0800e-10]}
    {'M2_C'}    {[4.0800e-10]}
    {'M3_C'}    {[4.0800e-10]}
    {'M4_C'}    {[4.0800e-10]}
    {'M5_C'}    {[1.5300e-09]}
    {'M6_C'}    {[1.5300e-09]}
    {'R1'  }    {[    0.0015]}
    {'R2'  }    {[    0.0050]}
    {'R3'  }    {[    0.0050]}
    {'R4'  }    {[    0.0020]}
    {'R5'  }    {[    0.0020]}
    {'R6'  }    {[    0.0022]}
    {'R7'  }    {[    0.0017]}
    {'R8'  }    {[   66.4000]}
    {'R9'  }    {[  304.4797]}
    {'C4'  }    {[5.6200e-12]}
    {'R10' }    {[1.0000e-03]}
    {'M1_R'}    {[    0.0036]}
    {'M2_R'}    {[    0.0036]}
    {'M3_R'}    {[    0.0036]}
    {'M4_R'}    {[    0.0036]}
    {'M5_R'}    {[    0.0160]}
    {'M6_R'}    {[    0.0160]}];

ts_Baxter = [0.70 0.001 0.3 0.001];
ts = (ts_Baxter./sum(ts_Baxter)).*Ts;

% The inital guess of time intervals % The inital guess of time intervals
% Assigned later dynamically


% List all of the numerical components in the netlist file for all
% FETs you must use the syntax used below:


% List out all char variables in the
Switch_Resistors = {
    'M1_R'
    'M2_R'
    'M3_R'
    'M4_R'
    'M5_R'
    'M6_R'
    };
% List of the switch sequency. Organized by: the FETs (column) vs time
% interval (rows) matching Switch_Resistors and ts respectivly

Switch_Names = {
    'M1'
    'M2'
    'M3'
    'M4'
    'M5'
    'M6'
    };



ON = 1;
OFF = 0;

ron = [0.0036 0.0036 0.0036 0.0036 0.0160 0.0160];

% List all the resistances of the diodes or FETS when they are on or
% off
SW_OFF = ones(1,6).*10000000;

%SW_ON = [M1_R,M2_R,etc...]
SW_ON = [ron(1) ron(2) ron(3) ron(4) ron(5) ron(6)];
SW = [SW_OFF;SW_ON;SW_ON];

Diode_Forward_Voltage = [1 1 1 1 1 1]'.*1.5;
u = [12 4 4 1.5 1.5 1.5 1.5 1.5 1.5]';
Order = [1 2 3 4];
%}
%% Analyze circuit

sim = SMPSim();
conv = sim.converter;
top = sim.topology;

% Set to run numerical parsesr

%conv.ts = ts;
conv.u = u;
top.order = Order;
top.Element_Properties = Numerical_Components;
top.Switch_Resistors = Switch_Resistors;
top.Switch_Resistor_Values = SW;
top.Switch_Sequence = swvec;
top.Fwd_Voltage = Diode_Forward_Voltage;
top.Switch_Names = Switch_Names;
%

top.loadCircuit(modelfile,swvec,1);



sim.u = u';
conv.setSwitchingPattern(1:size(swvec,1), ts)

Xss = sim.steadyState;
if(debug)
    sim.plotAllStates(1);
end

% ssOrder = plecs('get', circuitPath, 'StateSpaceOrder');
% outputs = ssOrder.Outputs;

outputs = top.outputLabels;

%% finalRun
% once everything seems to be error-free based on discrete time points,
% goes through once more with eigenvalue-based spacing to make sure no
% inter-sample violations are occuring.
finalRun = 0;

%% Symmetry check
% May be useful but not doing anything with it yet.  Can identify that DAB,
% etc. exhibit half-cycle symmetry
% TF = conv.checkForSymmetry;

%%
[Xf,ts,swinds] = timeSteppingPeriod(sim);

tic
while(1)
    Xss = sim.steadyState;
    
    %% Update constraints per the current switching vector
    %     [Cbnd, Dbnd, hyst, switchRef] = top.getPLECSConstraintMatrices(circuitPath);
    Cbnd = top.Cbnd; Dbnd = top.Dbnd; hyst = top.bndHyst; switchRef = top.switchRef;
    Cbnd = Cbnd(:,:,conv.swind);
    Dbnd = Dbnd(:,:,conv.swind);
    
    if(0)
        %% ~~Continuous time violation margin
        [ xs, t] = sim.SS_WF_Reconstruct;
        for i = 1:length(xs)-1
            swstate = find(t(i) <= cumsum(ts),1,'first');
            violationMargin(:,i) = Cbnd(:,:,swstate)*xs(:,i) + Dbnd(:,:,swstate)*us' - hyst(:,1) + hyst(:,2);
        end
        
        figure; plot(t(1:end-1),violationMargin)
        ylim([-100 100])
        legend(sim.switchNames) %% <- problem here, now.
    end
    
    %% Discrete timepoint violation margin
    % %     violateMarginStart = zeros(size(Cbnd,1) ,length(conv.swind));
    % %     targetValStart = violateMarginStart;
    % %     violateMarginEnd = zeros(size(Cbnd,1) ,length(conv.swind));
    % %     targetValEnd = violateMarginEnd;
    % %     for i = 1:length(conv.swind)
    % %         violateMarginStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us' - hyst(:,1) + hyst(:,2);
    % %         targetValStart(:,i) = Cbnd(:,:,i)*Xss(:,i) + Dbnd(:,:,i)*us' - hyst(:,1);
    % %         violateMarginEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us' - hyst(:,1) + hyst(:,2);
    % %         targetValEnd(:,i) = Cbnd(:,:,i)*Xss(:,i+1) + Dbnd(:,:,i)*us' - hyst(:,1);
    % %     end
    [violateMarginStart,violateMarginEnd,targetValStart,targetValEnd] = sim.checkDiscreteErr;
    errBefore = min(violateMarginStart,0);
    errAfter = min(violateMarginEnd,0);
    
    %% Check if we have the correct state to change times
    
    % if yes, proceed to jacobian.
    % if no, add states, then redo steady-state & error detection.
    % OR -- just do it anyway, let addUncontrolledSwitching report on whether
    % it actually did anything.
    
    %% Insert additional states as necessary
    % We only need to insert a new state if
    %   -- The interface has an error both before and after the switching
    %   -- OR it has one of the above, and the before and after switching
    %   positions aren't part of the modifiable uncontrolled times.
    
    [~,ints,~] = getIntervalts(conv);
    ints = ints';
    %     if size(ints,1) > size(ints,2)
    %         ints = ints';
    %     end
    
    %     errLocs = (errAfter<0 & circshift(errBefore,-1,2)<0);
    %     unadjustableLocs = (errAfter<0 | circshift(errBefore,-1,2)<0) & (ints ~= circshift(ints,-1));
    %     tLocs = errLocs | unadjustableLocs;
    %     insertAfter = any(tLocs,1);
    % The above only worked for always inserting AFTER
    % Note, this was to not inser again when the error was addressable by
    % the altready
    
    %     errLocs = (errAfter<0 & (errBefore<0));
    %     unadjustableLocs = (errAfter<0  & (ints ~= circshift(ints,-1))) | ...
    %         (errBefore<0 & (ints ~= circshift(ints,1)));
    %     adjustableButSameInterval = (errAfter<0  & ...
    %         ((ints == circshift(ints,-1)) & (conv.swind == circshift(conv.swind,-1)))) | ...
    %         (errBefore<0 & ...
    %         ((ints == circshift(ints,1)) & (conv.swind == circshift(conv.swind,1))));
    %     tLocs = errLocs | unadjustableLocs | adjustableButSameInterval;
    %     insertAt = any(tLocs,1);
    %
    %     adjType = zeros(size(Cbnd,1), length(conv.ts), 2);
    %     adjType(:,:,1) = (errAfter<0)+0;
    %     adjType(:,:,2) = (-1)*(errBefore<0);
    
    [tLocs,insertAt,adjType] = sim.findRequiredUncontrolledSwitching(violateMarginStart,violateMarginEnd);
    
    altered = 0;
    
    
    %     warning("Need something so it doesn't do corrections both before and after");
    for i = flip(find(insertAt))
        % addUncontrolledSwitching(obj, interval, beforeAfter, initialTime, switches, newStates)
        dt = min( min(conv.controlledts)/20 , conv.ts(i)/20);
        for j = 1:2
            if(any(adjType(:,i,j)))
                altered = altered | conv.addUncontrolledSwitching(i,(-1)^(j+1),dt,switchRef(tLocs(:,i),1),~switchRef(tLocs(:,i),2));
            end
        end
    end
    
    Xss = sim.steadyState;
    if(debug)
        sim.plotAllStates(1);
    end
    
    
    
    if ~altered
        
        if ~any(any(errBefore) | any(errAfter)) && finalRun == 0
            finalRun = 1;
            eigs2tis(conv);
            Xss = sim.steadyState;
            continue;
        elseif ~any(any(errBefore) | any(errAfter)) && finalRun == 1
            break;
        else
            finalRun = 0;
        end
        %% Correct based on discrete jacobian
        
        [Jt, J2t] = discreteJacobian(sim, 2);
        
        %% Check validity of jacobian function (test)
        if(debug2)
            [Jo] = discreteJacobian(sim, 1);
            func = @(x) getSSforJacobian(sim, x);
            [jac,err] = jacobianest(func,zeros(size(conv.ts)));
            jac = reshape(jac, size(Jo))*1e12;
            Jo
            jac./Jo
            %                 Jt = -jac
        end
        
        Jout = zeros(size(Cbnd,1),size(Jt,2),size(Jt,3));
        JoutStart = Jout;
        JoutEnd = Jout;
        %J(state, at time, time changed)
        for i = 1:size(Jt,3)
            for j = 1:size(Jt,2)
                Jout(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);
                
                % Cbnd is during the interval, but Jt is about the states at the
                % interface.
                JoutStart(:,j,i) = Cbnd(:,:,j)*Jt(:,j,i);
                endInt = circshift(1:size(Jt,2),-1);
                JoutEnd(:,j,i) = Cbnd(:,:,j)*Jt(:,endInt(j),i);
            end
        end
        
        
        %% Works with only one violation:
        % J = squeeze(Jout(:,intV,intV-1));
        % Err = violateMarginEnd
        % dt = -J\Err(:,intV);
        
        % % %         A = zeros(length(conv.ts),length(conv.ts));
        % % %         b = zeros(length(conv.ts),1);
        % % %         tindex = 1:length(conv.ts);
        % % %
        % % %         rows2elim = [];
        % % %         cols2elim = [];
        % % %
        % % %         for i = 1:length(j1)
        % % %             deltas = squeeze(JoutEnd(i1(i),j1(i),:));
        % % %             A(i,:) = deltas';
        % % %             b(i,1) = violateMarginEnd(i1(i),j1(i));
        % % %         end
        % % %         if(isempty(i)) i=0; end
        % % %         for j = 1:length(j2)
        % % %             deltas = squeeze(JoutStart(i2(j),j2(j),:));
        % % %             A(i+j,:) = deltas';
        % % %             b(i+j,1) = violateMarginStart(i2(j),j2(j));
        % % %         end
        % % %         if(isempty(j)) j=0; end
        % % %         for k = 1:(length(conv.ts) - (length(j1)+length(j2)))
        % % %             candidateSkips = setdiff(1:length(conv.ts), ...
        % % %                 [mod(j2-2,length(conv.ts))+1; j1]);
        % % %
        % % %             rows2elim = [rows2elim, i+j+k];  %Equations
        % % %             cols2elim = [cols2elim, candidateSkips(k)]; % Times
        % % %
        % % %
        % % %         %     A(i+j+k,:) = -NaN;
        % % %         %     A(:,candidateSkips(k)) = -NaN;
        % % %         %     b(i+j+k,:) = NaN;
        % % %         %     tindex(candidateSkips(k)) = NaN;
        % % %         end
        % % %
        % % %         % A(isnan(A)) = [];
        % % %         % b(isnan(b)) = [];
        % % %         % tindex(isnan(tindex)) = [];
        % % %
        % % %         A(rows2elim,:) = [];
        % % %         A(:,cols2elim) = [];
        % % %         b(rows2elim) = [];
        % % %         tindex(cols2elim) = [];
        % % %
        % % %         % for i = 1:length(intV)
        % % %         %     deltas = squeeze(Jout(stateV(i),intV(i),:));
        % % %         %     A(i,:) = deltas';
        % % %         %     b(i,1) = min(violateMarginEnd(stateV(i),intV(i)), violateMarginStart(stateV(i),intV(i)));
        % % %         % end
        % % %         % for j = [length(intV):length(ts)] - (length(intV)-1)
        % % %         %     candidateSkips = setdiff(1:length(ts), mod(intV-2,length(ts))+1);
        % % %         %     A(i+j,candidateSkips(j)) = 1;
        % % %         %     b(i+j,1) = 0;
        % % %         % end
        
        % % % %         [i1,j1] = find(violateMarginEnd < 0);
        % % % %         [i2,j2] = find(violateMarginStart < 0);
        % % % %         intV = [j1; j2]; stateV = [i1; i2];
        % % % %
        % % % %         A = zeros(length(j1)+length(j2),length(conv.ts));
        % % % %         b = zeros(length(j1)+length(j2),1);
        % % % %         tindex = 1:length(conv.ts);
        % % % %
        % % % %         for i = 1:length(j1)
        % % % %             deltas = squeeze(JoutEnd(i1(i),j1(i),:));
        % % % %             A(i,:) = deltas';
        % % % %             b(i,1) = targetValEnd(i1(i),j1(i));
        % % % %         end
        % % % %         if(isempty(i)), i=0; end
        % % % %         for j = 1:length(j2)
        % % % %             deltas = squeeze(JoutStart(i2(j),j2(j),:));
        % % % %             A(i+j,:) = deltas';
        % % % %             b(i+j,1) = targetValStart(i2(j),j2(j));
        % % % %         end
        % % % %         if(isempty(j)), j=0; end
        
        [i1,j1] = find(errBefore);
        [i2,j2] = find(errAfter);
        intV = [j1; j2]; stateV = [i1; i2];
        
        A = zeros(length(j1)+length(j2),length(conv.ts));
        b = zeros(length(j1)+length(j2),1);
        tindex = 1:length(conv.ts);
        
        for i = 1:length(j1)
            deltas = squeeze(JoutStart(i1(i),j1(i),:));
            A(i,:) = deltas';
            b(i,1) = targetValStart(i1(i),j1(i));
        end
        if(isempty(i)), i=0; end
        for j = 1:length(j2)
            deltas = squeeze(JoutEnd(i2(j),j2(j),:));
            A(i+j,:) = deltas';
            b(i+j,1) = targetValEnd(i2(j),j2(j));
        end
        if(isempty(j)), j=0; end
        
        unChangeable = isnan(sum(A,1));
        
        
        tsolve = zeros(size(conv.ts));
        
        %% Attempt: add zero net perturbation to time as a part of the
        %% equations -- MAY NOT  because some times are dropped.
        scaleF = norm(A(:,~unChangeable))/numel(A);
        [~, timeInts, ~] = conv.getIntervalts;
        A = [A; zeros(max(timeInts), size(A,2))];
        b = [b; zeros(max(timeInts), 1)];
        for i = 1:max(timeInts)
            if ~any(timeInts'==i & unChangeable)
                A(end-max(timeInts)+i,timeInts==i) = scaleF;
            end
        end
        
        A(:,unChangeable) = [];
        
        tsolve(~unChangeable) = -(A\b);
        
        if any(isnan(tsolve))
            tsolve(~unChangeable) = -pinv(A)*b;
        end
        
        [~,dtLims] = getDeltaT(sim.converter);
        
        tr = tsolve./dtLims;
        tr = max(tr(tr>1));
        
        if ~isempty(tr)
            tsolve = tsolve/tr;
        end
        
        %         deltaTs = zeros(size(conv.ts));
        %         deltaTs(tindex) = tsolve;
        %
        % %% for MRbuck initial iterations:
        % -violateMarginEnd(2,3)/Jt(3,4,3)
        
        %%HDSC only allow deadtimes to be adjusted
        %         A(:,1) = [];
        %         A(:,2) = [];
        %         deltaDTs = -A\b
        %         deltaTs = [0 deltaDTs(1) 0  deltaDTs(2)];
        
        oldts = conv.ts;
        for i=length(tsolve):-1:1
            if tsolve(i) ~= 0
                %                conv.adjustTiming(i, deltaTs(i));
                conv.adjustUncontrolledTiming(i, tsolve(i));
            end
        end
        % % %         conv.adjustUncontrolledTimingVector(1:length(tsolve), tsolve)
        newts = conv.ts;
        
        %         warning('Can get stuck trying to make intervals longer when it cannot');
        
        
        % for i = 1:length(intV)
        %     int = intV(i);
        %     J = squeeze(Jout(:,int,int-1));
        %     Err = violateMarginEnd;  %% THIS ISN't RIGHT, DEPENDS ON WHERE IT OCCURS
        %     dt(i) = -J\Err(:,int);
        %     conv.adjustTiming(intV(i), dt(i));
        % end
        
        
        % % Cbnd(:,:,3)*Jt(:,4,3)
        % % dt = -  Cbnd(:,:,3)*Xss(:,4) [2]   /  Cbnd(:,:,3)*Jt(:,4,3) [2]
        % %     =   -Err  /  Jacobian (assuming a state to correct from)
        % sim.ts(3) = sim.ts(3) + dt;
        % sim.ts(4) = sim.ts(4) - dt;
        
        if(debug)
            %             disp(tsolve)
            %             disp(sum(errBefore + errAfter, 'all'))
            Xss = sim.steadyState;
            sim.plotAllStates(10);
        end
        
        niter = niter+1;
        disp([niter sum(errBefore + errAfter, 'all')]);
        
        if(~any(tsolve))
            error('timing not modified');
            break;
        end
        
    end
end
toc

Xss = sim.steadyState;
sim.plotAllStates(10);

function Xs = getSSforJacobian(sim, newTs)
[tps] = sim.converter.validateTimePerturbations([1:length(newTs)], newTs/1e12);
Xs = sim.perturbedSteadyState(tps);
Xs = Xs(:,1:end-1);
end
