function [check] = VfwdIrev(obj)
%VfwdIrev This function checks to ensure the physical limitations of diodes
%are met in the circuit.
%   Detailed explanation goes here

% Determine if there is a violation of Reverse current or forward voltage
% in a diode:

% "For an ideal diode there should never be a reverse current." - Physics

% Remeber Representation of FET:

% Postive voltage and current from Drain to Source
%%% This is Cathode to Anode for the body diode!!!!

% Positive voltage and current from Anode to Cathode for Power diodes

%%%% Thus FET and Diode voltage and current are switched for diodes!!!!!

% For a FET there can be a reverse current if the FET is ON
check = 0;
ron = .05;
As = obj.As;
Bs = obj.Bs;
debug = 1;
Xs = obj.Xs;
u = obj.u;

% Variable Declaration:
Max_Diode_Reverse_Current = -0.005; % Negative for reverse current
Max_Diode_Forward_Voltage = 0.7; % Maximum Forward voltage of diode
Max_FET_Forward_Voltage = 0.5; % Maximum Forward voltage of diode in FET
Max_FET_Reverse_Current = 0.005; % Positive due to polarity of current measurement for FET diode

% Calculate the response of each time interval of the circuit
for i = 1:1:size(Xs,2)-1
    Y(:,:,i) = obj.Cs(:,:,i)*Xs(:,i+1)+obj.Ds(:,:,i)*u;
end

% Get measurement names and type of measurement Voltage (V) or Current (A)
[measurename,remain] = strtok(obj.Converter.Topology.Parser.OutputNamesCD(:,1));


%% Diode Current
% Find diode current measurement positions in C and D matrix
diode_current_pos = contains(measurename,'D') & contains(remain, ' A');

% Get diode current measurements for each time interval
Diode_current = Y.*diode_current_pos;

% Iterate through and check violation of reverse current in diode
for i = 1:1:size(Diode_current,3)
    violations = measurename(Diode_current(:,:,i)<Max_Diode_Reverse_Current,1);
    if ~isempty(violations)
        for j = 1:1:length(violations)
            fprintf('Reverse current violation of Diode %s at end of time interval %.0f \n', violations{j},i)
        end
    end
end


%% Diode Voltage

% Find diode voltage measurement positions in C and D matrix
diode_voltage_pos = contains(measurename,'D') & contains(remain, ' V');

% Get diode voltage measurements for each time interval
Diode_voltage = Y.*diode_voltage_pos;

% Iterate through and check violation of forward voltage in diode
for i = 1:1:size(Diode_voltage,3)
    violations = measurename(Diode_voltage(:,:,i)>Max_Diode_Forward_Voltage,1); % Positive because diode voltage is reported at anode to cathode
    if ~isempty(violations)
        for j = 1:1:length(violations)
            fprintf('Forward voltage violation of Diode %s during time interval %.0f \n', violations{j},i)
        end
    end
end

%% FET Voltage

% Find FET voltage measurement positions in C and D matrix
FET_voltage_pos = contains(measurename,'M') & contains(remain, ' V');

% Get FET voltage measurements for each time interval
FET_voltage = Y.*FET_voltage_pos;

% Iterate through and check violation of forward voltage in FET body diode
for i = 1:1:size(FET_voltage,3)
    violations = measurename(FET_voltage(:,:,i)<-Max_FET_Forward_Voltage,1); % Because FET Voltage is given from drain to source (Cathode to Anode)
    if ~isempty(violations)
        for j = 1:1:length(violations)
            fprintf('Forward voltage violation of FET''s diode %s during time interval %.0f \n', violations{j},i)
        end
    end
end


%% FET Current
ONorOFF = obj.Converter.Topology.Parser.ONorOFF;

% Find FET voltage measurement positions in C and D matrix
FET_current_pos = contains(measurename,'M') & contains(remain, ' A');

% Get FET voltage measurements for each time interval
FET_current = Y.*FET_current_pos;

% Iterate through and check violation of forward voltage in FET body diode

FET_Names = obj.Converter.Topology.Parser.StateNames(obj.Converter.Topology.Parser.DMpos(:,3)==1,1); % Get the name of FETs in the converter
[FET_names,~] = strtok(FET_Names,'_'); % Separate FET names from the '_C'
m = 0;
for k = 1:1:size(obj.Converter.Topology.Parser.DMpos,1) % Iterate through all state variables
    if obj.Converter.Topology.Parser.DMpos(k,3)==1 % If the index is a FET
        m = m+1; % Increment variable
        for i = 1:1:size(ONorOFF,2) % Iterate through all time intervals
            if ONorOFF(k,i) == -1 % If the FET is off during given time interval
                comparison = repmat(FET_names(m),size(measurename,1),1); % Use repmat to create copy of FET names that is equal in length to output names (Y)
                FET_current_pos = strcmp(measurename,comparison) & contains(remain, ' A'); % Compare chars to determine where the current measurement is for the given FET
                FET_current = Y.*FET_current_pos; % Solve for the measured value of the FET current
                violations = measurename(FET_current(:,:,i)>Max_FET_Reverse_Current,1); % Because FET Current is given from drain to source (Cathode to Anode)
                if ~isempty(violations)
                    for j = 1:1:length(violations)
                        fprintf('Reverse current violation of FET''s diode %s during time interval %.0f \n', violations{j},i)
                    end
                end
            end
            
        end
        
    end
end


%% Check to see if entering other state than the one listed in order

% Attempt at identifing fixing hard switching



global  I_solemnly_swear_that_I_am_up_to_no_good;

for i = 1:1:size(Xs,1) % Cycle through state variables
    for j = 2:1:size(Xs,2) % Cycle through time intervals
        j;% is time interval for Xss
        k = j-1; % k is time interval for everything else
        if ONorOFF(i,k) ~=0 % if FET or Diode
            
            Current_StateName = strsplit(obj.Converter.Topology.Parser.StateNames{i,1},'_'); % Separte
            Voltage_pos = strcmp(repmat(Current_StateName{1},10,1),measurename) & contains(remain, ' V'); % Find position of selected state variable voltage measurement
            Current_pos = strcmp(repmat(Current_StateName{1},10,1),measurename) & contains(remain, ' A'); % Find position of selected state variable current measurement
            Voltage = Y.*Voltage_pos; % Set Voltage of switching element
            Current = Y.*Current_pos; % Set Current of switching element
            
            if k == 1
                DeltaV = Voltage(:,:,k)-Voltage(:,:,end);
                DeltaC = Current(:,:,k)-Current(:,:,end);
                StartDeriv=As(i,:,k)*Xs(:,j-1)+Bs(i,:,k)*u;
                EndDeriv=As(i,:,k)*Xs(:,j)+Bs(i,:,k)*u;
            else
                DeltaV = Voltage(:,:,k)-Voltage(:,:,k-1);
                DeltaC = Current(:,:,k)-Current(:,:,k-1);
                StartDeriv=As(i,:,k)*Xs(:,j-1)+Bs(i,:,k)*u;
                EndDeriv=As(i,:,k)*Xs(:,j)+Bs(i,:,k)*u;
            end
            
            J = zeros(size(Xs,1),size(Xs,1));
            J(i,i) = 1;
            %  fun = @(t)((As(i,:,k)*expm(As(:,:,k).*t)*Xs(:,j)+expm(As(:,:,k).*t)*Bs(i,:,k)*u));
            
            I_solemnly_swear_that_I_am_up_to_no_good = [i,j,k];
            x0 = 4e-6;
            A = [];
            B = [];
            Aeq = [];
            Beq = [];
            LB = 0;
            UB = 5e-6;
            fun = @(t) t;
            
            options = optimoptions('fmincon','StepTolerance',1e-12,'ConstraintTolerance',1e-9);
            if debug
                options = optimoptions('fmincon','Display','off','StepTolerance',1e-65,'ConstraintTolerance',1e-3,'OptimalityTolerance',1e-5);
            end
            [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(fun,x0,A,B,Aeq,Beq,LB,UB,@obj.zero,options);
            
            % Plot dx/dt if debug
            if debug
                J = zeros(size(Xs,1),size(Xs,1));
                J(i,i) = 1;
                syms t
                stack = [];
                t1 = linspace(0,5e-6,1000);
                for counts = 1:1:length(t1)
                    eqn2 = J*(As(:,:,k)*expm(As(:,:,k).*t1(counts))*Xs(:,j)+expm(As(:,:,k).*t1(counts))*Bs(:,:,k)*u);
                    stack(end+1) = eqn2(1,:);
                end
                plot(t1,stack)
            end
            %{
            % Create symbolic equation 0 = Ax+Bu (Not the best way to do this)
            eqn = zeros(size(Xs,1),1) == J*(As(:,:,k)*expm(As(:,:,k).*t)*Xs(:,j)+expm(As(:,:,k).*t)*Bs(:,:,k)*u);
            the_answer=vpasolve(eqn,t,[0 5e-9]); % numerically solve in given time interval
            disp(the_answer)

            clear the_answer
            %}
            
            %%% Attempt to find local max or min of voltage for diode %%%
            %  time = 0;
            %  data = [];
            %             for n = 100:1:200
            %             dxdt = As(:,:,k)*expm(As(:,:,k)*time)*Xs(:,j)+expm(As(:,:,k)*time)*Bs(:,:,k)*u;
            %             time = time+1/n*dxdt(i,1);
            %             data(end+1) = time;
            %             end
            %end
            
            if obj.Converter.Topology.Parser.DMpos(i,2)==1 % if diode
                
                if ONorOFF(i,k) == 1 % if diode ON
                    if sum(Voltage(:,:,k)) < 0
                        fprintf('Violation of %s in state %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                        
                    end
                elseif ONorOFF(i,j-1) == -1 % if diode off
                    if sum(Current(:,:,k)) > 0
                        fprintf('Violation of %s in state %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                    end
                else
                    fprintf('Messed up\n')
                end
            elseif obj.Converter.Topology.Parser.DMpos(i,3)==1 % if FET
                if ONorOFF(i,j-1) == 1 % if FET ON
                    if sum(Voltage(:,:,j-1)) < 0
                        % This doesn't matter if the FET is on
                        % Is already checked above
                        % Can have bidirectional voltage and current
                        % fprintf('Violation of %s in state %.0f \n',obj.StateNames{i,1},j-1)
                    end
                    
                elseif ONorOFF(i,j-1) == -1 % if FET off
                    if sum(Current(:,:,j-1)) > 0
                        fprintf('Violation of %s in state %.0f \n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Find if FET turns on next
                    if size(ONorOFF,2)+1 == j
                        if ONorOFF(i,2) == 1
                            if Voltage < -10*ron*2
                                fprintf('Hard swithcing for %s in time interval j\n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                            end
                        end
                        
                    else
                        if ONorOFF(i,j) == 1
                            if Voltage < -10*ron*2
                                fprintf('Hard swithcing for %s in time interval j\n',obj.Converter.Topology.Parser.StateNames{i,1},j-1)
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    fprintf('Messed up\n')
                end
            else
                fprintf('Found a FET or Diode that wasn''t a FET or Diode\n Don''t panic\n')
            end
            
        end
        
    end
end

% Mischief Managed
clear I_solemnly_swear_that_I_am_up_to_no_good

end % That's all Folks
