function [check] = VfwdIrev(obj,Xs,u,order)
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

ron = .05;
Anum = obj.Anum;
Bnum = obj.Bnum;


% Variable Declaration:

Max_Diode_Reverse_Current = -0.005; % Negative for reverse current
Max_Diode_Forward_Voltage = 0.7; % Maximum Forward voltage of diode
Max_FET_Forward_Voltage = 0.5; % Maximum Forward voltage of diode in FET

% Calculate the response of each time interval of the circuit
for i = 1:1:size(Xs,2)-1
Y(:,:,i) = obj.Cnum(:,:,order(i))*Xs(:,i+1)+obj.Dnum(:,:,order(i))*u;
end

% Get measurement names and type of measurement Voltage (V) or Current (A)
[measurename,remain] = strtok(obj.OutputNamesCD(:,1));


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
            fprintf('Reverse current violation of Diode %s during time interval %.0f \n', violations{j},i)
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
    violations = measurename(Diode_voltage(:,:,i)>Max_Diode_Forward_Voltage,1);
    if ~isempty(violations)
        for j = 1:1:length(violations)
            fprintf('Forward voltage violation of Diode %s during time interval %.0f \n', violations{j},i)
        end
    end
end

%% Diode Voltage

% For a FET there can be a reverse current if the FET is ON (will put in later commit)

% Find FET voltage measurement positions in C and D matrix
FET_voltage_pos = contains(measurename,'M') & contains(remain, ' V');

% Get FET voltage measurements for each time interval
FET_voltage = Y.*FET_voltage_pos;

% Iterate through and check violation of forward voltage in FET body diode
for i = 1:1:size(FET_voltage,3)
    violations = measurename(FET_voltage(:,:,i)<-Max_FET_Forward_Voltage,1);
    if ~isempty(violations)
        for j = 1:1:length(violations)
            fprintf('Forward voltage violation of FET'' diode %s during time interval %.0f \n', violations{j},i)
        end
    end
end


%% Check to see if entering other state than the one listed in order

% Attempt at identifing fixing hard switching

ONorOFF = obj.ONorOFF;

for i = 1:1:size(Xs,1) % Cycle through state variables  
    for j = 2:1:size(Xs,2) % Cycle through time intervals 
        j;% is time interval for Xss
        k = j-1; % k is time interval for everything else
        if ONorOFF(i,k) ~=0 % if FET or Diode
            
            Current_StateName = strsplit(obj.StateNames{i,1},'_'); % Separte
            Voltage_pos = strcmp(repmat(Current_StateName{1},10,1),measurename) & contains(remain, ' V'); % Find position of selected state variable voltage measurement
            Current_pos = strcmp(repmat(Current_StateName{1},10,1),measurename) & contains(remain, ' A'); % Find position of selected state variable current measurement
            Voltage = Y.*Voltage_pos; % Set Voltage of switching element
            Current = Y.*Current_pos; % Set Current of swithcing element
            
            if k == 1
                DeltaV = Voltage(:,:,k)-Voltage(:,:,end);
                DeltaC = Current(:,:,k)-Current(:,:,end);
                StartDeriv=Anum(i,:,k)*Xs(:,j-1)+Bnum(i,:,k)*u;
                EndDeriv=Anum(i,:,k)*Xs(:,j)+Bnum(i,:,k)*u;
            else
                DeltaV = Voltage(:,:,k)-Voltage(:,:,k-1);
                DeltaC = Current(:,:,k)-Current(:,:,k-1);
                StartDeriv=Anum(i,:,k)*Xs(:,j-1)+Bnum(i,:,k)*u;
                EndDeriv=Anum(i,:,k)*Xs(:,j)+Bnum(i,:,k)*u;
            end
            
            % Initally set start time equal to zero
            time = 0;
            data = [];
            %%% Attempt to find local max or min of voltage for diode %%%
            
            for n = 100:1:200
            dxdt = Anum(:,:,k)*expm(Anum(:,:,k)*time)*Xs(:,j)+expm(Anum(:,:,k)*time)*Bnum(:,:,k)*u;
            time = time+1/n*dxdt(i,1);
            data(end+1) = time;
            end
            
            
            
            if obj.DMpos(i,2)==1 % if diode
                
                if ONorOFF(i,k) == 1 % if diode ON
                    if sum(Voltage(:,:,k)) < 0
                        fprintf('Violation of %s in state %.0f \n',obj.StateNames{i,1},j-1)

                    end
                elseif ONorOFF(i,j-1) == -1 % if diode off
                    if sum(Current(:,:,k)) > 0
                        fprintf('Violation of %s in state %.0f \n',obj.StateNames{i,1},j-1)
                    end
                else
                    fprintf('Messed up\n')
                end
            elseif obj.DMpos(i,3)==1 % if FET
                if ONorOFF(i,j-1) == 1 % if FET ON
                    if sum(Voltage(:,:,j-1)) < 0
                        % This doesn't matter if the FET is on
                        % Is already checked above
                        % Can have bidirectional voltage and current
                        % fprintf('Violation of %s in state %.0f \n',obj.StateNames{i,1},j-1)
                    end
         
                elseif ONorOFF(i,j-1) == -1 % if FET off
                    if sum(Current(:,:,j-1)) > 0
                        fprintf('Violation of %s in state %.0f \n',obj.StateNames{i,1},j-1)
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Find if FET turns on next
                    if size(ONorOFF,2)+1 == j
                        if ONorOFF(i,2) == 1
                            if Voltage < -10*ron*2
                                fprintf('Hard swithcing for %s in time interval j\n')
                            end
                        end
     
                    else
                        if ONorOFF(i,j) == 1
                            if Voltage < -10*ron*2
                                fprintf('Hard swithcing for %s in time interval j\n')
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    fprintf('Messed up\n')
                end
            else
                fprintf('Found a FET or Diode that wasn''t a FET or Diode\n')
            end

        end

    end
end




end
