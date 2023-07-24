function [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrix(obj,circuitPath,foreced_Refresh)
%get ConstraintMatrics get Cbnd and Dbnd from LTSpice model
%   Note that the PLECs model must have probes added to all switching
%   device currents and voltages (double-click -> assertions -> (+))

if nargin == 1
    circuitPath = 'A cool circuit';
end

if nargin < 3
    forceRefresh = 0;
end


switchNames = obj.switchLabels;
outputs = obj.outputLabels;

swvec = obj.swseq;

Cbnd = zeros(size(obj.Cs));
Dbnd = zeros(size(obj.Ds));
hyst = zeros(size(obj.Ds,1), 2);
switchRef = zeros(size(obj.Ds,1), 2);

if isempty(obj.devType) || forceRefresh
    % Only re-run this once per file
    obj.isFET = zeros(size(switchNames));
    obj.isDiode = obj.isFET;
    obj.Vf = obj.isFET;
    for i = 1:length(switchNames)
        Cell_to_extract = switchNames{i};
        obj.devType{i} = Cell_to_extract(1);
        obj.isFET(i) = strcmp(obj.devType(i), 'M');
        obj.isDiode(i) = strcmp(obj.devType(i), 'D');
        if obj.isFET(i)
            obj.Vf(i) = 0; % Hard coded Vf for now
        end
    end
end

% For finding parallel devices
combCD3 = reshape(cat(2, obj.Cs, obj.Ds), size(obj.Cs,1),...
    (size(obj.Cs,2) + size(obj.Ds,2))*size(obj.Cs,3) );



for i = 1:length(switchNames)
    
    switchSignals = strncmp(pad(switchNames{i},length(switchNames{i})+1),outputs,length(switchNames{i})+1);
    assert(sum(switchSignals) >= 2, ['Switching device ' switchNames{i} ' does not have probes attached to both device current and device voltage']);
    
    currents = endsWith(outputs,'A');
    voltages = endsWith(outputs,'V');
    
    devCurrent = find(switchSignals & currents);
    devVoltage = find(switchSignals & voltages);
    devType = cellfun(@(s)s(1:1),switchNames(i));
    
    
    %% Find Parallel devices
    %This is to make sure we don't get an error about a reverse
    %diode voltage when a parallel MOSFET turns on.  Although that
    %is still an error that needs to be fixed, ignoring it will help the
    %loop converge faster by making the correction in the previous
    %interval
    [~,index_A1,~] = intersect(combCD3,combCD3(devVoltage,:),'rows');
    [~,index_A2,~] = intersect(combCD3,-combCD3(devVoltage,:),'rows');
    
    findInds = union(index_A1, index_A2);
    switchInds = zeros(size(findInds));
    for j = 1:length(findInds)
        sameVolts = strcmp(switchNames, outputs{findInds(j)}(1:end-2));
        switchInds(j) = find(sameVolts)';
    end
    
    if obj.isDiode(i)
        %on-state current > 0
        Cbnd(devCurrent, :, swvec(:,i)==1) = obj.Cs(devCurrent, :, swvec(:,i)==1);
        Dbnd(devCurrent, :, swvec(:,i)==1) = obj.Ds(devCurrent, :, swvec(:,i)==1);
        hyst(devCurrent,:)=[0, 10e-6];
        switchRef(devCurrent,:) = [i, 1];
        
        
        %off-state -(voltage) > 0
        %            Cbnd(devVoltage, :, swvec(:,i)==0) = -obj.topology.Cs(devVoltage, :, swvec(:,i)==0);
        %            Dbnd(devVoltage, :, swvec(:,i)==0) = -obj.topology.Ds(devVoltage, :, swvec(:,i)==0);
        Cbnd(devVoltage, :, all(swvec(:,switchInds)==0,2)) = -obj.Cs(devVoltage, :, all(swvec(:,switchInds)==0,2));
        Dbnd(devVoltage, :, all(swvec(:,switchInds)==0,2)) = -obj.Ds(devVoltage, :, all(swvec(:,switchInds)==0,2));
        
        
        % Vf = evalin('base', plecs('get', [modelFile '/' devStr], 'Vf'));
        hyst(devVoltage,:)=[-1.5, 0.25];
        switchRef(devVoltage,:) = [i, 0];
        
        
    elseif obj.isFET(i) && ~(any(obj.isFET(switchInds)) && any(obj.isDiode(switchInds)))
        %off-state voltage > 0
        Cbnd(devVoltage, :, swvec(:,i)==0) = obj.Cs(devVoltage, :, swvec(:,i)==0);
        Dbnd(devVoltage, :, swvec(:,i)==0) = obj.Ds(devVoltage, :, swvec(:,i)==0);
        hyst(devVoltage,:)=[-1, 0.1];
        switchRef(devVoltage,:) = [i, 0];
        
        %off-state -(current) > 0
        %Note: this is currently useless -- current is == zero always
        %when switch is off.
        %            Cbnd(devCurrent, :, swvec(:,i)==0) = -obj.Cs(devCurrent, :, swvec(:,i)==0);
        %            Dbnd(devCurrent, :, swvec(:,i)==0) = -obj.Ds(devCurrent, :, swvec(:,i)==0);
        %            hyst(devCurrent,:)=[-1, .5];
        
    end
end


zeroRows = switchRef(:,1) == 0;
Cbnd(zeroRows,:,:) = [];
Dbnd(zeroRows,:,:) = [];
hyst(zeroRows,:) = [];
switchRef(zeroRows,:) = [];

% eliminate redundant constraints
combBnd = reshape(cat(2, Cbnd, Dbnd), size(Cbnd,1), (size(Cbnd,2) + size(Dbnd,2))*size(Cbnd,3) );
[~, ia, ~] = unique(combBnd,'rows');

Cbnd = Cbnd(ia,:,:);
Dbnd = Dbnd(ia,:,:);
hyst = hyst(ia,:);
switchRef = switchRef(ia,:);

obj.Cbnd = Cbnd;
obj.Dbnd = Dbnd;
obj.bndHyst = hyst;
obj.switchRef = switchRef;

end

