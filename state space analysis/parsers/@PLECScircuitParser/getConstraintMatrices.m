function [Cbnd, Dbnd, hyst, switchRef] = getConstraintMatrices(obj,circuitPath)
%getPLECSConstraintMatrics get Cbnd and Dbnd from PLECs model
%   Note that the PLECs model must have probes added to all switching
%   device currents and voltages (double-click -> assertions -> (+))

    if nargin == 1
        circuitPath = obj.sourcefn;
    end

    S = strsplit(circuitPath,'/');
    modelFile = S{1};

    ssOrder = plecs('get', circuitPath, 'StateSpaceOrder');
    switchNames = ssOrder.Switches;
    outputs = ssOrder.Outputs;
    
    swvec = obj.topology.swseq;
    
    Cbnd = zeros(size(obj.topology.Cs));
    Dbnd = zeros(size(obj.topology.Ds));
    hyst = zeros(size(obj.topology.Ds,1), 2);
    switchRef = zeros(size(obj.topology.Ds,1), 2);
    
    isFET = zeros(size(switchNames));
    isDiode = isFET;
    for i = 1:length(switchNames)
        devType = plecs('get',[modelFile '/' switchNames{i}], 'Type');
        isFET(i) = strcmp(devType, 'MosfetWithDiode') || strcmp(devType, 'Mosfet');
        isDiode(i) = strcmp(devType, 'Diode');
    end
    
    % For finding parallel devices
    combCD3 = reshape(cat(2, obj.topology.Cs, obj.topology.Ds), size(obj.topology.Cs,1),...
    	(size(obj.topology.Cs,2) + size(obj.topology.Ds,2))*size(obj.topology.Cs,3) );

    for i = 1:length(switchNames)
       switchSignals = strncmp(switchNames{i},outputs,length(switchNames{i}));
       assert(sum(switchSignals) >= 2, ['Switching device ' switchNames{i} ' does not have probes attached to both device current and device voltage']);

       currents = endsWith(outputs,'Am');
       voltages = endsWith(outputs,'Vm');

       devCurrent = find(switchSignals & currents);
       devVoltage = find(switchSignals & voltages);

       
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
           sameVolts = strfind(switchNames, outputs{findInds(j)}(1:end-3));
           switchInds(j) = find(~cellfun(@isempty,sameVolts))';
       end
       
       if isDiode(i)
           %on-state current > 0
           Cbnd(devCurrent, :, swvec(:,i)==1) = obj.topology.Cs(devCurrent, :, swvec(:,i)==1);
           Dbnd(devCurrent, :, swvec(:,i)==1) = obj.topology.Ds(devCurrent, :, swvec(:,i)==1);
           hyst(devCurrent,:)=[0, obj.diodeCurrentThreshold];
           switchRef(devCurrent,:) = [i, 1];
          
           %off-state -(voltage) > 0
%            Cbnd(devVoltage, :, swvec(:,i)==0) = -obj.topology.Cs(devVoltage, :, swvec(:,i)==0);
%            Dbnd(devVoltage, :, swvec(:,i)==0) = -obj.topology.Ds(devVoltage, :, swvec(:,i)==0);
           Cbnd(devVoltage, :, all(swvec(:,switchInds)==0,2)) = -obj.topology.Cs(devVoltage, :, all(swvec(:,switchInds)==0,2));
           Dbnd(devVoltage, :, all(swvec(:,switchInds)==0,2)) = -obj.topology.Ds(devVoltage, :, all(swvec(:,switchInds)==0,2));
           
           devStr = outputs{devVoltage};
           devStr = devStr(1:find(devStr==':',1)-1);
           Vf = evalin('base', plecs('get', [modelFile '/' devStr], 'Vf'));
           hyst(devVoltage,:)=[-Vf, Vf/10];
           switchRef(devVoltage,:) = [i, 0];
          
       elseif isFET(i) && ~(any(isFET(switchInds)) && any(isDiode(switchInds)))
           % If a FET and there is no parallel diode
           %off-state voltage > 0
           %% I got rid of this because:
           % if you have an antiparallel diode, it was turning on the
           % MOSFET instead, then allowing negative curent (throught the
           % MOSFET).  I don't have a seperate indication for MOSFET
           % conduction vs body diode conduction at the moment.
           Cbnd(devVoltage, :, swvec(:,i)==0) = obj.topology.Cs(devVoltage, :, swvec(:,i)==0);
           Dbnd(devVoltage, :, swvec(:,i)==0) = obj.topology.Ds(devVoltage, :, swvec(:,i)==0);
           hyst(devVoltage,:)=[-obj.bodyDiodeVoltage, obj.bodyDiodeVoltageHysteresis];
           switchRef(devVoltage,:) = [i, 0];

           %off-state -(current) > 0
           %Note: this is currently useless -- current is == zero always
           %when switch is off.
%            Cbnd(devCurrent, :, swvec(:,i)==0) = -obj.topology.Cs(devCurrent, :, swvec(:,i)==0);
%            Dbnd(devCurrent, :, swvec(:,i)==0) = -obj.topology.Ds(devCurrent, :, swvec(:,i)==0);   
%            hyst(devCurrent,:)=[-1, .5];

       end
    end
    
    
    %Remove any zero rows
%     for i = size(Cbnd,1):-1:1
%         if ~(any(Cbnd(i,:,:)))
%             Cbnd(i,:,:) = [];
%             Dbnd(i,:,:) = [];
%             hyst(i,:) = [];
%             switchRef(i,:) = [];
%         end
%     end

    zeroRows = switchRef(:,1) == 0;
    Cbnd(zeroRows,:,:) = [];
    Dbnd(zeroRows,:,:) = [];
    hyst(zeroRows,:) = [];
    switchRef(zeroRows,:) = [];
    
    % eliminate redundant constraints
    combBnd = reshape(cat(2, Cbnd, Dbnd), size(Cbnd,1), (size(Cbnd,2) + size(Dbnd,2))*size(Cbnd,3) );
    [c, ia, ic] = unique(combBnd,'rows');
    
    Cbnd = Cbnd(ia,:,:);
    Dbnd = Dbnd(ia,:,:);
    hyst = hyst(ia,:);
    switchRef = switchRef(ia,:);

    
    obj.topology.Cbnd = Cbnd;
    obj.topology.Dbnd = Dbnd;
    obj.topology.bndHyst = hyst;
    obj.topology.switchRef = switchRef;

end

