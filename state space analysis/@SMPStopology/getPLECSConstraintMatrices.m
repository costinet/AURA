function [Cbnd, Dbnd, hyst, switchRef] = getPLECSConstraintMatrices(obj,circuitPath)
%getPLECSConstraintMatrics get Cbnd and Dbnd from PLECs model
%   Note that the PLECs model must have probes added to all switching
%   device currents and voltages (double-click -> assertions -> (+))

    S = strsplit(circuitPath,'/');
    modelFile = S{1};

    ssOrder = plecs('get', circuitPath, 'StateSpaceOrder');
    switchNames = ssOrder.Switches;
    outputs = ssOrder.Outputs;
    
    swvec = obj.swseq;
    
    Cbnd = zeros(size(obj.Cs));
    Dbnd = zeros(size(obj.Ds));
    hyst = zeros(size(obj.Ds,1), 2);
    switchRef = zeros(size(obj.Ds,1), 2);

    for i = 1:length(switchNames)
       switchSignals = strncmp(switchNames{i},outputs,length(switchNames{i}));
       assert(sum(switchSignals) >= 2, ['Switching device ' switchNames{i} ' does not have probes attached to both device current and device voltage']);

       currents = endsWith(outputs,'Am');
       voltages = endsWith(outputs,'Vm');

       devCurrent = find(switchSignals & currents);
       devVoltage = find(switchSignals & voltages);
       devType = plecs('get',[modelFile '/' switchNames{i}], 'Type');

       if strcmp(devType, 'Diode')
           %on-state current > 0
           Cbnd(devCurrent, :, swvec(:,i)==1) = obj.Cs(devCurrent, :, swvec(:,i)==1);
           Dbnd(devCurrent, :, swvec(:,i)==1) = obj.Ds(devCurrent, :, swvec(:,i)==1);
           hyst(devCurrent,:)=[0, 10e-6];
           switchRef(devCurrent,:) = [i, 1];

           %off-state -(voltage) > 0
           Cbnd(devVoltage, :, swvec(:,i)==0) = -obj.Cs(devVoltage, :, swvec(:,i)==0);
           Dbnd(devVoltage, :, swvec(:,i)==0) = -obj.Ds(devVoltage, :, swvec(:,i)==0);
           
           devStr = outputs{devVoltage};
           devStr = devStr(1:find(devStr==':',1)-1);
           Vf = evalin('base', plecs('get', [modelFile '/' devStr], 'Vf'));
           hyst(devVoltage,:)=[-Vf, Vf/10];
           switchRef(devVoltage,:) = [i, 0];
           

       elseif strcmp(devType, 'MosfetWithDiode') || strcmp(devType, 'Mosfet')
           %off-state voltage > 0
           Cbnd(devVoltage, :, swvec(:,i)==0) = obj.Cs(devVoltage, :, swvec(:,i)==0);
           Dbnd(devVoltage, :, swvec(:,i)==0) = obj.Ds(devVoltage, :, swvec(:,i)==0);
           hyst(devVoltage,:)=[-.5, .25];
           switchRef(devVoltage,:) = [i, 0];

           %off-state -(current) > 0
           %Note: this is currently useless -- current is == zero always
           %when switch is off.
%            Cbnd(devCurrent, :, swvec(:,i)==0) = -obj.Cs(devCurrent, :, swvec(:,i)==0);
%            Dbnd(devCurrent, :, swvec(:,i)==0) = -obj.Ds(devCurrent, :, swvec(:,i)==0);   
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

    
    obj.Cbnd = Cbnd;
    obj.Dbnd = Dbnd;
    obj.bndHyst = hyst;

end

