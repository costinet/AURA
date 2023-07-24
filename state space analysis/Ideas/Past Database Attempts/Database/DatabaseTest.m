%% Setup

clear all 
close all
clc
addpath(genpath(pwd))


%% Create Transistor Object

% GUIs for entering data will open automatically in the order:
%   1. Datasheet table parameters 
%       Standard static characteristics, no conditions, no experimental data.
%   2. Additional table parameters 
%       Other static characteristics, may have conditions, may be experimental.
%   3. Graph parameters 
%       Dynamic characteristics, may have conditions, may be experimental.
%       Use WebPlotDigitizer to get data points - http://web.eecs.utk.edu/~dcostine/personal/PowerDeviceLib/DigiTest/index.html
%       Use the dropdown at the top to change parameters or add a new parameter.
%       Use the 'Previous Page' and 'Next/Add Page' buttons to browse/add multiple copies of data for the parameter given in the dropdown.
%
% Add data and save, then close the GUI window to open the next window. 
% Observe unit prefixes carefully!
% Change the device number to enter a new device.
deviceNumber = 'EPC2014C';
t = Transistor(deviceNumber);


%% View Transistor Datasheet

% Datasheet contains all types of information in one page and is non-editable.
t.datasheet;

%% Get Table Parameter
disp(newline)
devices = Transistor.listDevices;
for i = 1:numel(devices)
    device = devices{i};
    paramValue = Transistor.getTableParam(device,'id','max');
    disp([device ' Id max. = ' num2str(paramValue) ' (A)']);
end

%% Get FOM Using Get Table Parameter
disp(newline)
devices = Transistor.listDevices;
for i = 1:numel(devices)
    device = devices{i};
    rds_typ = Transistor.getTableParam(device,'rds_on','typ');
    qg_typ = Transistor.getTableParam(device,'Qg','typ');
    disp(['FOM for ' device ' is ' num2str(rds_typ*qg_typ)])
end

%% Get Graph Data
devices = Transistor.listDevices;
figure
set(gcf, 'Position',  [500, 300, 1200, 400])

subplot(1,2,1)
for i = 1:numel(devices)
    device = devices{i};
    graphData = Transistor.getGraphData(device,'Coss','Vds');
    plot1 = graphData{1};
    hold on
    xlim([0 inf]);
    set(gca, 'YScale', 'log');
    plot(plot1(:,1),1e12.*plot1(:,2))
end
legend(devices)
ylabel('Coss (pF)')
xlabel('Vds (V)')
title('Output Capacitance Curves for Various Devices')

subplot(1,2,2)
for i = 1:numel(devices)
    device = devices{i};
    graphData = Transistor.getGraphData(device,'Vgs','Qg');
    plot1 = graphData{1};
    hold on
    xlim([0 inf]);
    plot(1e9.*plot1(:,1),plot1(:,2))
end
legend(devices)
ylabel('Vgs (V)')
xlabel('Qg (nC)')
title('Gate Charge Curves for Various Devices')

%% Get Graph Point
disp(newline)
devices = Transistor.listDevices;
for i = 1:numel(devices)
    device = devices{i};
    [yValue, ConditionsStruct] = Transistor.getGraphPoint(device,'Id','Vds',.7);
    strings = Transistor.ConditionsStructToStr(ConditionsStruct);
    for j = 1:numel(yValue)
        if isnan(yValue(j))
            str = ['For ' device ', there is no data on Id when Vds = .5(V)'];
        else
            str = ['For ' device ', Id = ' num2str(yValue(j)) '(A) when Vds = .5(V)'];
        end
        if ~isempty(strings{j})
            str = [str ' and when ' strings{j}];
        end
        disp(str)
    end
    disp(' ')
end

%% Compare Graph Parameter Curves with Conditions for a Given Device
devices = Transistor.listDevices;
device = devices{end};
figure;
ax = gca;
Transistor.showGraphData(device,'Id','Vds', ax);

%% Find Devices with Vds >= 60V and Id < 50A
disp(newline)
validDevices1 = Transistor.filter('Vds','max','>=',60);
validDevices2 = Transistor.filter('Rds_on','typ','<',3e-3);
validDevices3 = intersect(validDevices1,validDevices2);
disp(['Vds (max) >= 60V: ' newline strjoin(validDevices1,', ') newline])
disp(['Rds_on (typ) < 3mΩ: ' newline strjoin(validDevices2,', ') newline])
disp(['Both: ' newline strjoin(validDevices3,', ') newline])


%% Filter Devices & Get Qoss @ Given Vds
disp(newline)
devices = Transistor.filter('Vds', 'max', '>=', 50);
for i = 1:numel(devices)
    device = devices{i};
    Qoss = Transistor.getChargeFromCoss(device,50);
    for j = 1:numel(Qoss)
        disp([device ' Qoss=' num2str(Qoss{j}*1e9) ' (nC)']);
    end
end

%% Other Methods
disp(newline)
methods(t)
