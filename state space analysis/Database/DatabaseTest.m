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
t = Transistor('EPC2024');


%% View Transistor Datasheet

% Datasheet contains all types of information in one page and is non-editable.
t.datasheet;


%% Other Methods

% View the other public methods. There are push/pull methods for merging databases that will be used later.
methods(t)

% getChargeFromCoss finds the output capicitance charge at a given Vds.
charge = t.getChargeFromCoss(20);
disp(['Charge= ' num2str(charge*1e9) ' (nC)']);

