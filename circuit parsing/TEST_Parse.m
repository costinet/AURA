%% TEST Parse Class

% This script tests the circuit parsing code given the filename,
% testcase and measurements input

clear

% Set path
testDir = [erase(mfilename('fullpath'), mfilename) 'NetLists'];
addpath(testDir);
testDir = [erase(mfilename('fullpath'), mfilename) 'TestCases'];
addpath(testDir);

%% User Input

% Select .net file
filename = 'Buck.net';
% Current options for filename:
% Boost.net
% Buck.net
% Buck2.net
% Dickson.net
% Flyback.net
% Forward.net

% Select test case file
testcase = 'TEST_ABCD_Buck';
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


%% Run functions
parse = NetListParse();
parse.initialize(filename,Voltage,Current);
parse.ABCD();
testfun = str2func(testcase);
testfun(parse);

    
    
    