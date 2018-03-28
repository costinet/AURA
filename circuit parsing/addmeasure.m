function [NL,NLnets,MeasureFlags] = addmeasure(NL,NLnets)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% List of Voltages wanted to measure

% Example:
% Measure Voltage from Node 1 to Node 2 and Node 4 to Node 3 would be
% written as:
% 
% Voltage = [1 2
%            4 3];

Voltage = [2 1];



% List of Currentes to Measusre
% Currents listed as going from Node 1 to Node 2 and Node 4 to Node 3
% should be listed as:
% 
% Current = [1 2
%     4 3];

Current = [4 1];

D = size(Voltage);

Flag = zeros(size(NL,1),1);

for i = 1:1:D(1)

Match = NL(:,2:3) == Voltage(i,:);
Matches = Match(:,1).*Match(:,2);

Flag(Matches==1,:) = true;



end

NL(Matches==1,:) = [];

