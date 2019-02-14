function [Xs] = CorrectXs(obj,keep_SS,Xs)
%CORRECTXS Takes the given Xs and corrects dependent states to reflect
%actual values that contain KVL and KCL rules
%   Must come after StateVarIndex to get the index of dependent states
%   outputs
%
%   Can input Xs put need to have associated Ds, Cs, and u in class 
%
%     %%%%%%   %      %  %%%%%%%    %%%%%%
%    %      %  %      %  %      %  %      %
%    %      %  %      %  %      %  %      %
%    %%%%%%%%  %      %  %%%%%%%   %%%%%%%%
%    %      %  %      %  %%        %      %
%    %      %  %      %  % %       %      %
%    %      %  %      %  %  %      %      %
%    %      %  %      %  %   %     %      %
%    %      %   %    %   %    %    %      %
%    %      %    %%%%    %     %   %      %

if nargin==1
    As = obj.As;
    Bs = obj.Bs;
    Cs = obj.Cs;
    Ds = obj.Ds;
    Xs = obj.Xs;
    u = obj.u;
    keep_SS = false;
    
elseif nargin == 2
    As = obj.As;
    Bs = obj.Bs;
    Cs = obj.Cs;
    Ds = obj.Ds;
    Xs = obj.Xs;
    u = obj.u;
    
elseif nargin == 3
    As = obj.As;
    Bs = obj.Bs;
    Cs = obj.Cs;
    Ds = obj.Ds;
    u = obj.u;
end

OutputNames = obj.Converter.Topology.Parser.OutputNames;
% DependentNames = obj.Converter.Topology.Parser.DependentNames;
StateNumbers = obj.Converter.Topology.Parser.StateNumbers;

%% Reconstruction of Dependent variables
% Dont need anymore due to fix of SS_Soln.m
% Dependent variables are calculated with independent variables

%Xs(end+1,:) = zeros(size(DependentNames,1),size(Xs,2));


% For now diode is 8th row in C and D will char match dependent variables
% to find where they are in output and calculate
start = size(OutputNames,1);
for i = 2:1:size(Xs,2)
    Xs(start+1:end,i) = Cs(StateNumbers(start+1:end),:,i-1)*Xs(:,i)+Ds(StateNumbers(start+1:end),:,i-1)*u;
    if ~keep_SS
        Xs(:,1) = Xs(:,end);
    end
end
if nargin==1
    obj.Xs_circuit = Xs;
    obj.Xs = Xs;
end

end % That's all Folks

