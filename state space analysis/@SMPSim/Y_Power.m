function [Y_Power,Cds_Loss_Voltage] = Y_Power(obj)
%Y_Power Calculateds the power consumed for all of the elements in the circuit.
%   The results are palced as a class in the SMPSim class matrix

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\

%% Calculate <V*I> for all outputs

As = obj.As;
Bs = obj.Bs;
Cs = obj.Cs;
Ds = obj.Ds;
ts = obj.ts;
u = obj.u;

ns = size(As(:,:,1),1);
pi = length(u);
mo = size(Cs(:,:,1),1);


Xov = obj.Xs(:,1);
XTsv = [Xov;u;zeros(ns,1)];

for i = 1:1:size(As,3)
    Aa(:,:,i) = [As(:,:,i), Bs(:,:,i);
        zeros(pi,ns+pi)];
    Aav(:,:,i) = [Aa(:,:,i), zeros(ns+pi,ns);
        eye(ns), zeros(ns,ns+pi)];
    if i == 1
        XTsv(:,i) = expm(Aav(:,:,i)*ts(i))*XTsv;
    else
        XTsv(:,i) = expm(Aav(:,:,i)*ts(i))*XTsv(:,i-1);
    end
end


X(:,:) = XTsv(ns+pi+1:end,:);

for i = 1:1:size(As,3)
    if i == 1
        Y(:,i) = Cs(:,:,i)*X(:,1)+Ds(:,:,i)*u*ts(i);
    else
        Y(:,i) = Cs(:,:,i)*(X(:,i)-X(:,i-1))+Ds(:,:,i)*u*ts(i);
    end
    
end

Ytot = sum(Y,2)/sum(ts);

%Y = Cs(:,:,end)*X+Ds(:,:,end)*u;

split = size(Y,1)/2;

Y_Power = Ytot(1:split).*Ytot(split+1:end);


%% Find the Coss values of all device FETs

% this is the voltage difference from time j to j+1

ONorOFF = obj.Converter.Topology.Parser.ONorOFF;
Cds_Loss_Voltage = zeros(size(ONorOFF,1),1);
for i = 1:1:size(ONorOFF,1)
    for j = 1:1:size(ONorOFF,2)
        if j == size(ONorOFF,2)
            if ONorOFF(i,1)==2 && ONorOFF(i,j)==-1
                Cds_Loss_Voltage(i,j) = obj.Xs(i,j+1);
            end
        else
            if ONorOFF(i,j+1)==2 && ONorOFF(i,j)==-1
                Cds_Loss_Voltage(i,j) = obj.Xs(i,j+1);
            end
        end
    end 
end


%% Finding the Power loss due to Resistors








end

