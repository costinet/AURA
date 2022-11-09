function [newts,newswinds] = format1DtimingVector(obj,ts,swinds)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    warning("format1DtimingVector still isn't right -- it may count BD turnoff as a new controlled switching action")

    top = obj.converter.topology;
    switchedDevs = diff(top.swseq(swinds,:),1);
    diodes = ~cellfun(@isempty,regexp(top.switchLabels, 'D.*'));
    activeSwitched = switchedDevs(:,~diodes);
    newInt = sum(abs(activeSwitched),2)~=0;

    newts = zeros(length(newInt)+1, sum(newInt)+1);
    newswinds = zeros(length(newInt)+1, sum(newInt)+1);

    newts(1) = ts(1);
    newswinds(1) = swinds(1);
    col = 1;
    row = 2;
    for i = 1:length(newInt)
        if newInt(i) == 1
            col = col + 1;
            row = 1;
        end
        newts(row,col) = ts(i+1);
        newswinds(row,col) = swinds(i+1);
        row = row+1;
    end

    newts(all(newswinds == 0,2),:)=[];
    newswinds(all(newswinds == 0,2),:)=[];

end