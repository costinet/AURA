function [Combinations] = BDCheck(obj,Combinations,state,switches)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%     _   _   _  ____    _
%    / \ | | | |/ _  |  / \
%   / _ \| | | | (_| | / _ \
%  / ___ | |_| |> _  |/ ___ \
% /_/   \_\___//_/ |_/_/   \_\


BD_state = obj.BD_state;
BD_OFF_state = obj.BD_OFF_state;


for i = 1:1:length(Combinations)
    for j = 1:1:size(state,2)
        if sum(repmat(BD_state(Combinations(i),j),1,length(Combinations))==Combinations)==0
            Combinations(end+1) = BD_state(Combinations(i),j);
        end
        if sum(repmat(BD_OFF_state(Combinations(i),j),1,length(Combinations))==Combinations)==0
            Combinations(end+1) = BD_OFF_state(Combinations(i),j);
        end
    end
    
end

end % That's all Folks

