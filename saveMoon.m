function path = saveMoon(t,t0_dist,path,keplerN)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rrN = arrayfun(@(t)TritonR(t,keplerN), (t-t0_dist)','UniformOutput',false);
rrN = cell2mat(rrN)';
dlmwrite(path,rrN,'precision',10)
end

