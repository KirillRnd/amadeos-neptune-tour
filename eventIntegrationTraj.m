function [value, isterminal, direction] = stopByDist(s, y, r0,dist)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r=y(1:3);
value = dist - norm(r-r0);
isterminal = 1;
direction = 0;
end

