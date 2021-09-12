function [value, isterminal, direction] = stopByDist(t, y, r0, dist)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r=y(1:3);
%rN_e = planetEphemeris(t0-t/24/3600,'SolarSystem','Neptune','430');
%rN_e=rN_e'*1e3;
value = dist-norm(r0-r);

isterminal = 1;
direction = 0;
end