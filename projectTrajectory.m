function [azimuth,elevation] = projectTrajectory(t,r)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
r = rotationNeptune(t, r);
x=r(1);
y=r(2);
z=r(3);
[azimuth,elevation,r] = cart2sph(x,y,z);
end

