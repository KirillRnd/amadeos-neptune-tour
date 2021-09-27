function res = partialIntegrationJupiter(t,y)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
res = zeros(1,6)';
r=y(1:3);
V=y(4:6);
res(1:3)=V;

mugJ=126686534*1e9; %м3/с2
gJupiter=-mugJ*r/norm(r)^3;
res(4:6)=gJupiter;
%mugNeptune=6.809e15;
end

