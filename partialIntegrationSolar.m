function res = partialIntegrationSolar(t,y)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
res = zeros(1,6)';
r=y(1:3);
V=y(4:6);
res(1:3)=V;
mugSun=132712.43994*(10^6)*(10^(3*3));
res(4:6)=-mugSun*r/norm(r)^3;
end

