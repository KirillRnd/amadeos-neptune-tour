function res = partialIntegrationSolar(t,y)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
res = zeros(1,6)';
r=y(1:3);
V=y(4:6);
res(1:3)=V;

mugSun=132712.43994*(10^6)*(10^(3*3));
t_Neptune=juliandate(2050,3,31);
rNS= planetEphemeris(t_Neptune+t/24/3600,'Neptune','Sun')';
rSun=rNS*1e3-r;
gSun=-mugSun*rSun/norm(rSun)^3;

res(4:6)=gSun;
end

