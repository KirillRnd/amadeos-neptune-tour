function res = partialIntegration(t,y,body)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
res = zeros(1,6)';
r=y(1:3);
V=y(4:6);
res(1:3)=V;
[gx_zonal, gy_zonal, gz_zonal] = gravityzonal(r', 'Neptune');
gNeptune=[gx_zonal; gy_zonal; gz_zonal];



%Влияние солнца
mugSun=132712.43994*(10^6)*(10^(3*3));
t_Neptune=juliandate(2050,3,31);
rNS= planetEphemeris(t_Neptune+t/24/3600,'Neptune','Sun')';
rSun=rNS*1e3-r;
gSun=-mugSun*rSun/norm(rSun)^3;

%Влияние Тритона
keplerT = [61.315381532    354532.843 0.00000000 111.935349 308.921483   0.000000 219.270622];%31.3.2050
rTr = TritonR(t, keplerT)-r;
mugTr=1427.6*1e9;%https://ssd.jpl.nasa.gov/?sat_phys_par
gTr=-mugTr*rTr/norm(rTr)^3;

%res(4:6)=gNeptune+gTr;%+gSun;
res(4:6)=gNeptune+gTr+gSun;
%mugNeptune=6.809e15;
%res(4:6)=-mugNeptune*r/norm(r)^3;
end

