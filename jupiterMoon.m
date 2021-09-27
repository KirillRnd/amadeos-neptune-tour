function RT = jupiterMoon(t,kepler)
%Положение Тритона в СК Нептуна в данный момент
%t должно отсчитываться от даты прилёта в секундах!
n=pi*kepler(1)/180/(24*3600);
a = kepler(2)*1e3;
e = kepler(3);
p=a*(1-e^2);
i=pi*kepler(4)/180;
M0=pi*kepler(5)/180;%На указанную дату
omega=pi*kepler(6)/180;
Omega=pi*kepler(7)/180;

M=M0+n*t;
E = EAnomaly(e, M);
nu=2*atan(tan(E/2)*sqrt((1+e)/(1-e)));%Только для Тритона

%Rflat=[a*cos(nu); a*sin(nu); 0];
%RT=rotmZYX*Rflat;
mugJ=126686534*1e9;
[RT, VT] = orb2rv(p,e,i,Omega,omega,nu,mugJ);
end
%Mean motion (deg/day), semi-major axis (km), eccentricity, inclination, mean anomaly,
%argument of pericenter, longitude of the ascending node are calculated independently.
%Год  Месяц День       n, deg/day      a, km       e          i, deg     M, deg   omega, deg Omega, deg  
%2049  6 10.000000    61.315381532    354532.843 0.00000000 111.880638 299.189584   0.000000 219.102919