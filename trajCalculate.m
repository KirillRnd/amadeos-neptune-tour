load('D:\Копия с Pavillon\MATLAB\MVT_1\EJN\БУТ\gelio23.1.mat')
res = TEMP(:,20);

rN=mvt2icrf(res(17:19))*1e3;
VN=mvt2icrf(res(26:28))*1e3;


m1=       res(2);      % Номер строки в T_base для узла первой планеты
n1=       res(3);      % Номер столбца в T_base для узла первой планеты
r1=       mvt2icrf(res(4:6))*1e3;    % Радиус-вектор первой планеты маршрута
vinf1=    res(7)*1e3;      % Гиперболический избыток скорости при старте
r2=       mvt2icrf(res(8:10))*1e3;   % Радиус-вектор второй планеты маршрута
m2=        res(11);     % Номер строки в T_base для узла второй планеты
n2=        res(12);    % Номер столбца в T_base для узла второй планеты
vInGel2=   mvt2icrf(res(13:15))*1e3;  % Гелиоцентрическая скорость встречи со второй планетой маршрута
T12=      res(16);     % Время перелета между первой и второй планетами маршрута

% Далее периодичность в структуре

r3=        mvt2icrf(res(17:19))*1e3;  % Радиус-вектор третьей планеты маршрута
m3=        res(20);     % Номер строки в T_base для узла третьей планеты
dv2=       res(21);     % Затраты характеристической скорости у второй планеты маршрута
vOutGel2=  mvt2icrf(res(22:24))*1e3;  % Гелиоцентрическая скорость отбытия от второй планеты маршрута
n3=        res(25);     % Номер столбца в T_base для узла третьей планеты
vInGel3=   mvt2icrf(res(26:28))*1e3;  % Гелиоцентрическая скорость встречи c третьей планетой маршрута
T23=       res(29);


load('D:\Копия с Pavillon\MATLAB\MVT_1\EJN\БУТ\T_base.mat')
[year,month,day] = year2date(T_Earth(m1,n1));
t_Earth=juliandate(year,month,day);

[year,month,day] = year2date(T_Jupiter(m2,n2));
t_Jupiter=juliandate(year,month,day);

[year,month,day] = year2date(T_Neptune(m3,n3));
t_Neptune=juliandate(year,month,day);

[r_e,V_e] = planetEphemeris(t_Earth,'Sun','Earth','430');
r_e=r_e'*1e3;
V_e=V_e'*1e3;
v1 = (norm(V_e)+vinf1)*V_e/norm(V_e);
tspan12 = linspace(t_Earth,t_Jupiter,100)*24*3600;
%tspan12 = linspace(0,T12*365*24*3600,10000);
[t12, y12] = ode45(@(t,y) partialIntegrationSolar(t,y),tspan12,[r2;-vInGel2]);
rr12=y12(end:-1:1,1:3);
VV12=y12(end:-1:1,4:6);

tspan23 = linspace(t_Jupiter,t_Neptune,100)*24*3600;
[t23, y23] = ode45(@(t,y) partialIntegrationSolar(t,y),tspan23,[r2;vOutGel2]);
rr23=y23(:,1:3);
VV23=y23(:,4:6);
ae = 149597870700;
figure(3)

rr=[rr12;rr23];
VV=[VV12;VV23];
t=[t12;t23];

plot3(rr(:,1)/ae,rr(:,2)/ae,rr(:,3)/ae,'k');
% plot3(rr12(:,1)/ae,rr12(:,2)/ae,rr12(:,3)/ae);
hold on;
plot3(rr12(end,1)/ae,rr12(end,2)/ae,rr12(end,3)/ae,'+k');
rrE = arrayfun(@(t)planetEphemeris(t,'Sun','Earth','430'), t/24/3600,'UniformOutput',false);
rrE = 1e3*cell2mat(rrE);
plot3(rrE(:,1)/ae,rrE(:,2)/ae,rrE(:,3)/ae,'--b');

rrP = arrayfun(@(t)planetEphemeris(t,'Sun','Mars','430'), t/24/3600,'UniformOutput',false);
rrP = 1e3*cell2mat(rrP);
plot3(rrP(:,1)/ae,rrP(:,2)/ae,rrP(:,3)/ae,'--b');

rrJ = arrayfun(@(t)planetEphemeris(t,'Sun','Jupiter','430'), t/24/3600,'UniformOutput',false);
rrJ = 1e3*cell2mat(rrJ);
plot3(rrJ(:,1)/ae,rrJ(:,2)/ae,rrJ(:,3)/ae,'--b');

rrP = arrayfun(@(t)planetEphemeris(t,'Sun','Saturn','430'), t/24/3600,'UniformOutput',false);
rrP = 1e3*cell2mat(rrP);
plot3(rrP(:,1)/ae,rrP(:,2)/ae,rrP(:,3)/ae,'--b');

rrP = arrayfun(@(t)planetEphemeris(t,'Sun','Uranus','430'), t/24/3600,'UniformOutput',false);
rrP = 1e3*cell2mat(rrP);
plot3(rrP(:,1)/ae,rrP(:,2)/ae,rrP(:,3)/ae,'--b');

rrP = arrayfun(@(t)planetEphemeris(t,'Sun','Neptune','430'), t/24/3600,'UniformOutput',false);
rrP = 1e3*cell2mat(rrP);
plot3(rrP(:,1)/ae,rrP(:,2)/ae,rrP(:,3)/ae,'--b');

axis equal;
hold off;
box off;
xlabel('x, a.e.')
ylabel('y, a.e.')
% plot3(rr23(:,1)/ae,rr23(:,2)/ae,rr23(:,3)/ae);
% hold off;

dlmwrite('r-gelio.csv',rr,'precision',10)
dlmwrite('V-gelio.csv',VV,'precision',10)
dlmwrite('t-gelio.csv',t,'precision',10)