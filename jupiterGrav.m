load('D:\Копия с Pavillon\MATLAB\MVT_1\EJN\БУТ\gelio23.1.mat')
res = TEMP(:,20);

rN=mvt2icrf(res(17:19))*1e3;
VN=mvt2icrf(res(26:28))*1e3;
angle2vectors = @(u,v) atan2(norm(cross(u,v)),dot(u,v));
load('MVT/T_base.mat')
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

[year,month,day] = year2date(T_Jupiter(m2,n2));
t_Jupiter=juliandate(year,month,day);

[r_J,V_J] = planetEphemeris(t_Jupiter,'Sun','Jupiter','430');
r_J=r_J'*1e3;
V_J=V_J'*1e3;

Vhyp1=vInGel2-V_J;
Vhyp2=vOutGel2-V_J;

dirVper = Vhyp1+Vhyp2;
dirVper=dirVper/norm(dirVper);

dirRper = Vhyp1-Vhyp2;
dirRper=dirRper/norm(dirRper);

beta=angle2vectors(Vhyp1,Vhyp2);

mugJ=126686534*1e9; %м3/с2
a=mugJ/(Vhyp1'*Vhyp1);
RJ=69911000;
e=1/sin(beta/2);
rp=a*(e-1);
p=a*(e^2-1);

V_hyp_per = dirVper*sqrt(mugJ/p)*(e+1);
r_hyp_per = dirRper*rp;

dist=RJ*70;
Rrel=[0; 0; 0;];
options = odeset('Events',@(t, y)stopByDist(t,y,Rrel,dist));
options = odeset(options,'RelTol',1e-13);
tspan=linspace(0,40*24*3600,10000);
[t1, y] = ode45(@(t,y) partialIntegrationJupiter(t,y),tspan,[r_hyp_per;-V_hyp_per],options);
rr1=y(end:-1:1,1:3);
VV1=-y(end:-1:1,4:6);

[t2, y] = ode45(@(t,y) partialIntegrationJupiter(t,y),tspan,[r_hyp_per;V_hyp_per],options);
rr2=y(:,1:3);
VV2=y(:,4:6);
t2=t2+t1(end);
%в пределах 30 радиусов Юпитера
rr=[rr1;rr2];
VV=[VV1;VV2];

t=[t1;t2];
figure(1);
%plot3(0, 0, 0, 'b--o')

RA_North=2*pi*268.06/360; %https://en.wikipedia.org/wiki/Poles_of_astronomical_bodies
DE_North=2*pi*64.50/360;

[PoleNx,PoleNy,PoleNz]= sph2cart(RA_North,DE_North,4);
[x,y,z] = sphere(50);
surf(x, y, z);
set(gca,'FontSize',14)
hold on;

plot3(rr(:, 1)/RJ, rr(:, 2)/RJ, rr(:, 3)/RJ, 'g', 'LineWidth', 1);
plot3([PoleNx; -PoleNx], [PoleNy; -PoleNy], [PoleNz; -PoleNz], 'k', 'LineWidth', 1);
%Время пролёта в центре массива совпадает с датой 22.1.2033
tMoons = t-t1(end);
%Ио
keplerT = [203.319432880    421941.192 0.00425971  25.488489 293.700025 113.745069 358.148522];%22.1.2033
rrT = arrayfun(@(t)jupiterMoon(t,keplerT), tMoons','UniformOutput',false);
rrT = cell2mat(rrT)';
plot3(rrT(:, 1)/RJ, rrT(:, 2)/RJ, rrT(:, 3)/RJ, 'm', 'LineWidth', 2.5);
%Европа
keplerT = [101.373921510    671043.288 0.00965902  25.110625 323.144797 300.031607 358.573774];%22.1.2033
rrT = arrayfun(@(t)jupiterMoon(t,keplerT), tMoons','UniformOutput',false);
rrT = cell2mat(rrT)';
plot3(rrT(:, 1)/RJ, rrT(:, 2)/RJ, rrT(:, 3)/RJ, 'm', 'LineWidth', 2.5);
%Ганимед
keplerT = [50.318631422   1070425.532 0.00200786  25.625514  86.907240  15.029019 357.904160];%22.1.2033
rrT = arrayfun(@(t)jupiterMoon(t,keplerT), tMoons','UniformOutput',false);
rrT = cell2mat(rrT)';
plot3(rrT(:, 1)/RJ, rrT(:, 2)/RJ, rrT(:, 3)/RJ, 'm', 'LineWidth', 2.5);
%Каллисто
keplerT = [21.583172346   1882040.909 0.00739007  25.232828 243.194338  17.617524 358.198304];%22.1.2033
rrT = arrayfun(@(t)jupiterMoon(t,keplerT), tMoons','UniformOutput',false);
rrT = cell2mat(rrT)';
plot3(rrT(:, 1)/RJ, rrT(:, 2)/RJ, rrT(:, 3)/RJ, 'm', 'LineWidth', 2.5);

axis equal;
hold off;

rr_J_SO = arrayfun(@(t,x,y,z)rotationJupiter(t,[x,y,z]),...
    t(:),rr(:, 1),rr(:, 2),rr(:, 3),'UniformOutput',false);
rr_J_SO = cell2mat(rr_J_SO')';
dlmwrite('r-J-gelio.csv',rr_J_SO,'precision',10)
dlmwrite('t-J-gelio.csv',t,'precision',10)