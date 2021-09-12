%Общая структура моделирования
load('D:\Копия с Pavillon\MATLAB\MVT\EJN\БУТ\gelio23.1.mat')
res = TEMP(:,20);

rN=mvt2icrf(res(17:19))*1e3;
VN=mvt2icrf(res(26:28))*1e3;

T12=res(16);
T23=res(29);
m3=res(20);
n3=res(25);

load('D:\Копия с Pavillon\MATLAB\MVT\EJN\БУТ\T_base.mat')
[year,month,day] = year2date(T_Neptune(m3,n3));
t_Neptune=juliandate([year,month,day]);
%t_Neptune = T12+T23;
%t_Neptune=juliandate(2049,6,10);

[rN_e, VN_e] = planetEphemeris(t_Neptune,'SolarSystem','Neptune','430');
rN_e=rN_e'*1e3;
VN_e=VN_e'*1e3;

%Интегрирование последнего куска траектории немного назад
tspan = [0, T23*365*24*3600/10];
%Размер сферы влияния
dist=1e9;
options = odeset('Events',@(t, y)stopByDist(t,y,rN_e,dist));
[t, y] = ode45(@(t,y) partialIntegrationSolar(t,y),tspan,[rN_e;-VN-VN_e],options);
r_border=y(end,1:3)';
V_border=-y(end,4:6)';
r0=r_border-rN_e;
V0=V_border;
t0_dist=t(end);
date = t_Neptune-t0_dist/24/3600;
%Начальные параметры по прибытию в систему Нептуна

ae = 149597870700;
%r0 = [-1 0 0]'*1 * 10 ^(6+3); %м/с
phi=pi/15;
%V0 = [cos(phi) sin(phi) 0]'*5.5301 *(10^3); %м/с


mug=6.809e15;
h0=V0'*V0/2-mug/norm(r0);
%r0T = [1 0 0]' * 354759e3; %м/с
%V0T = [0 1 0]' * 4.39e3; %м/с

%vInGel3 =[0.6650;5.6788;-0.5557];

y0= [r0; V0];
tspan=[0, t0_dist*0.75];
%Массив импульсов с направлениями и временами
psi=3*pi/4;
val_DV = 3800;
%dV=[cos(psi) sin(psi) 0]'*val_DV;
%tdV=[24*3600*1.6];
dV=[];
tdV=[];
%Интегрирование со знанием импульсов
[t, y] = complexIntegration(y0, dV, tdV, tspan);
%Визуализация

r0T = [1 0 0]' * 354759e3; %м/с
V0T = [0 1 0]' * 4.39e3; %м/с
y0T= [r0T; V0T];
tspanT=[0, 24*3600*5.87];
%Интегрирование Тритона для визуализации
[tT, yT] = complexIntegration(y0T, [], [], tspanT);
rT=yT(:,1:3);
rr=y(:,1:3);
figure(1);
%plot3(0, 0, 0, 'b--o')
R_Neptune = 24622000/ae;
[x,y,z] = sphere;
surf(R_Neptune*x, R_Neptune*y, R_Neptune*z);
set(gca,'FontSize',14)
hold on;
plot3(rr(:, 1)/ae, rr(:, 2)/ae, rr(:, 3)/ae, 'g', 'LineWidth', 1);

plot3(rT(:, 1)/ae, rT(:, 2)/ae, rT(:, 3)/ae, 'b', 'LineWidth', 2.5);
axis equal
%plot3(r0(1),r0(2),r0(3), 'go', 'LineWidth', 2.5);
title('Траектория КА')
xlabel('x, a.e.')
ylabel('y, a.e.')
zlabel('z, a.e.')
view(0,90)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
hold off;