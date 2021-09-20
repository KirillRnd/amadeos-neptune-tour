%Общая структура моделирования 8.5.2049
load('D:\Копия с Pavillon\MATLAB\MVT_1\EJN\БУТ\gelio23.1.mat')
res = TEMP(:,20);

rN=mvt2icrf(res(17:19))*1e3;
VN=mvt2icrf(res(26:28))*1e3;

T12=res(16);
T23=res(29);
m3=res(20);
n3=res(25);

load('D:\Копия с Pavillon\MATLAB\MVT_1\EJN\БУТ\T_base.mat')
[year,month,day] = year2date(T_Neptune(m3,n3));
% year=2049;
% month=5;
% day=8;
t_Neptune=juliandate(year,month,day);
%T = 2.631286788824651e+04; %период орбиты на расттоянии 2RN
T = 4.423802058240195e+05; %период орбиты на расттоянии a=(25RN+1RN+1e7)/2
tdV_1 = 1.539377738365407e+07;
SunLine = planetEphemeris(t_Neptune,'Neptune','Sun','430');
SunLine=SunLine/norm(SunLine);
%t_Neptune = T12+T23;
%t_Neptune=juliandate(2049,6,10);

[rN_e, VN_e] = planetEphemeris(t_Neptune,'Sun','Neptune','430');
rN_e=rN_e'*1e3;
VN_e=VN_e'*1e3;

%Интегрирование последнего куска траектории немного назад
tspan = [0, T23*365*24*3600/10];
%Размер сферы влияния
dist=87e9;

Rrel=[0; 0; 0;];
options = odeset('Events',@(t, y)stopByDist(t,y,Rrel,dist));
options = odeset(options,'RelTol',1e-13);
Vrel=VN-VN_e;
[t, y] = ode45(@(t,y) partialIntegrationSolar(t,y),tspan,[Rrel;-Vrel],options);
r_border=y(end,1:3)';
V_border=-y(end,4:6)';
r0=r_border;
V0=V_border;
t0_dist=t(end);
date = t_Neptune-t0_dist;
%Начальные параметры по прибытию в систему Нептуна

ae = 149597870700;
%r0 = [-1 0 0]'*1 * 10 ^(6+3); %м/с

%V0 = [cos(phi) sin(phi) 0]'*5.5301 *(10^3); %м/с


mug=6.809e15;
h0=V0'*V0/2-mug/norm(r0);
%r0T = [1 0 0]' * 354759e3; %м/с
%V0T = [0 1 0]' * 4.39e3; %м/с

%vInGel3 =[0.6650;5.6788;-0.5557];

%Определяем кеплеровы параметры целевой орбиты.

% h=cross(SunLine,r0);
% h=h/norm(h);
% a = 24622000*40;
% e = 0.5;
% rPer=SunLine*a*(1-e);
% VPerDir=cross(SunLine,h);
% VPerVal=sqrt(2*(mug/a+mug/norm(rPer)));
% VPer=VPerDir*VPerVal;
% [a_1,eMag_1,i_1,O_1,o_1,nu_1,truLon_1,argLat_1,lonPer_1,p_1] = rv2orb(rPer',VPer',mug);
% 
% angle2vectors = @(u,v) atan2(norm(cross(u,v)),dot(u,v));
% %SunLineProj = [SunLine(1) SunLine(2) 0];
% %i=angle2vectors(SunLine,SunLineProj);
% i=2.7831;
% %Omega=angle2vectors([1 0 0],r0);
% Omega=3.0920;
% omega=5.3698;
%omega=2*pi-angle2vectors(SunLine,r0);

%keplerTarget = [61.315381532;    a/1000; e; 180*i/pi; 0;   180*omega/pi; -180*Omega/pi];%2049  6 10.000000  геоэкваториальные


y0= [r0; V0];

%Массив импульсов с направлениями и временами

val_DV1 = 8.9;%затраты хар. скорости м/с
dir_dV1 = -cross(V0,cross(SunLine,V0));
dir_dV1=dir_dV1/norm(dir_dV1);

val_DV2 = 1249;%1126;%затраты хар. скорости м/с
dir_dV2=-1.0e+04 *[1.778878827943030  -1.155881837047985  -0.487994013230544];
dir_dV2=dir_dV2/norm(dir_dV2);

val_DV3 = 0;
phi=7.4*pi/6;
psi=0*pi/6;
dir_dV3 =[cos(phi)*cos(psi)  sin(phi)*cos(psi)  sin(psi)];
dir_dV3=dir_dV3/norm(dir_dV3);

dV=[dir_dV1*val_DV1;[0, 0, 0];dir_dV2*val_DV2;];
tdV=[t0_dist*0.1,t0_dist*0.98,tdV_1];
lsp=[100;100;10000;10000];

% dV=[dir_dV1*val_DV1;[0, 0, 0];dir_dV2*val_DV2;dir_dV3*val_DV3;];
% tdV=[t0_dist*0.1,t0_dist*0.9,tdV_1,tdV_1+0.964*T];
% lsp=[100;100;10000;10000;10000];

%tspan=[0, tdV_1+100*T];
tspan=[0,tdV_1+T];
% dV=[];
% tdV=[];
% lsp=[100];
%Интегрирование со знанием импульсов
[t, y] = complexIntegration(y0, dV, tdV, tspan, lsp);
rr=y(:,1:3);
VV=y(:,4:6);
%Визуализация
%http://www.sai.msu.ru/neb/nss/html/multisat/nssreq8hr.htm
keplerT = [61.315381532;    354532.843; 0.00000000; 111.874553;  77.689068;   0.000000; 219.084078];%8.5.2049
%keplerT = [61.315381532;    354532.843; 0.00000000; 111.869215; 101.218918;   0.000000; 219.067518];%2049  4 9  геоэкваториальные
%keplerT = [61.315381532;    354532.843; 0.00000000; 111.880638; 299.189584;   0.000000; 219.102919];%2049  6 10.000000  геоэкваториальные
%kepler = [61.315381532;    354532.843; 0.00000000; 128.929671; 318.004206;   0.000000; 228.795867];%2049  6 10.000000  геоэклиптические


% r0T = [-150957.08108; -306953.91765; -93813.85880]*1e3; %MULTISAT МГУ Geo-planetocentric 2049 6 10
% V0T = [-337961.32962; 132768.17325; 109407.55504]*1e3/(24*3600); %MULTISAT МГУ Geo-planetocentric 2049 6 10
% 
% r0T = [-61454.95752;-296564.43426;-184614.76788]*1e3; %MULTISAT МГУ Geo-ecliptic 2049 6 10
% V0T = [-300389.30109;-74237.70314;219249.38420]*1e3/(24*3600); %MULTISAT МГУ Geo-ecliptic 2049 6 10

r0T = [-61454.95752;-198657.01137;-287347.29513]*1e3; %MULTISAT МГУ Geo-equatorial 2049 6 10
V0T = [-300389.30109;-155324.15744;171627.31470]*1e3/(24*3600); %MULTISAT МГУ Geo-equatorial 2049 6 10

   
   
%r0T = [1 0 0]' * 354759e3; %м/с
%V0T = [0 1 0]' * 4.39e3; %м/с
%y0T= [r0T; V0T];
%tspanT=[0, 24*3600*5.9];
%Интегрирование Тритона для визуализации
%[tT, yT] = complexIntegration(y0T, [], [], tspanT);
%rT=yT(:,1:3);

RA_North=pi*299.36/180; %https://en.wikipedia.org/wiki/Poles_of_astronomical_bodies
DE_Nort=pi*43.46/180;



figure(1);
%plot3(0, 0, 0, 'b--o')
RN=24622000;
R_Neptune = 24622000/ae;
[PoleNx,PoleNy,PoleNz]= sph2cart(RA_North,DE_Nort,R_Neptune*4);
[x,y,z] = sphere(50);
surf(R_Neptune*x, R_Neptune*y, R_Neptune*z);
set(gca,'FontSize',14)
hold on;


plot3([R_Neptune*20*SunLine(1) R_Neptune*25*SunLine(1)], ...
    [R_Neptune*20*SunLine(2) R_Neptune*25*SunLine(2)],...
    [R_Neptune*20*SunLine(3) R_Neptune*25*SunLine(3)], 'y-x', 'LineWidth', 3);

plot3([-R_Neptune*50*SunLine(1) -R_Neptune*80*SunLine(1)], ...
    [-R_Neptune*50*SunLine(2) -R_Neptune*80*SunLine(2)],...
    [-R_Neptune*50*SunLine(3) -R_Neptune*80*SunLine(3)], 'y-x', 'LineWidth', 3);


plot3(rr(:, 1)/ae, rr(:, 2)/ae, rr(:, 3)/ae, 'g', 'LineWidth', 1);
plot3([PoleNx; -PoleNx], [PoleNy; -PoleNy], [PoleNz; -PoleNz], 'k', 'LineWidth', 1);

%plot3(rr(100, 1)/ae, rr(100, 2)/ae, rr(100, 3)/ae, 'rO', 'LineWidth', 1.5);
ind1=lsp(1)+lsp(2)+1;
ind2=lsp(1)+lsp(2)+lsp(3);
minimum = min(vecnorm(rr(ind1:ind2,:),2,2));
ind_per = find(minimum == vecnorm(rr(ind1:ind2,:),2,2))+lsp(1)+lsp(2);
plot3(rr(ind_per, 1)/ae, rr(ind_per, 2)/ae, rr(ind_per, 3)/ae, 'rO', 'LineWidth', 1.5);
%ind_3_man=lsp(1)+lsp(2)+lsp(3)+lsp(4);
%plot3(rr(ind_3_man, 1)/ae, rr(ind_3_man, 2)/ae, rr(ind_3_man, 3)/ae, 'rO', 'LineWidth', 1.5);

%plot3([rr(1, 1)/ae;rr(1, 1)/ae+V0(1)*1e-7], [rr(1, 2)/ae;rr(1, 2)/ae+V0(2)*1e-7], [rr(1, 3)/ae;rr(1, 3)/ae+V0(3)*1e-7], 'r', 'LineWidth', 1);

%plot3([rr(1, 1)/ae;rr(1, 1)/ae+dV(1)*1e-7], [rr(1, 2)/ae;rr(1, 2)/ae+dV(2)*1e-7], [rr(1, 3)/ae;rr(1, 3)/ae+dV(3)*1e-7], 'r', 'LineWidth', 1);



rrT = arrayfun(@(t)TritonR(t,keplerT), linspace(0,5.90*24*3600,50),'UniformOutput',false);
rrT = cell2mat(rrT)';

%plot3(rT(:, 1)/ae, rT(:, 2)/ae, rT(:, 3)/ae, 'b', 'LineWidth', 2.5);

plot3(rrT(:, 1)/ae, rrT(:, 2)/ae, rrT(:, 3)/ae, 'b', 'LineWidth', 2.5);

%rT_1 = TritonR(tdV(4),keplerT);
%plot3(rT_1(1)/ae, rT_1(2)/ae, rT_1(3)/ae, 'bO', 'LineWidth', 2.5);
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
xlim([-50*R_Neptune 65*R_Neptune])
ylim([-50*R_Neptune 65*R_Neptune])
zlim([-40*R_Neptune 40*R_Neptune])
box off;
hold off;