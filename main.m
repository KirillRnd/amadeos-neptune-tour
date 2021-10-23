%Общая структура моделирования 31.3.2050
load('MVT/gelio23.1.mat')
res = TEMP(:,20);

rN=mvt2icrf(res(17:19))*1e3;
VN=mvt2icrf(res(26:28))*1e3;

T12=res(16);
T23=res(29);
m3=res(20);
n3=res(25);

load('MVT/T_base.mat')
[year,month,day] = year2date(T_Neptune(m3,n3));
% year=2049;
% month=5;
% day=8;
RT=1353400;%Радиус Тритона
t_Neptune=juliandate(year,month,day);
%T = 2.631286788824651e+04; %период орбиты на расттоянии 2RN
T0 =  1.640023631767377e+06;
T = 4.292471197610088e+05; %период орбиты на расттоянии a=(25RN+minimum)/2
T2 = 2.722295586988940e+06;
TTr = 2.327317237000118e+06;
Tend = 1.763575681396873e+06;
TTr4_1=5.87*4*24*3600;
TTr3_1=5.87*3*24*3600;
TTr2_1=5.87*2*24*3600;
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
[t, y] = ode45(@(t,y) partialIntegrationSolarNeptune(t,y),tspan,[Rrel;-Vrel],options);
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
angle2vectors = @(u,v) atan2(norm(cross(u,v)),dot(u,v));
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

val_DV1 = 8.59;%затраты хар. скорости м/с
dir_dV1 = -cross(V0,cross(SunLine,V0));
dir_dV1=dir_dV1/norm(dir_dV1);

%Первый тормозной импульс
val_DV2 = -950;%1126;%затраты хар. скорости м/с
dir_dV2=1.0e+04 *[1.781998893321638  -1.142666122750968  -0.484724116406552];
dir_dV2=dir_dV2/norm(dir_dV2);

tdV_2 = 1.539525149705983e+07;

%Второй тормозной импульс
%val_DV3 = -3.084123550418008e+02;
val_DV3 = -3.084123550418008e+02;
dir_dV3 = 1.0e+04 *[1.665011548038263  -1.141106172321206  -0.486895771533629];
dir_dV3=dir_dV3/norm(dir_dV3);

tdV_3 = tdV_2+T0;

%Поднятие до хвоста магнитопаузы первый импульс
val_DV4 = 1568;
dir_dV4 = 1.0e+02 *[-8.213716499245859   5.963945086903120   2.730439197810426];
dir_dV4=dir_dV4/norm(dir_dV4);
t_10T = 2.132602916624223e+07;
tdV_4 = t_10T+0.5*T;
%Поднятие до хвоста магнитопаузы второй импульс
val_DV5 = 1434;
dir_dV5 =1.0e+03 *[-2.139196402652551,   1.366800171828381,   0.638846286935880];
dir_dV5=dir_dV5/norm(dir_dV5);

T_int = 6.762583512360826e+05;
tdV_5 = tdV_4+T_int;
%Корректировка наклонения, чтобы точно в хвост магнитосферы попадать
val_DV6 = 52;
dir_dV6 =[-0.103171264325036   0.897977155531577  -2.269032961984546];
dir_dV6=dir_dV6/norm(dir_dV6);
T_long=2.804723854961529e+06;
tdV_6 = tdV_5+0.9*T_long;
%Спускаемся до Орбиты Тритона
val_DV7 = -263.2;
dir_dV7 =1.0e+03 *[1.170254764220188  -0.855716254550331  -0.391215977525669];
dir_dV7=dir_dV7/norm(dir_dV7);

T_long10_ap=5.170850471987467e+07;
tdV_7 = T_long10_ap+T_long*0.97;

%Устраиваем рандеву
val_DV8 = +24.56;
dir_dV8 =1.0e+03 *[-4.002190456180498   3.523715235508035   1.586704707928871];
dir_dV8=dir_dV8/norm(dir_dV8);

T_corr=1.243502576186567e+06;
T_enc=2.0590*T_corr;
tdV_8 = tdV_7+T_corr;

%Резонанс 4:1
val_DV9 = -82.75;
dir_dV9 =1.0e+03 *[-4.175713651185919   3.138992832089351   1.470355312252891];
dir_dV9=dir_dV9/norm(dir_dV9);

tdV_9 = tdV_8+T_enc;

%Резонанс 4:1 второй
val_DV10 = -234;
dir_dV10 =1.0e+03 *[-4.220086735027001   3.263528068258492   1.550851347342127];
dir_dV10=dir_dV10/norm(dir_dV10);

tdV_10 = tdV_9+TTr4_1;

%Тёмная зона
val_DV11 = -4;
dir_dV11 =1.0e+03 *[1.110697712707899  -0.744824370192724  -0.326725931032282];
dir_dV11=dir_dV11/norm(dir_dV11);

tdV_11 = tdV_10+TTr2_1;
% val_DVres = -111.7;
% dir_dVres =1.0e+03 *[  -4.474323061004518   3.067352615996629   1.350263920399283];
% dir_dVres=dir_dVres/norm(dir_dVres);
% 
% tdVres = tdV_5+0.833*TTr;
% 
val_DVend = -950;
dir_dVend =1.0e+02 *[9.299300014176843  -4.824496278194602  -2.935576262556574];
dir_dVend=dir_dVend/norm(dir_dVend);

tdVend = tdV_11+2.386833104641877e+06;

% dV=[dir_dV1*val_DV1;[0, 0, 0];dir_dV2*val_DV2;dir_dV3*val_DV3;dir_dV4*val_DV4;dir_dV5*val_DV5;dir_dVres*val_DVres;dir_dVend*val_DVend;];
% tdV=[t0_dist*0.1,t0_dist*0.98,tdV_2,tdV_3,tdV_4,tdV_5,tdVres,tdVend];
% lsp=[100;100;10000;10000;10000;10000;10000;10000;10000];
% 
% dV=[dir_dV1*val_DV1;[0, 0, 0];dir_dV2*val_DV2;dir_dV3*val_DV3;];
% tdV=[t0_dist*0.1,t0_dist*0.98,tdV_2,tdV_3];
% lsp=[100;100;10000;10000;10000;];
dV=[dir_dV1*val_DV1;[0, 0, 0];dir_dV2*val_DV2;dir_dV3*val_DV3;[0, 0, 0];...
    dir_dV4*val_DV4;dir_dV5*val_DV5;dir_dV6*val_DV6;[0, 0, 0];...
    dir_dV7*val_DV7;dir_dV8*val_DV8;dir_dV9*val_DV9;dir_dV10*val_DV10;dir_dV11*val_DV11;dir_dVend*val_DVend;];
tdV=[t0_dist*0.1,t0_dist*0.98,tdV_2,tdV_3,t_10T,tdV_4,tdV_5,tdV_6,T_long10_ap,tdV_7,tdV_8,tdV_9,tdV_10,tdV_11,tdVend];
lsp=[100;100;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;];%на один больше числа импульсов


% dV=[dir_dV1*val_DV1;[0, 0, 0];dir_dV2*val_DV2;dir_dV3*val_DV3;[0, 0, 0];...
%     dir_dV4*val_DV4;dir_dV5*val_DV5;dir_dV6*val_DV6;[0, 0, 0];...
%     dir_dV7*val_DV7;dir_dV8*val_DV8;dir_dV9*val_DV9;dir_dV10*val_DV10;dir_dV11*val_DV11;dir_dVend*val_DVend;];
% tdV=[t0_dist*0.1,t0_dist*0.98,tdV_2,tdV_3,t_10T,tdV_4,tdV_5,tdV_6,T_long10_ap,tdV_7,tdV_8,tdV_9,tdV_10];
% lsp=[100;100;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;];%на один больше числа импульсов


% dV=[dir_dV7*val_DV7;dir_dV8*val_DV8;dir_dV9*val_DV9;dir_dV10*val_DV10;dir_dV11*val_DV11;dir_dVend*val_DVend;];
% tdV=[tdV_7,tdV_8,tdV_9,tdV_10,tdV_11,tdVend];
% lsp=[10000;10000;10000;10000;10000;10000;10000;];%на один больше числа импульсов



rr_10T=1.0e+07 *[1.993305859334368;2.221608322826301;0.791584046491936];
VV_10T=1.0e+04 *[1.588094205170606;-1.178573767259988;-0.536972261163208];

rr_10longT=1.0e+09 *[1.175293787087023;1.045202900038859; 0.389891072126957];
VV_10longT=1.0e+03 *[1.031460770410799;-0.994856568367942;-0.443671157865987];

%y0=[rr_10longT;VV_10longT];

% 
% dV=[dir_dVres*val_DVres;dir_dVend*val_DVend;];
% tdV=[tdVres,tdVend];
% lsp=[10000;10000;10000];
% 
% r0=1.0e+08 *[  -2.207005142152339  -2.586411449673645  -0.982428788551360];
% V0=1.0e+03 *[  -4.424602685343921   3.034713680448864   1.335840878747880];
% y0=[r0,V0]';
%tspan=[0,tdVend+TTr3_1];

tspan=[0,tdVend+TTr3_1];

% dV=[];
% tdV=[];
% lsp=[100];
%Интегрирование со знанием импульсов
[t, y] = complexIntegration(y0, dV, tdV, tspan, lsp);
rr=y(:,1:3);
VV=y(:,4:6);
%Визуализация
%http://www.sai.msu.ru/neb/nss/html/multisat/nssreq8hr.htm
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
R_Neptune = 1;
[PoleNx,PoleNy,PoleNz]= sph2cart(RA_North,DE_Nort,R_Neptune*4);
[x,y,z] = sphere(50);
surf(R_Neptune*x, R_Neptune*y, R_Neptune*z,'HandleVisibility','off');
set(gca,'FontSize',14)
hold on;


plot3([R_Neptune*20*SunLine(1) R_Neptune*25*SunLine(1)], ...
    [R_Neptune*20*SunLine(2) R_Neptune*25*SunLine(2)],...
    [R_Neptune*20*SunLine(3) R_Neptune*25*SunLine(3)], 'y-x', 'LineWidth', 3,'DisplayName','Магнитопауза');

plot3([-R_Neptune*50*SunLine(1) -R_Neptune*80*SunLine(1)], ...
    [-R_Neptune*50*SunLine(2) -R_Neptune*80*SunLine(2)],...
    [-R_Neptune*50*SunLine(3) -R_Neptune*80*SunLine(3)], 'y-x', 'LineWidth', 3,'DisplayName','Хвост магнитосферы');

ind_phase_2=40200;
ind_phase_3=80200;
%Первая фаза
plot3(rr(1:ind_phase_2, 1)/RN, rr(1:ind_phase_2, 2)/RN, rr(1:ind_phase_2, 3)/RN, 'g', 'LineWidth', 1,'DisplayName','Первый этап');
%Вторая фаза
plot3(rr(ind_phase_2:ind_phase_3, 1)/RN, rr(ind_phase_2:ind_phase_3, 2)/RN, rr(ind_phase_2:ind_phase_3, 3)/RN, 'r', 'LineWidth', 1,'DisplayName','Второй этап');
%Третья фаза
plot3(rr(ind_phase_3:end, 1)/RN, rr(ind_phase_3:end, 2)/RN, rr(ind_phase_3:end, 3)/RN, 'b', 'LineWidth', 1,'DisplayName','Третий этап');
plot3([PoleNx; -PoleNx], [PoleNy; -PoleNy], [PoleNz; -PoleNz], 'k', 'LineWidth', 1,'DisplayName','Ось вращения Нептуна');


%plot3(rr(100, 1)/ae, rr(100, 2)/ae, rr(100, 3)/ae, 'rO', 'LineWidth', 1.5);
ind1=lsp(1)+lsp(2)+1;
ind2=lsp(1)+lsp(2)+lsp(3);
minimum = min(vecnorm(rr(ind1:ind2,:),2,2));
ind_per = find(minimum == vecnorm(rr(ind1:ind2,:),2,2))+lsp(1)+lsp(2);
%plot3(rr(ind_per, 1)/RN, rr(ind_per, 2)/RN, rr(ind_per, 3)/RN, 'rO', 'LineWidth', 1.5);
% ind_3_man=lsp(1)+lsp(2)+lsp(3)+lsp(4);
% plot3(rr(ind_3_man, 1)/RN, rr(ind_3_man, 2)/RN, rr(ind_3_man, 3)/RN, 'rO', 'LineWidth', 1.5);
% ind_4_man=lsp(1)+lsp(2)+lsp(3)+lsp(4)+lsp(5);
% plot3(rr(ind_4_man, 1)/RN, rr(ind_4_man, 2)/RN, rr(ind_4_man, 3)/RN, 'rO', 'LineWidth', 1.5);
% ind_5_man=lsp(1)+lsp(2)+lsp(3)+lsp(4)+lsp(5)+lsp(6);
% plot3(rr(ind_5_man, 1)/RN, rr(ind_5_man, 2)/RN, rr(ind_5_man, 3)/RN, 'rO', 'LineWidth', 1.5);
% ind_6_man=lsp(1)+lsp(2)+lsp(3)+lsp(4)+lsp(5)+lsp(6)+lsp(7);
% plot3(rr(ind_6_man, 1)/RN, rr(ind_6_man, 2)/RN, rr(ind_6_man, 3)/RN, 'rO', 'LineWidth', 1.5);
%ind_7_man=lsp(1)+lsp(2)+lsp(3)+lsp(4)+lsp(5)+lsp(6)+lsp(7)+lsp(8);
%ind_7_man=40000;
%plot3(rr(ind_7_man, 1)/RN, rr(ind_7_man, 2)/RN, rr(ind_7_man, 3)/RN, 'rO', 'LineWidth', 1.5);

%plot3([rr(1, 1)/ae;rr(1, 1)/ae+V0(1)*1e-7], [rr(1, 2)/ae;rr(1, 2)/ae+V0(2)*1e-7], [rr(1, 3)/ae;rr(1, 3)/ae+V0(3)*1e-7], 'r', 'LineWidth', 1);

%plot3([rr(1, 1)/ae;rr(1, 1)/ae+dV(1)*1e-7], [rr(1, 2)/ae;rr(1, 2)/ae+dV(2)*1e-7], [rr(1, 3)/ae;rr(1, 3)/ae+dV(3)*1e-7], 'r', 'LineWidth', 1);

%plot3(rr(end, 1)/RN, rr(end, 2)/RN, rr(end, 3)/RN, 'kO', 'LineWidth', 1.5);

%Тритон
keplerT = [61.315381532    354532.843 0.00000000 111.935349 308.921483   0.000000 219.270622];%31.3.2050
rrT = arrayfun(@(t)TritonR(t,keplerT), linspace(0,5.90*24*3600,500),'UniformOutput',false);
rrT = cell2mat(rrT)';
plot3(rrT(:, 1)/RN, rrT(:, 2)/RN, rrT(:, 3)/RN, 'b', 'LineWidth', 2.5,'DisplayName','Орбита Тритона');

rrT_real = arrayfun(@(t)TritonR(t,keplerT), t'-t0_dist,'UniformOutput',false);
rrT_real = cell2mat(rrT_real)';
min(vecnorm(rrT_real-rr,2,2))-RT
rT_1 = TritonR(tspan(2),keplerT);
%plot3(rrT_real(end,1)/RN, rrT_real(end,2)/RN, rrT_real(end,3)/RN, 'bO', 'LineWidth', 2.5);

%Нереида - нет эфемерид на дату

%Наяда
keplerN3 = [1221.774432464     48233.100 0.00032800  46.434045 236.010000 266.936289  22.820489];%31.3.2050
rrN3 = arrayfun(@(t)TritonR(t,keplerN3), linspace(0,2*24*3600,1000),'UniformOutput',false);
rrN3 = cell2mat(rrN3)';
plot3(rrN3(:, 1)/RN, rrN3(:, 2)/RN, rrN3(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');

%Таласса
keplerN4 = [1155.188396002     50069.200 0.00015600  46.859520  18.089245 354.055519  29.206508];%31.3.2050
rrN4 = arrayfun(@(t)TritonR(t,keplerN4), linspace(0,2*24*3600,1000),'UniformOutput',false);
rrN4 = cell2mat(rrN4)';
plot3(rrN4(:, 1)/RN, rrN4(:, 2)/RN, rrN4(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');

%Деспина
keplerN5 = [1074.933449933     52531.300 0.00013900  47.125758 168.020812 313.172633  29.286075];%31.3.2050
rrN5 = arrayfun(@(t)TritonR(t,keplerN5), linspace(0,2*24*3600,1000),'UniformOutput',false);
rrN5 = cell2mat(rrN5)';
plot3(rrN5(:, 1)/RN, rrN5(:, 2)/RN, rrN5(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');

%Галатея

keplerN6 = [839.456823581     61945.100 0.00012000  47.040284 179.420028 231.310371  29.328417];%31.3.2050
rrN6 = arrayfun(@(t)TritonR(t,keplerN6), linspace(0,2*24*3600,1000),'UniformOutput',false);
rrN6 = cell2mat(rrN6)';
plot3(rrN6(:, 1)/RN, rrN6(:, 2)/RN, rrN6(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');

%2004_N1

kepler2004_N1 = [378.887763843    105284.000 0.00000000  47.056000   7.635286   0.000000  29.373000];%31.3.2050
rrN2004_N1 = arrayfun(@(t)TritonR(t,kepler2004_N1), linspace(0,2*24*3600,1000),'UniformOutput',false);
rrN2004_N1 = cell2mat(rrN2004_N1)';
plot3(rrN2004_N1(:, 1)/RN, rrN2004_N1(:, 2)/RN, rrN2004_N1(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');

%Ларисса

keplerN7 = [648.892682712     73545.700 0.00138600  47.269028 304.463616 214.121357  29.114807];%31.3.2050
rrN7 = arrayfun(@(t)TritonR(t,keplerN7), linspace(0,2*24*3600,1000),'UniformOutput',false);
rrN7 = cell2mat(rrN7)';
plot3(rrN7(:, 1)/RN, rrN7(:, 2)/RN, rrN7(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');

%Протей

keplerN8 = [320.732634042    117646.000 0.00051000  47.626297 320.055732 124.639512  29.341616];%31.3.2050
rrN8 = arrayfun(@(t)TritonR(t,keplerN8), linspace(0,2*24*3600,1000),'UniformOutput',false);
rrN8 = cell2mat(rrN8)';
plot3(rrN8(:, 1)/RN, rrN8(:, 2)/RN, rrN8(:, 3)/RN, 'm', 'LineWidth', 1,'DisplayName','Орбиты спутников');


axis equal
%plot3(r0(1),r0(2),r0(3), 'go', 'LineWidth', 2.5);
title('Траектория КА')
xlabel('x, RN')
ylabel('y, RN')
zlabel('z, RN')
view(0,90)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlim([-50 65])
ylim([-50 65])
zlim([-40 40])
box off;
hold off;
legend;

%Проекция на Нептун
% figure(2);
% 
% ind_pr_1=20200;%Второй виток
% ind_pr_2=22400;%Второй виток
% [rr_proj_az,rr_proj_el] = arrayfun(@(t,x,y,z)projectTrajectory(t,[x,y,z]),...
%     t(ind_pr_1:ind_pr_2),rr(ind_pr_1:ind_pr_2, 1),rr(ind_pr_1:ind_pr_2, 2),rr(ind_pr_1:ind_pr_2, 3),'UniformOutput',false);
% rr_proj_az = cell2mat(rr_proj_az)';
% rr_proj_el = cell2mat(rr_proj_el)';
% plot(rr_proj_az(:), rr_proj_el(:), 'k.', 'LineWidth', 0.2);
% 
% xlim([-pi pi])
% ylim([-pi pi])
% ind_pr_1=20000;
% ind_pr_2=80000;
% %Вычисление в СО Нептуна
% rr_necc = rr(ind_pr_1:ind_pr_2, :);
% rr_N_SO = arrayfun(@(t,x,y,z)rotationNeptune(t,[x,y,z]),...
%     t(ind_pr_1:ind_pr_2),rr_necc(:,1),rr_necc(:,2),rr_necc(:,3),'UniformOutput',false);
% rr_N_SO = cell2mat(rr_N_SO')';
% t_N_SO = t(ind_pr_1:ind_pr_2)-t(ind_pr_1);
% t_N_SO(end)/24/3600
% %dlmwrite('r-N-gelio.csv',rr_N_SO,'precision',10)
% %dlmwrite('t-N-gelio.csv',t_N_SO,'precision',10)
% 
% n = cross( VV(ind_pr_2,:), rr(ind_pr_2,:));
% n = n/norm(n);
% 
% 
% 
% rr_necc_norm = vecnorm(rr_necc, 2, 2);
% tanalpha=1/65.5;
% 
% n_plus=n-tanalpha*rr_necc(1,:)/rr_necc_norm(1);
% n_plus=n_plus/norm(n_plus);
% n_minus=-n-tanalpha*rr_necc(1,:)/rr_necc_norm(1);
% n_minus=n_minus/norm(n_minus);
% 
% rr_necc_plus = rr_necc./rr_necc_norm;
% rr_necc_minus = rr_necc./rr_necc_norm;
% r_cos_plus=rr_necc_plus*n_plus';
% r_cos_minus=rr_necc_minus*n_minus';
% 
% for i = 1:length(r_cos_plus)
%     rr_necc_plus(i,:)=rr_necc_plus(i,:)+n_plus*r_cos_plus(i);
%     rr_necc_minus(i,:)=rr_necc_minus(i,:)+n_minus*r_cos_minus(i);
% end
% 
% rr_necc_plus_norm = vecnorm(rr_necc_plus, 2, 2);
% rr_necc_minus_norm = vecnorm(rr_necc_minus, 2, 2);
% 
% rr_necc_plus = rr_necc_plus./rr_necc_plus_norm;
% rr_N_SO_plus = arrayfun(@(t,x,y,z)rotationNeptune(t,[x,y,z]),...
%     t(ind_pr_1:ind_pr_2),rr_necc_plus(:,1),rr_necc_plus(:,2),rr_necc_plus(:,3),'UniformOutput',false);
% rr_N_SO_plus = cell2mat(rr_N_SO_plus')';
% dlmwrite('dir-plus-second-part.csv',rr_N_SO_plus,'precision',10)
% 
% rr_necc_minus = rr_necc_minus./rr_necc_plus_norm;
% rr_N_SO_minus = arrayfun(@(t,x,y,z)rotationNeptune(t,[x,y,z]),...
%     t(ind_pr_1:ind_pr_2),rr_necc_minus(:,1),rr_necc_minus(:,2),rr_necc_minus(:,3),'UniformOutput',false);
% rr_N_SO_minus = cell2mat(rr_N_SO_minus')';
% dlmwrite('dir-minus-second-part.csv',rr_N_SO_minus,'precision',10)
% 
% angle2vectors(rr_N_SO(1,:),rr_N_SO_minus(1,:))
angles=zeros(1,length(tdV));
for i = 1:length(tdV)
    t_c=tdV(i);
    %t_connection = t_Neptune+(t(end)-t0_dist)/24/3600;
    t_connection = t_Neptune+(t_c-t0_dist)/24/3600;
    r_Earth_Neptune = planetEphemeris(t_connection,'Earth','Neptune','430');
    r_Earth_Sun = planetEphemeris(t_connection,'Earth','Sun','430');
    angle_SEN = 180*angle2vectors(r_Earth_Neptune,r_Earth_Sun)/pi;
    angles(i)=angle_SEN;
end
