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
TTr = 5.87*24*3600;
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
% initIntParams = struct;
% initIntParams.tspan = tspan;
% initIntParams.y0 = [Rrel;-Vrel];
% initIntParams.options = options;
% 
% hashInt = DataHash(initIntParams);
% fileHash = strcat('hashedInt/',hashInt,'.mat')
% try
%     sInt = load(fileHash);
%     t = sInt.t;
%     y = sInt.y;
% catch ME
%     if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
%         [t, y] = ode45(@(t,y) partialIntegrationSolarNeptune(t,y),initIntParams.tspan,initIntParams.y0,initIntParams.options);
%         save(fileHash,'t','y')
%     end   
% end
[t, y] = hashedIntegration(@(t,y) partialIntegrationSolarNeptune(t,y),tspan,[Rrel;-Vrel],options);


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
%Касаемся орбиты Тритона
val_DV5 = 270;
dir_dV5 =1.0e+03 *[-2.139196402652551,   1.366800171828381,   0.638846286935880];
dir_dV5=dir_dV5/norm(dir_dV5);

T_int = 6.762583512360826e+05;
tdV_5 = tdV_4+T_int;
%Энкаунтер с Тритоном
val_DV6 = 220.655;
dir_dV6 = 1.0e+03 *[3.640413305369320,  -3.005436955999429,  -1.357248293412794];
dir_dV6 = dir_dV6/norm(dir_dV6);

tdV_6 = tdV_5+3.5113*7.97e+05;
%Энкаунтер с Тритоном 2
val_DV7 = -13.212;
dir_dV7 = 1.0e+03 *[3.731151654288378  -3.445677784199575  -1.801600192258331];
dir_dV7 = dir_dV7/norm(dir_dV7);
tdV_7 = tdV_6+1.987*TTr;
%Энкаунтер с Тритоном 3
val_DV8 = 94.932;
dir_dV8 = 1.0e+03 *[3.589201605140489  -3.227523224302143  -1.734252863949240];
dir_dV8 = dir_dV8/norm(dir_dV8);
tdV_8 = tdV_7+2.998*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV9 = -2;
dir_dV9 = 1.0e+03 *[-1.536732032827026   1.645669390059489   0.695481357706951];
dir_dV9 = dir_dV9/norm(dir_dV9);
tdV_9 = tdV_8+TTr;
%Энкаунтер с Тритоном 4
val_DV10 = -37.031;
dir_dV10 = 1.0e+03 *[3.501853128951009  -3.332688519397423  -1.650696766597817];
dir_dV10 = dir_dV10/norm(dir_dV10);
tdV_10 = tdV_9+1.0008*TTr;

%Энкаунтер с Тритоном 5
val_DV11 = 94.962;
dir_dV11 = 1.0e+03 *[3.470507534290897  -3.467064704155914  -1.814914108907523];
dir_dV11 = dir_dV11/norm(dir_dV11);
tdV_11 = tdV_10+3.0021*TTr;

%Промежуточный манёвр коррекции для понижения перигея
val_DV12 = -2;
dir_dV12 = 1.0e+03 *[-1.379903382950022   1.760913384615123   0.714260792752529];
dir_dV12 = dir_dV12/norm(dir_dV12);
tdV_12 = tdV_11+1*TTr;

%Энкаунтер с Тритоном 6
val_DV13 = -85.215;
dir_dV13 = 1.0e+03 *[3.386189131116046  -3.423505398112130  -1.700954677300097];
dir_dV13 = dir_dV13/norm(dir_dV13);

tdV_13 = tdV_12+1.000*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV14 = -2;
dir_dV14 = 1.0e+03 *[-0.964983404031219   1.260277975967273   0.598201893055367];
dir_dV14 = dir_dV14/norm(dir_dV14);
tdV_14 = tdV_13+1.5*TTr;
%Энкаунтер с Тритоном 7
val_DV15 = 76.04;
dir_dV15 = 1.0e+03 *[3.554103382839291  -3.599150576443777  -1.849969492476835];
dir_dV15 = dir_dV15/norm(dir_dV15);

tdV_15 = tdV_14+1.501*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV16 = -2;
dir_dV16 = 1.0e+03 *[-1.200905580679059   1.853577372788878   0.746802427197531];
dir_dV16 = dir_dV16/norm(dir_dV16);
tdV_16 = tdV_15+1*TTr;
%Энкаунтер с Тритоном 8
val_DV17 = -35.31;
dir_dV17 = 1.0e+03 *[3.228560366961085  -3.758710742956700  -1.909290424775327];
dir_dV17 = dir_dV17/norm(dir_dV17);
tdV_17 = tdV_16+1.0017*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV18 = -2;
dir_dV18 = 1.0e+03 *[-0.749543351815548   1.363351309851669   0.622482208939737];
dir_dV18 = dir_dV18/norm(dir_dV18);
tdV_18 = tdV_17+1.5*TTr;
%Энкаунтер с Тритоном 9
val_DV19 = 68.735;
dir_dV19 = 1.0e+03 *[3.177607207961839  -3.716045190976224  -1.895540698866740];
dir_dV19 = dir_dV19/norm(dir_dV19);

tdV_19 = tdV_18+1.501*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV20 = -2;
dir_dV20 = 1.0e+03 *[-0.822339547209598   1.978184882115210   0.781341795704561];
dir_dV20 = dir_dV20/norm(dir_dV20);
tdV_20 = tdV_19+1*TTr;

%Энкаунтер с Тритоном 10
val_DV21 = -34.697;
dir_dV21 = 1.0e+03 *[2.922979928438882  -3.903984629634861  -1.974449515787414];
dir_dV21 = dir_dV21/norm(dir_dV21);
tdV_21 = tdV_20+1.001*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV22 = -3;
dir_dV22 = 1.0e+03 *[-0.822339547209598   1.978184882115210   0.781341795704561];
dir_dV22 = dir_dV22/norm(dir_dV22);
tdV_22 = tdV_21+1.5*TTr;
%Резкое понижение апогея
val_DV23 = -520;
dir_dV23 = 1.0e+03 *[2.866912394809400  -3.910147152744014  -1.988935669339936];
dir_dV23 = dir_dV23/norm(dir_dV23);
tdV_23 = tdV_22+1.5018*TTr;
%Энкаунтер с Тритоном на ближней стороне
val_DV24 = 4.05;
dir_dV24 = 1.0e+03 *[-3.975721997871332   1.966582308976535   0.973540177473503];
dir_dV24 = dir_dV24/norm(dir_dV24);
tdV_24 = tdV_23+0.37*TTr;

%Энкаунтер с Тритоном 2
val_DV25 = 194.963;
dir_dV25 = 1.0e+03 *[-4.007290489042338   1.968176427925677   1.338686065826040];
dir_dV25 = dir_dV25/norm(dir_dV25);
tdV_25 = tdV_24+1.13*TTr;
%Энкаунтер с Тритоном 3
val_DV26 = -34.668;
dir_dV26 = 1.0e+03 *[-4.502674513664267   1.969211944640743   1.694649410314734];
dir_dV26 = dir_dV26/norm(dir_dV26);
tdV_26 = tdV_25+2.001*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV27 = -2;
dir_dV27 = 1.0e+03 *[1.564344859757592  -0.034744351549018  -0.260039332373449];
dir_dV27 = dir_dV27/norm(dir_dV27);
tdV_27 = tdV_26+1.5*TTr;
%Энкаунтер с Тритоном 4
val_DV28 = 8.665;
dir_dV28 = 1.0e+03 *[-4.199742203424900   1.995019751200849   1.502180345913365];
dir_dV28 = dir_dV28/norm(dir_dV28);
tdV_28 = tdV_27+1.501*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV29 = -2;
dir_dV29 = 1.0e+03 *[1.564344859757592  -0.034744351549018  -0.260039332373449];
dir_dV29 = dir_dV29/norm(dir_dV29);
tdV_29 = tdV_28+1*TTr;

%Энкаунтер с Тритоном 5
val_DV30 = -79.815;
dir_dV30 = 1.0e+03 *[-4.669084050601275   1.531077325488367   1.751148650797024];
dir_dV30 = dir_dV30/norm(dir_dV30);
tdV_30 = tdV_29+1.001*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV31 = -2;
dir_dV31 = 1.0e+03 *[1.485198036408711   0.223323924845030  -0.181624846688641];
dir_dV31 = dir_dV31/norm(dir_dV31);
tdV_31 = tdV_30+1.5*TTr;
%Энкаунтер с Тритоном 6
val_DV32 = 33.93;
dir_dV32 = 1.0e+03 *[-4.610930688004008   1.259585259534775   1.663642318795774];
dir_dV32 = dir_dV32/norm(dir_dV32);
tdV_32 = tdV_31+1.502*TTr;
%Промежуточный манёвр коррекции для понижения перигея
val_DV33 = -2;
dir_dV33 = 1.0e+03 *[1.844931574494393   0.577789646214428  -0.064494179847868];
dir_dV33 = dir_dV33/norm(dir_dV33);
tdV_33 = tdV_32+1.0*TTr;
%Энкаунтер с Тритоном 7
val_DV34 = -71.18;
dir_dV34 = 1.0e+03 *[-4.886357311446234   1.095625322364601   1.874105144803822];
dir_dV34 = dir_dV34/norm(dir_dV34);
tdV_34 = tdV_33+1.001*TTr;

%Промежуточный манёвр коррекции для понижения перигея
val_DV35 = -2;
dir_dV35 = 1.0e+03 *[1.359869263491436   0.432434979231523  -0.100929076135750];
dir_dV35 = dir_dV35/norm(dir_dV35);
tdV_35 = tdV_34+1.5*TTr;

%Энкаунтер с Тритоном 8
val_DV36 = -36.38;
dir_dV36 = 1.0e+03 *[-4.714238345967947   0.746113909300086   1.781136954203919];
dir_dV36 = dir_dV36/norm(dir_dV36);
tdV_36 = tdV_35+1.501*TTr;

%Промежуточный манёвр коррекции для понижения перигея
val_DV37 = -2;
dir_dV37 = 1.0e+03 *[1.597682296888205   0.810245211671238   0.031714554878516];
dir_dV37 = dir_dV37/norm(dir_dV37);
tdV_37 = tdV_36+1.0*TTr;

%Энкаунтер с Тритоном 9
val_DV38 = 34.66;
dir_dV38 = 1.0e+03 *[-4.849444019006159   0.380798450010796   1.963948941729103];
dir_dV38 = dir_dV38/norm(dir_dV38);
tdV_38 = tdV_37+1.001*TTr;

%Промежуточный манёвр коррекции для понижения перигея
val_DV39 = -5;
dir_dV39 = 1.0e+03 *[1.145027804185982   0.636967746385523   0.012195155725596];
dir_dV39 = dir_dV39/norm(dir_dV39);
tdV_39 = tdV_38+1.5*TTr;

%Энкаунтер с Тритоном 10
val_DV40 = 43.13;
dir_dV40 = 1.0e+03 *[-4.693656963305115   0.152598896172514   1.899475367422481];
dir_dV40 = dir_dV40/norm(dir_dV40);
tdV_40 = tdV_39+1.501*TTr;

%Промежуточный манёвр коррекции для понижения перигея
val_DV41 = -3;
dir_dV41 = 1.0e+03 *[1.286248369996431   0.957740122165524   0.161407782726512];
dir_dV41 = dir_dV41/norm(dir_dV41);
tdV_41 = tdV_40+1.0*TTr;

%Энкаунтер с Тритоном 11
val_DV42 = -37.44;
dir_dV42 = 1.0e+03 *[-4.859480978750873  -0.140033663489182   2.127578082770063];
dir_dV42 = dir_dV42/norm(dir_dV42);
tdV_42 = tdV_41+1.001*TTr;

%Промежуточный манёвр коррекции для понижения перигея
val_DV43 = -5;
dir_dV43 = 1.0e+02 *[9.743009007621941   7.195609827225232   0.950261663383251];
dir_dV43 = dir_dV43/norm(dir_dV43);
tdV_43 = tdV_42+1.5*TTr;

%Уход от Тритона
val_DV44 = 200;
dir_dV44 = 1.0e+03 *[-4.942918239541855  -0.345181963694507   2.504921130466467];
dir_dV44 = dir_dV44/norm(dir_dV44);
tdV_44 = tdV_43+1.500*TTr;

%Повышение перигея орбиты
val_DV45 = 300;
dir_dV45 = 1.0e+02 *[3.865663443364614   3.211216088614391   0.567747237792116];
dir_dV45 = dir_dV45/norm(dir_dV45);
tdV_45 = tdV_44+5.5*TTr;

%Коррекция наклонения
val_DV46 = 70;
dir_dV46 = 1.0e+00 *[1.080462483200164  -1.637966332194327   1.925823563862027];
dir_dV46 = dir_dV46/norm(dir_dV46);
tdV_46 = tdV_45+7*TTr;

TLastEll = 80*24*3600;
%Падение на Нептун
val_DV47 = -850;
dir_dV47 = 1.0e+02 *[8.307335450106432   5.379627996746832  -0.380315690546450];
dir_dV47 = dir_dV47/norm(dir_dV47);
tdV_47 = tdV_46+4.5*TLastEll;

tdV_48 = tdV_47+0.5*TLastEll;

dV=[dir_dV1*val_DV1;[0, 0, 0];dir_dV2*val_DV2;dir_dV3*val_DV3;[0, 0, 0];...
    dir_dV4*val_DV4;dir_dV5*val_DV5;dir_dV6*val_DV6;dir_dV7*val_DV7;dir_dV8*val_DV8;...
    dir_dV9*val_DV9;dir_dV10*val_DV10;dir_dV11*val_DV11;dir_dV12*val_DV12;dir_dV13*val_DV13;...
    dir_dV14*val_DV14;dir_dV15*val_DV15;dir_dV16*val_DV16;dir_dV17*val_DV17;dir_dV18*val_DV18;...
    dir_dV19*val_DV19;dir_dV20*val_DV20;dir_dV21*val_DV21;dir_dV22*val_DV22;dir_dV23*val_DV23;...
    dir_dV24*val_DV24;dir_dV25*val_DV25;dir_dV26*val_DV26;dir_dV27*val_DV27;dir_dV28*val_DV28;...
    dir_dV29*val_DV29;dir_dV30*val_DV30;dir_dV31*val_DV31;dir_dV32*val_DV32;dir_dV33*val_DV33;...
    dir_dV34*val_DV34;dir_dV35*val_DV35;dir_dV36*val_DV36;dir_dV37*val_DV37;dir_dV38*val_DV38;...
    dir_dV39*val_DV39;dir_dV40*val_DV40;dir_dV41*val_DV41;dir_dV42*val_DV42;dir_dV43*val_DV43;...
    dir_dV44*val_DV44;dir_dV45*val_DV45;dir_dV46*val_DV46;dir_dV47*val_DV47;];
tdV=[t0_dist*0.1,t0_dist*0.98,tdV_2,tdV_3,t_10T,tdV_4,tdV_5,tdV_6,tdV_7,tdV_8,tdV_9,tdV_10,tdV_11,...
    tdV_12,tdV_13,tdV_14,tdV_15,tdV_16,tdV_17,tdV_18,tdV_19,tdV_20,tdV_21,tdV_22,tdV_23,tdV_24,...
    tdV_25,tdV_26,tdV_27,tdV_28,tdV_29,tdV_30,tdV_31,tdV_32,tdV_33,tdV_34,tdV_35,tdV_36,tdV_37,...
    tdV_38,tdV_39,tdV_40,tdV_41,tdV_42,tdV_43,tdV_44,tdV_45,tdV_46,tdV_47];
lsp=[100;100;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;...
    10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;...
    10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;...
    10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;10000;];%на один больше числа импульсов

disp(['Затраты хар. скорости ', num2str(sum(vecnorm(dV,2,2)),'%.2f\n'), 'м/с'])
tspan=[0,tdV_48];

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
ind_phase_3=230200;
ind_phase_4=440200;
%plot3(rr(end-10000:end, 1)/RN, rr(end-10000:end, 2)/RN, rr(end-10000:end, 3)/RN, 'g', 'LineWidth', 1,'DisplayName','Вся траектория');
%plot3(rr(:, 1)/RN, rr(:, 2)/RN, rr(:, 3)/RN, 'g', 'LineWidth', 1,'DisplayName','Вся траектория');
%plot3(rr(end-10000:end, 1)/RN, rr(end-10000:end, 2)/RN, rr(end-10000:end, 3)/RN, 'r', 'LineWidth', 1,'DisplayName','Вся траектория');
%Первая фаза
plot3(rr(1:ind_phase_2, 1)/RN, rr(1:ind_phase_2, 2)/RN, rr(1:ind_phase_2, 3)/RN, 'k', 'LineWidth', 1,'DisplayName','Первый этап');
%Вторая фаза
plot3(rr(ind_phase_2:ind_phase_3, 1)/RN, rr(ind_phase_2:ind_phase_3, 2)/RN, rr(ind_phase_2:ind_phase_3, 3)/RN, 'r', 'LineWidth', 1,'DisplayName','Второй этап');
%Третья фаза
plot3(rr(ind_phase_3:ind_phase_4, 1)/RN, rr(ind_phase_3:ind_phase_4, 2)/RN, rr(ind_phase_3:ind_phase_4, 3)/RN, 'b', 'LineWidth', 1,'DisplayName','Третий этап');
%Четвёртая фаза
plot3(rr(ind_phase_4:end, 1)/RN, rr(ind_phase_4:end, 2)/RN, rr(ind_phase_4:end, 3)/RN, 'm', 'LineWidth', 1,'DisplayName','Четвёртый этап');

plot3([PoleNx; -PoleNx], [PoleNy; -PoleNy], [PoleNz; -PoleNz], 'k', 'LineWidth', 1,'DisplayName','Ось вращения Нептуна');


%plot3(rr(100, 1)/ae, rr(100, 2)/ae, rr(100, 3)/ae, 'rO', 'LineWidth', 1.5);
%ind1=lsp(1)+lsp(2)+1;
%ind2=lsp(1)+lsp(2)+lsp(3);
%minimum = min(vecnorm(rr(ind1:ind2,:),2,2));
%ind_per = find(minimum == vecnorm(rr(ind1:ind2,:),2,2))+lsp(1)+lsp(2);
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
d_tr = min(vecnorm(rrT_real(end-10000:end,:)-rr(end-10000:end,:),2,2))-RT;
%d_tr = min(vecnorm(rrT_real(1:end,:)-rr(1:end,:),2,2))-RT;
disp(['Минимальное расстояние до Тритона ', num2str(d_tr/1000,'%.2f\n'), 'km.'])

rT_1 = TritonR(tspan(2),keplerT);

%plot3(rrT_real(end,1)/RN, rrT_real(end,2)/RN, rrT_real(end,3)/RN, 'bO', 'LineWidth', 2.5);

[x,y,z] = sphere(20);
r_surf=RT/RN;
%surf(rrT_real(end,1)/RN+r_surf*x, rrT_real(end,2)/RN+r_surf*y, rrT_real(end,3)/RN+r_surf*z,'HandleVisibility','off');

%Нереида - нет эфемерид на дату

% %Наяда
% keplerN3 = [1221.774432464     48233.100 0.00032800  46.434045 236.010000 266.936289  22.820489];%31.3.2050
% rrN3 = arrayfun(@(t)TritonR(t,keplerN3), linspace(0,2*24*3600,1000),'UniformOutput',false);
% rrN3 = cell2mat(rrN3)';
% plot3(rrN3(:, 1)/RN, rrN3(:, 2)/RN, rrN3(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');
% 
% %Таласса
% keplerN4 = [1155.188396002     50069.200 0.00015600  46.859520  18.089245 354.055519  29.206508];%31.3.2050
% rrN4 = arrayfun(@(t)TritonR(t,keplerN4), linspace(0,2*24*3600,1000),'UniformOutput',false);
% rrN4 = cell2mat(rrN4)';
% plot3(rrN4(:, 1)/RN, rrN4(:, 2)/RN, rrN4(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');
% 
% %Деспина
% keplerN5 = [1074.933449933     52531.300 0.00013900  47.125758 168.020812 313.172633  29.286075];%31.3.2050
% rrN5 = arrayfun(@(t)TritonR(t,keplerN5), linspace(0,2*24*3600,1000),'UniformOutput',false);
% rrN5 = cell2mat(rrN5)';
% plot3(rrN5(:, 1)/RN, rrN5(:, 2)/RN, rrN5(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');
% 
% %Галатея
% 
% keplerN6 = [839.456823581     61945.100 0.00012000  47.040284 179.420028 231.310371  29.328417];%31.3.2050
% rrN6 = arrayfun(@(t)TritonR(t,keplerN6), linspace(0,2*24*3600,1000),'UniformOutput',false);
% rrN6 = cell2mat(rrN6)';
% plot3(rrN6(:, 1)/RN, rrN6(:, 2)/RN, rrN6(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');
% 
% %2004_N1
% 
% kepler2004_N1 = [378.887763843    105284.000 0.00000000  47.056000   7.635286   0.000000  29.373000];%31.3.2050
% rrN2004_N1 = arrayfun(@(t)TritonR(t,kepler2004_N1), linspace(0,2*24*3600,1000),'UniformOutput',false);
% rrN2004_N1 = cell2mat(rrN2004_N1)';
% plot3(rrN2004_N1(:, 1)/RN, rrN2004_N1(:, 2)/RN, rrN2004_N1(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');
% 
% %Ларисса
% 
% keplerN7 = [648.892682712     73545.700 0.00138600  47.269028 304.463616 214.121357  29.114807];%31.3.2050
% rrN7 = arrayfun(@(t)TritonR(t,keplerN7), linspace(0,2*24*3600,1000),'UniformOutput',false);
% rrN7 = cell2mat(rrN7)';
% plot3(rrN7(:, 1)/RN, rrN7(:, 2)/RN, rrN7(:, 3)/RN, 'm', 'LineWidth', 1,'HandleVisibility','off');
% 
% %Протей
% 
% keplerN8 = [320.732634042    117646.000 0.00051000  47.626297 320.055732 124.639512  29.341616];%31.3.2050
% rrN8 = arrayfun(@(t)TritonR(t,keplerN8), linspace(0,2*24*3600,1000),'UniformOutput',false);
% rrN8 = cell2mat(rrN8)';
% plot3(rrN8(:, 1)/RN, rrN8(:, 2)/RN, rrN8(:, 3)/RN, 'm', 'LineWidth', 1,'DisplayName','Орбиты спутников');
% 

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
xlim([-100 80])
ylim([-80 100])
zlim([-100 100])
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
ind_pr_1=1;
ind_pr_2=470000;
t_Neptune=juliandate(2050,3,31);
t_N_SO_sec=t(ind_pr_1:ind_pr_2)-t0_dist;
% %Вычисление в СО Нептуна
rr_necc = rr(ind_pr_1:ind_pr_2, :);
rr_N_SO = arrayfun(@(t,x,y,z)rotationNeptune(t,[x,y,z]),...
    t_N_SO_sec,rr_necc(:,1),rr_necc(:,2),rr_necc(:,3),'UniformOutput',false);
rr_N_SO = cell2mat(rr_N_SO')';
t_N_SO_days =t_Neptune+(t(ind_pr_1:ind_pr_2)-t0_dist)/24/3600;
%t_N_SO(end)/24/3600
dlmwrite('r-N-full.csv',rr_N_SO,'precision',10)
dlmwrite('t-full-days.csv',t_N_SO_days,'precision',10)

rrT_real = arrayfun(@(t)TritonR(t,keplerT), t_N_SO_sec','UniformOutput',false);
rrT_real = cell2mat(rrT_real)';

rrTrelative=rr_necc-rrT_real;
rr_T_SO = arrayfun(@(t,x,y,z)rotationTriton(t,[x,y,z]),...
    t_N_SO_sec,rrTrelative(:,1),rrTrelative(:,2),rrTrelative(:,3),'UniformOutput',false);
rr_T_SO = cell2mat(rr_T_SO')';
% 
dlmwrite('r-T-full.csv',rr_T_SO,'precision',10)
% 
% %Направление на Солнце
rr_N_S = planetEphemeris(t_N_SO_days,'Neptune','Sun');
rr_N_SO = arrayfun(@(t,x,y,z)rotationNeptune(t,[x,y,z]),...
     t_N_SO_sec,rr_N_S(:,1),rr_N_S(:,2),rr_N_S(:,3),'UniformOutput',false);
rr_N_SO = cell2mat(rr_N_SO')';
rr_T_SO = arrayfun(@(t,x,y,z)rotationTriton(t,[x,y,z]),...
    t_N_SO_sec,rr_N_S(:,1),rr_N_S(:,2),rr_N_S(:,3),'UniformOutput',false);
rr_T_SO = cell2mat(rr_T_SO')';
n = cross( VV(ind_pr_2,:), rr(ind_pr_2,:));
n = n/norm(n);
% 
% 
% 
ind_pr_1=450200;
ind_pr_2=470200;
rr_necc = rr(ind_pr_1:ind_pr_2, :);
rr_necc_norm = vecnorm(rr_necc, 2, 2);
tanalpha=1/135.5;

n_plus=n-tanalpha*rr_necc(1,:)/rr_necc_norm(1);
n_plus=n_plus/norm(n_plus);
n_minus=-n-tanalpha*rr_necc(1,:)/rr_necc_norm(1);
n_minus=n_minus/norm(n_minus);

rr_necc_plus = rr_necc./rr_necc_norm;
rr_necc_minus = rr_necc./rr_necc_norm;
r_cos_plus=rr_necc_plus*n_plus';
r_cos_minus=rr_necc_minus*n_minus';

for i = 1:length(r_cos_plus)
    rr_necc_plus(i,:)=rr_necc_plus(i,:)+n_plus*r_cos_plus(i);
    rr_necc_minus(i,:)=rr_necc_minus(i,:)+n_minus*r_cos_minus(i);
end

rr_necc_plus_norm = vecnorm(rr_necc_plus, 2, 2);
rr_necc_minus_norm = vecnorm(rr_necc_minus, 2, 2);

rr_necc_plus = rr_necc_plus./rr_necc_plus_norm;
rr_N_SO_plus = arrayfun(@(t,x,y,z)rotationNeptune(t,[x,y,z]),...
    t(ind_pr_1:ind_pr_2),rr_necc_plus(:,1),rr_necc_plus(:,2),rr_necc_plus(:,3),'UniformOutput',false);
rr_N_SO_plus = cell2mat(rr_N_SO_plus')';
dlmwrite('dir-plus-second-part.csv',rr_N_SO_plus,'precision',10)

rr_necc_minus = rr_necc_minus./rr_necc_plus_norm;
rr_N_SO_minus = arrayfun(@(t,x,y,z)rotationNeptune(t,[x,y,z]),...
    t(ind_pr_1:ind_pr_2),rr_necc_minus(:,1),rr_necc_minus(:,2),rr_necc_minus(:,3),'UniformOutput',false);
rr_N_SO_minus = cell2mat(rr_N_SO_minus')';
dlmwrite('dir-minus-second-part.csv',rr_N_SO_minus,'precision',10)
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
