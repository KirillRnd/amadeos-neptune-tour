t_Earth=juliandate(2030,1,7);
t_Jupyter=juliandate(2032,6,27);


[r_E, V_E]=planetEphemeris(t_Earth,'SolarSystem','Earth','430');
r_E=r_E'*1e3;
V_E=V_E'*1e3;

[r_J, V_J]=planetEphemeris(t_Jupyter,'SolarSystem','Jupiter','430');
r_J=r_J'*1e3;
V_J=V_J'*1e3;

r1vec = r_E;  % Радиус-вектор точки вылета
r2vec = r_J;  % Радиус-вектор точки прилета

% гравитационный параметр центрального тела
muC = 132712.43994*(10^6)*(10^(3*3));

a=(norm(r_E)+norm(r_J))/2;

%tf = 2*pi*sqrt(a^3/muC);      % Время полета
tf = 2.5*365*24*3600;
m = 0;  % Число витков вокруг центрального тела
        % m = 0: без витков
        % m > 0, m < 0: с abs(m) числом витков, знак определяет одно из двух возможных решений
    
longway = -1; % специальный параметр, используется только его знак:
        % в зависимости от знака будет выбрана короткая или длинная дуга

[V1, V2, extremal_distances, exitflag] = klambert(r1vec, r2vec, tf, m, muC, longway);

norm(V2)/1000

dV1 = norm((V1-V_E))/1000;
dV2 = norm((V2-V_J))/1000;

i=180*asin(norm((V2-V_J))/norm(V_J))/pi;