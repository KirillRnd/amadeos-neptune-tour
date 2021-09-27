function r_Jupiter = rotationJupiter(t, r)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
RA_North=2*pi*268.06/360; %https://en.wikipedia.org/wiki/Poles_of_astronomical_bodies
DE_North=2*pi*64.50/360;
% RA_North=pi/3;
% DE_North=pi/7;
T_Jupiter=9.925*3600;%секунды
w_Jupiter = 2*pi/(T_Jupiter); %рад/с
ang0=0;
ang=w_Jupiter*t+ang0;
%[PoleNx,PoleNy,PoleNz]= sph2cart(RA_North,DE_North,1);
tmpq=my_eul_to_quat(RA_North, DE_North, ang,"ZXZs");
%Atr=quat2dcm(tmpq);
%r_neptune=Atr*r;
r_Jupiter = quatrotate(tmpq,r)';
x=r_Jupiter(2);
y=r_Jupiter(3);
z=r_Jupiter(1);

r_Jupiter=[x; y; z];
%q = quaternion(cos(ang/2), PoleNx*sin(ang/2), PoleNy*sin(ang/2), PoleNz*sin(ang/2));
end

