function r_neptune = rotationNeptune(t, r)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
RA_North=2*pi*299.36/360; %https://en.wikipedia.org/wiki/Poles_of_astronomical_bodies
DE_North=2*pi*43.46/360;
% RA_North=pi/3;
% DE_North=pi/7;
T_neptune=0.6653*24*3600;%секунды
w_neptune = 2*pi/(T_neptune); %рад/с
ang0=0;
ang=w_neptune*t+ang0;
%[PoleNx,PoleNy,PoleNz]= sph2cart(RA_North,DE_North,1);
tmpq=my_eul_to_quat(RA_North, DE_North, ang,"ZXZs");
%Atr=quat2dcm(tmpq);
%r_neptune=Atr*r;
r_neptune = quatrotate(tmpq,r)';
x=r_neptune(2);
y=r_neptune(3);
z=r_neptune(1);

r_neptune=[x; y; z];
%q = quaternion(cos(ang/2), PoleNx*sin(ang/2), PoleNy*sin(ang/2), PoleNz*sin(ang/2));
end

