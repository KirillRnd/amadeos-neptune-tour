function r_triton = rotationTriton(t, r)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
RA_North=2.253271008544791; %из положения Тритона
DE_North=-0.381889177677173;
% RA_North=pi/3;
% DE_North=pi/7;
T_Triton=-5.87*24*3600;%секунды
w_Triton = 2*pi/(T_Triton); %рад/с
ang0=pi*308.921483/180;
ang=w_Triton*t+ang0;
%[PoleNx,PoleNy,PoleNz]= sph2cart(RA_North,DE_North,1);
tmpq=my_eul_to_quat(RA_North, DE_North, ang,"ZXZs");
%Atr=quat2dcm(tmpq);
%r_neptune=Atr*r;
r_triton = quatrotate(tmpq,r)';
x=r_triton(2);
y=r_triton(3);
z=r_triton(1);

r_triton=[x; y; z];
%q = quaternion(cos(ang/2), PoleNx*sin(ang/2), PoleNy*sin(ang/2), PoleNz*sin(ang/2));
end

