function r_icrf = mvt2icrf(r_mvt)

phi = (23+26/60+21.406/60/60)/180*pi;

M = [1,         0,       0;
     0,  cos(phi),-sin(phi);
     0,  sin(phi), cos(phi)];
 
r_icrf = M*r_mvt;

end