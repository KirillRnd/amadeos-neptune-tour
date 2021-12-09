ra=1.355182519042239e+02*RN;
rp=27.057823528276312*RN;

a=(ra+rp)/2;
e=(ra-rp)/(ra+rp);
p=a*(1-e^2)
V_ell_ap = sqrt(mug/p)*(1-e);
V_ell_per = sqrt(mug/p)*(1+e);

TTr_4=2*pi*sqrt(a^3/mug)
TTr_4/24/3600