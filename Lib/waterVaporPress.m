function pwv = waterVaporPress(T, hum)
%	water vapor pressure in kPa
%	T in K
%	hum - rel. humid. 0<hum<1
teT = 300.0./T;
pwv = hum.*(2.408e11*teT.^5.*exp(-22.644*teT))/10.0;