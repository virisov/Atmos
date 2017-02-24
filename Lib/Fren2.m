function r = Fren2(eps, tet)
%	power Fresnel coefficients
[rv, rh] = Fren1(eps, tet);
r = abs([rv, rh]).^2;
