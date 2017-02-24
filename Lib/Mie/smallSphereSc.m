function [Sc,Ab] = smallSphereSc(eps,ka)
%
e = (eps-1)/(eps+2);
Sc = (8*pi/3)*ka^6*(abs(e))^2;
Ab = 4*pi*ka^3*imag(e);