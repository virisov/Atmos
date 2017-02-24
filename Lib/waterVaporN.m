function n = waterVaporN(freq, T, press, pwv, c_H2O)
%	freq in GHz
%	T in K
%	press - air pressure in kPa
%	e - wv press., kPa
pd = press - pwv;
teT = 300.0./T;
k = 1:size(c_H2O,1);
s = S_H2O(k,teT,pwv,c_H2O);
f = F_H2O(k,freq,teT,pd,pwv,c_H2O);
z = sum(s.*f);
if z<0, z=0; end
nc = (1.13e-6*pwv.*pd.*teT.^3 + 3.57e-5*pwv.^2*teT.^10.8)*freq;
n = z + nc;

function y = S_H2O(k, teT, pwv, c_H2O)
y = c_H2O(k,2).*(pwv.*teT.^3.5*exp(c_H2O(k,3).*(1-teT)));

function y = F_H2O(k, freq, teT, pd, pwv, c_H2O)
g = c_H2O(k,4)*(pd.*teT.^0.6 + 4.8*pwv.*teT.^1.1)*1e-3;
y = lineShape(freq,c_H2O(k,1),g,0);