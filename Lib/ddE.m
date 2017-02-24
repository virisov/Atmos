function eps = ddE(sal, tempC, freq)
%
% model dielectric permittivity
%
% W.Ellison et al."New permittivity of seawater",
% Radio Sci.,vol.33, no.3, pp.639-648, 1998.
% sal in [0/00 (0.001)], temp in [C], freq in [GHz]

tm = cumprod(tempC*ones(1,5));

c1 = 0.086374 + 0.030606*tm(1) - 0.0004121*tm(2);
c2 = 0.077454 + 0.001687*tm(1) + 0.00001937*tm(2);
sigma = c1 + c2*sal;

a1 = 81.8200 - 6.0503e-2*tm(1) - 3.1661e-2*tm(2) + 3.1097e-3*tm(3) - 1.1791e-4*tm(4) + 1.4838e-6*tm(5);
a2 = 0.12544 + 9.4037e-3*tm(1) - 9.5551e-4*tm(2) + 9.0888e-5*tm(3) - 3.6011e-6*tm(4) + 4.7130e-8*tm(5);

eps0 = a1 + a2*sal;

b1 = 17.303    - 0.66651*tm(1) + 5.1482e-3*tm(2) + 1.2145e-3*tm(3) - 5.0325e-5*tm(4) + 5.8272e-7*tm(5);
b2 = -6.272e-3 + 2.357e-4*tm(1) + 5.075e-4*tm(2) - 6.3983e-5*tm(3) + 2.4630e-6*tm(4) - 3.0676e-8*tm(5);
   
tau = (b1 + b2*sal)*1e-12;

epsInf = 6.4587 - 0.04203*tm(1) - 0.0065881*tm(2) + 0.00064924*tm(3) - 1.2328e-5*tm(4) + 5.0433e-8*tm(5);

omega = 2*pi*freq*1e9;
denom = 1 + (omega*tau).^2;
epsR = epsInf + (eps0 - epsInf)./denom;
eStar = 8.8419e-12;  % [F/m]
epsI = omega*tau*(eps0 - epsInf)./denom + sigma./(eStar*omega);

eps = complex(epsR, epsI);
