function eps = waterEps(freq, T)
%	freq in GHz
%	T in K
t = 300.0./T - 1;
e0 = 77.66 + 103.3*t;
e1 = 0.0671*e0;
e2 = 3.52;
g1 = 20.2 - 146.0*t + 316.0*t.^2;
g2 = 39.8*g1;
z = (e0-e1)./(1 + i*g1./freq) + (e1-e2)./(1 + i*g2./freq);
eps = e0 - z;
