function eps = iceEps(freq, T)
%	freq in GHz, T in K
teT = 300.0./T;
ai = (teT-0.171).*exp(17.0 - 22.1*teT);
bi = 1e-5*((0.233./(1-0.993./teT)).^2 + 6.33./teT - 1.31);
eps = 3.15 + i*(ai./freq + bi.*freq);
