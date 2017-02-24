function [g, dgdt] = abs_liquidH2O(temperature, water_content, frequency)
% INPUTS:
% temperature = TEMPERATURE (K)
% water_content  = WATER_CONTENT (g/m^3)
% frequency  = FREQUENCY (GHz)
% OUTPUTS:
% g    = ABSORPTION COEFFICIENT DUE TO LIQUID WATER  (NEPER/KM)
% dgdt = DERIVATIVE OF g WITH RESPECT TO TEMPERATURE

Rho = 1e6;  %    g/m^3
DT = 1; % K
wk = 2*pi*frequency/0.3;     % rad/m
eps = waterEps(frequency,temperature);
g = 3e3*wk.*(water_content/Rho).*imag((eps-1)./(eps+2));   % 1/km
eps1 = waterEps(frequency,temperature+DT);
dgdt = (3e3*wk.*(water_content/Rho).*imag((eps1-1)./(eps1+2)) - g)/DT;
end
