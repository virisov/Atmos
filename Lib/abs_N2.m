function [g, dgdt, dgdp] = abs_N2(temperature, pressure, frequency)
%  adopted from DOTLTR code VGI 8/17/2006
% INPUTS:
% temperature = TEMPERATURE (K)
% pressure  = PRESSURE (MB)
% frequency  = FREQUENCY (GHZ)
% OUTPUTS:
% g    = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR  (NEPER/KM)
% dgdt = DERIVATIVE OF g WITH RESPECT TO TEMPERATURE
% dgdp = DERIVATIVE OF g WITH RESPECT TO PRESSURE
% The derivatives were programmed and tested by R. Hill Aug. 2003

Th = 300./temperature;
g = (5.87e-14*frequency.^2).*(pressure.^2).*(Th.^4.5);
dgdt = g.*(-4.5./temperature);
dgdp = 2*g./pressure;

return
