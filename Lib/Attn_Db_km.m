function [GammaDB, Gamma] = Attn_Db_km(freq,N)
%	attenuation in dB per km and neper per m
GammaDB = 0.1820*freq*N;
wk = 2*pi*freq/30.0;    % rad/cm
Gamma = 0.2*wk*N*1e-3; % m^-1
