function [press, T, potT] = standardAtm(h, tCsurf, G)
%	standard atmosphere
%	altitude h in km
%	press in kPa
%	T and potT in K
if nargin<3
    G = 9.80665;
end
if nargin<2
    tCsurf = 15.0;
end

%Lps = [-6.5, 0, 1, 2.8, 0, -2.8, -2];  % C/km
%Hgeopot = [-0.6, 11, 20, 32, 47, 51, 71];  % km
%Tb = [19, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5]; % K


P0 = 1013.25;
T0 = 273.16 + tCsurf;
Lapse = 6.5;
Rd = 287.04;
Alfa = G/(1e-3*Lapse*Rd);
T = T0 - Lapse*h;
press = P0*(T/T0).^Alfa;
potT = T.*(P0./press).^0.286;
press = press/10; % *0.1 mB->kPa
end
