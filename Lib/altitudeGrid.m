function z = altitudeGrid(Nz,zMiMa)
% exponentially increasing grid in Z starting from 0 to zMax=zMiMa(2);
% 
% zMin = zMiMa(1) = z(2) 
%
% old values:
%  N = 800;
%  
b = zMiMa(2)/zMiMa(1);  % zMax/zMin
u = fzero(@(u)((exp(u*Nz)-1)/(exp(u)-1)-b),5e-3); % find exponential distribution of z
z = zMiMa(1)*(exp(u*[0:Nz]) - 1)/(exp(u)-1);
