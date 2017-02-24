function s = lineShape(freq, fk, gamma, delta)
%	define the shape of individual line
g2 = gamma.^2;
fmf = fk - freq;
fpf = fk + freq;
z = (gamma - delta.*fmf)./(fmf.^2 + g2) + (gamma - delta.*fpf)./(fpf.^2 + g2);
s = z.*freq./fk;
