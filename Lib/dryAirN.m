function n = dryAirN(freq, T, press, e, c_O2)
%	dry air attenuation
%	freq in GHz
%	teT = 300/T, T in K;
%	press in kPa
%	e - wv press., kPa
teT = 300.0./T;
k = 1:size(c_O2,1);
s = S_O2(k,teT,press,c_O2);
f = F_O2(k,freq,teT,press,e,c_O2);
z = sum(s.*f);
if z<0, z=0; end
g = 4.8e-3*(press + 1.1*e)*teT.^0.8;
if g<=0 
   trm1 = 0;
else
   trm1 = 6.14e-4/(g + freq^2./g);
end
trm2 = 1.4e-10*(1 - 1.2e-5*freq^1.5)*press*teT.^1.5;
n = (trm1 + trm2).*freq*press*teT.^2 + z;

function y = F_O2(k,freq,teT,press,e,c_O2)
gamma = c_O2(k,4).*(press.*teT.^(0.8 - c_O2(k,7)) + 1.1*e.*teT)*1e-3;
delta = c_O2(k,5).*(press.*teT.^c_O2(k,6))*1e-3;
%delta = (c_O2(k,5) + c_O2(k,6)*teT).*press.*teT.^0.8*1e-3;
y = lineShape(freq,c_O2(k,1),gamma,delta);

function y = S_O2(k,teT,press,c_O2)
y = c_O2(k,2).*press.*teT.^3.*exp(c_O2(k,3).*(1-teT))*1e-6;


