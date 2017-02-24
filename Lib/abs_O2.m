function [g, dgdt,dgdp, dgdv] = abs_O2(temperature, pressure, vapor_density, frequency)
%  adopted for MATLAB from DOTLTR code VGI 8/17/2006
% Real function o2abs  calculates the oxygen (o2) resonant and nonresonant absorption.
%
% See notes below for description
% Converted June 2003 from Pascal to Fortran by
% Ron Richter ETL SET 
%
% pwr 10/28/88
% modified:  AJG 3/22/91 to allow for adjustment of:
%            1, O2 linewidth temperature exponent XX
%            2, 118 line strength
%            3, O2-N2 low temperature exponent
%            4, O2-N2 temperature exponent breakpoint
%            AJG 4/3/91 to include non-overlapping SMMW lines
%
% INPUTS:
% temperature (K)
% pressure (mb)
% vapor_density (g/m^3) (water vapor density - enters linewidth calculation
% due to greater broadening efficiency of H2O)
% frequency (GHz)
% OUTPUTS:
% g = air absorption coefficient due to oxygen , in (Np/km)
% dgdt = DERIVATIVE OF g WITH RESPECT TO TEMPERATURE
% dgdp = DERIVATIVE OF g WITH RESPECT TO PRESSURE
% dgdv = DERIVATIVE OF g WITH RESPECT TO water vapor density
%
% frequencies, strength, linewidths from:
%    Liebe, Radio Science v.20, p. 1069 (1985)
% except for n=1 from:
%    Setzer&Pickett, J. Chem. Phys. v.67, p.340 (1977)
% interference from SY3, with LAGM=.16
% reference for equations:
%    P.W. Rosenkranz, JQSRT v.39, p.287 (1988)
% line are arranged 1-,1+,3-,3+,etc.
% m.k.: frequency of 773.8387 o2 abs line was changed to 773.8397GHz
% as a result of comparison between JPL online
% catalog and Rosenkranz abs line centers

F = [118.7503d0, 56.2648d0, 62.4863d0, 58.4466d0, 60.3061d0, 59.5910d0, ...
      59.1642d0, 60.4348d0, 58.3239d0, 61.1506d0, 57.6125d0, 61.8002d0, ...
      56.9682d0, 62.4112d0, 56.3634d0, 62.9980d0, 55.7838d0, 63.5685d0, ...
      55.2214d0, 64.1278d0, 54.6712d0, 64.6789d0, 54.1300d0, 65.2241d0, ...
      53.5957d0, 65.7648d0, 53.0669d0, 66.3021d0, 52.5424d0, 66.8368d0, ...
      52.0214d0, 67.3696d0, 51.5034d0, 67.9009d0];

% width in GHz/bar
w300 = [1.630d0, 1.646d0, 1.468d0, 1.449d0, 1.382d0, 1.360d0,         ...
        1.319d0, 1.297d0, 1.266d0, 1.248d0, 1.221d0, 1.207d0,         ...
        1.181d0, 1.171d0, 1.144d0, 1.139d0, 1.110d0, 1.108d0,         ...
        1.079d0, 1.078d0, 1.050d0, 1.050d0, 1.020d0, 1.020d0,         ...
        1.000d0, 1.000d0, 0.970d0, 0.970d0, 0.940d0, 0.940d0,         ...
        0.920d0, 0.920d0, 0.890d0, 0.890d0];

S300 = [0.2936d-14, 0.8079d-15, 0.2480d-14, 0.2228d-14,               ...
        0.3351d-14, 0.3292d-14, 0.3721d-14, 0.3891d-14,               ...
        0.3640d-14, 0.4005d-14, 0.3227d-14, 0.3715d-14,               ...
        0.2627d-14, 0.3156d-14, 0.1982d-14, 0.2477d-14,               ...
        0.1391d-14, 0.1808d-14, 0.9124d-15, 0.1230d-14,               ...
        0.5603d-15, 0.7842d-15, 0.3228d-15, 0.4689d-15,               ...
        0.1748d-15, 0.2632d-15, 0.8898d-16, 0.1389d-15,               ...
        0.4264d-16, 0.6899d-16, 0.1924d-16, 0.3229d-16,               ...
        0.8191d-17, 0.1423d-16];

% interference in 1/bar
U = [  -0.0247d0,    0.2881d0,    -0.4290d0,     0.6848d0,           ...
       -0.7170d0,    0.8266d0,    -0.6032d0,     0.5664d0,           ...
       -0.2635d0,    0.1731d0,    -0.2414d0,     0.1738d0,           ...
       -0.0556d0,   -0.0048d0,    -0.0596d0,     0.0134d0,           ...
       -0.0920d0,    0.0541d0,    -0.1151d0,     0.0814d0,           ...
       -0.0706d0,    0.0415d0,    -0.0314d0,     0.0069d0,           ...
       -0.0066d0,   -0.0143d0,     0.0252d0,    -0.0428d0,           ...
        0.0579d0,   -0.0726d0,     0.0883d0,    -0.1002d0,           ...
        0.1165d0,   -0.1255d0];

V = [   0.0003d0,   -0.0069d0,     0.0238d0,    -0.0647d0,           ...
        0.0916d0,   -0.1413d0,     0.1858d0,    -0.2323d0,           ...
        0.2686d0,   -0.3039d0,     0.3536d0,    -0.3797d0,           ...
        0.4104d0,   -0.4277d0,     0.4750d0,    -0.4860d0,           ...
        0.5025d0,   -0.5079d0,     0.5514d0,    -0.5525d0,           ...
        0.5520d0,   -0.5520d0,     0.5520d0,    -0.5520d0,           ...
        0.5520d0,   -0.5520d0,     0.5520d0,    -0.5520d0,           ...
        0.5520d0,   -0.5520d0,     0.5520d0,    -0.5520d0,           ...
        0.5520d0,   -0.5520d0];
     
% {data for SMMW lines (from Liebe)
% center frequency in GHz}
F0 = [368.4983, 424.7631, 487.2494, 715.3931, 773.8397, 834.1453];
% A1 IN KHZ/MB}
A1 = [6.79d-6, 6.38d-5, 2.35d-5, 9.96d-6, 6.71d-5, 1.8d-5];
% A2 (UNITLESS)}
A2 = [0.020, 0.011, 0.011, 0.089, 0.079, 0.079];
% A3 IN GHZ/MB}
A3 = [1.92d-3, 1.916d-3,1.92d-3, 1.81d-3, 1.81d-3, 1.81d-3];

wb300 = 0.48;
XX    = 0.8;
TN2   = 210.0;
YY    = 0.8;
AA    = 1.0;
BB    = 0.0;
CC    = 0.0;

freq2 = frequency.^2;
  
%RC = 217;
RC = 216.49;

PRESWV   = vapor_density.*temperature/RC;
dPRESWV_t = vapor_density/RC;
dPRESWV_v = temperature/RC;
% dPRESWV_p = 0

PRESDA   = pressure - PRESWV;
dPRESDA_t = - dPRESWV_t;
dPRESDA_v = - dPRESWV_v;
% dPRESDA_p = 1.

TH   = 300./temperature;
dTH_t = -TH./temperature;

B   =   TH.^XX;
dB_t = (XX.*B./TH).*dTH_t;  % =XX*B/TH*( -TH/temperature)

%  AAA   = AA * (  (pressure/1013.0)**BB) * (    TH**CC) % note that AAA=1.
   AAA   = 1;
  dAAA_t = AA*(  (pressure/1013.0).^BB) .* ( (CC*TH.^CC./TH).*dTH_t ); % note that dAAA_t=0. because CC=0
% dAAA_p = AA * ( BB*(pressure/1013.0)**(BB-1.) /1013.0 * (TH**CC) Note that dAAA_p=0. because BB=0
  dAAA_p = AA*((BB*(pressure/1013.0).^BB )./pressure).*(TH.^CC);
% dAAA_v = 0

TH1   = 1.0 - TH;
dTH1_t = -dTH_t;

TH16   = 6.89526d-3 *  TH1;
dTH16_t = 6.89526d-3 * dTH1_t;

TH02 = TH.^0.2; 
THp2   =  PRESDA.*TH02 + 1.1*PRESWV.*TH;
dTHp2_t = dPRESDA_t.*TH02 + PRESDA.*((0.2*TH02./TH).*dTH_t) + 1.1*(dPRESWV_t.*TH + PRESWV.*dTH_t);
dTHp2_v = dPRESDA_v.*TH02 + 1.1*dPRESWV_v.*TH;
dTHp2_p = TH02; % dPRESDA_p = 1. , dTH_p = 0                    dPRESWV_p = 0

TH30 = TH.^3; 
TH3   =  PRESDA   .* TH30;
dTH3_t = dPRESDA_t .* TH30 +  PRESDA.*(3*(TH30./TH).*dTH_t);
dTH3_v = dPRESDA_v .* TH30;
dTH3_p =              TH30; % dPRESDA_p = 1.
  
Fp04 = 0.04191*frequency;
G = B;
dG_t = dB_t;
idx = find(temperature < TN2);
if ~isempty(idx)
   G(idx)   = ((300/TN2)^XX) * ((TN2./temperature(idx)).^YY);
	dG_t(idx) = G(idx).*(-YY./temperature(idx));
end
 DEN =   0.001*(0.22*  PRESDA.*B   + 0.78*PRESDA.*G + 1.1*PRESWV.*TH);
dDEN_t = 0.001*(0.22*(dPRESDA_t.*B + PRESDA.*dB_t) + 0.78*(dPRESDA_t.*G + PRESDA.*dG_t) + ...
                 1.1*(dPRESWV_t.*TH + PRESWV.*dTH_t) );
dDEN_v = 0.001*(0.22*(dPRESDA_v.*B ) + 0.78*(dPRESDA_v.*G ) + 1.1*(dPRESWV_v.*TH) );
dDEN_p = 0.001*(0.22*B + 0.78*G );% dPRESDA_p = 1, % dPRESWV_p = 0

 DFNR   = wb300 * DEN;
dDFNR_t = wb300 * dDEN_t;
dDFNR_v = wb300 * dDEN_v;
dDFNR_p = wb300 * dDEN_p;

denomSUM1   =  TH.*(freq2 + DFNR.^2);
ddenomSUM1_t = dTH_t.*(freq2 + DFNR.^2) + TH.*(2*DFNR.*dDFNR_t);
ddenomSUM1_v =  TH.*(2*DFNR.*dDFNR_v);
ddenomSUM1_p =  TH.*(2*DFNR.*dDFNR_p);

 SUM1   = 1.6d-17*freq2.*DFNR./denomSUM1;
dSUM1_t = 1.6d-17*freq2.*(dDFNR_t - DFNR.*ddenomSUM1_t./denomSUM1)./denomSUM1;
dSUM1_v = 1.6d-17*freq2.*(dDFNR_v - DFNR.*ddenomSUM1_v./denomSUM1)./denomSUM1;
dSUM1_p = 1.6d-17*freq2.*(dDFNR_p - DFNR.*ddenomSUM1_p./denomSUM1)./denomSUM1;
% outputs from here show that derivatives that contribute to SUM1, dSUM1, which are
% dDFNR, ddenomSUM1, dTH_t, dDEN, dPRESDA, dB_t, dG_t, dPRESWV, are all correct. 
 
for K = 1:34
   if mod(K,2)==1
      BFAC = exp(K*(K+1)*TH16);
      dBFAC_t = BFAC.*(K*(K+1)*dTH16_t);
   end

   DF  = w300(K) *  DEN;
	dDF_t = w300(K) * dDEN_t;
	dDF_v = w300(K) * dDEN_v;
	dDF_p = w300(K) * dDEN_p;

   DF2  = DF.^2;
	dDF2_t = 2.*DF.*dDF_t;
	dDF2_v = 2.*DF.*dDF_v;
	dDF2_p = 2.*DF.*dDF_p;


   Y = DEN.*(U(K) + V(K)*TH);
	dY_t = dDEN_t.*(U(K) + V(K)*TH) + DEN.*(V(K) * dTH_t);
   dY_v = dDEN_v.*(U(K) + V(K)*TH);
   dY_p = dDEN_p.*(U(K) + V(K)*TH);

   STR = S300(K)*BFAC; %these 3 lines would be more efficient within an ELSE of the IF below
	dSTR_t = S300(K)*dBFAC_t;
	dSTR_p = 0;
   
   if K==1
      STR = STR.*AAA;
      dSTR_t = dSTR_t.*AAA + STR.*dAAA_t;
      dSTR_p =               STR.*dAAA_p;
   end
   f1 = frequency - F(K);
   den = (f1.^2 + DF2);
   SF1 = (DF + f1.*Y)./den;
   dSF1_t = ((dDF_t + f1.*dY_t) - (DF + f1.*Y).*dDF2_t./den)./den;
   dSF1_v = ((dDF_v + f1.*dY_v) - (DF + f1.*Y).*dDF2_v./den)./den;
   dSF1_p = ((dDF_p + f1.*dY_p) - (DF + f1.*Y).*dDF2_p./den)./den;
   
   f2 = frequency + F(K);
   den = (f2.^2 + DF2);
   SF2 = (DF - f2.*Y)./den;
   dSF2_t = ((dDF_t - f2.*dY_t) - (DF + f2.*Y).*dDF2_t./den)./den;
   dSF2_v = ((dDF_v - f2.*dY_v) - (DF + f2.*Y).*dDF2_v./den)./den;
   dSF2_p = ((dDF_p - f2.*dY_p) - (DF + f2.*Y).*dDF2_p./den)./den;
   
   den = freq2./F(K).^2;
   SUM1 = SUM1 + STR.*(SF1 + SF2).*den;
   dSUM1_t = dSUM1_t + (dSTR_t.*(SF1 + SF2) + STR.*(dSF1_t + dSF2_t)).*den;
   dSUM1_v = dSUM1_v + (                      STR.*(dSF1_v + dSF2_v)).*den;
   dSUM1_p = dSUM1_p + (dSTR_p.*(SF1 + SF2) + STR.*(dSF1_p + dSF2_p)).*den;

end

dSUM1_t = (0.5034d12/pi)*(dSUM1_t.*TH3 + SUM1.*dTH3_t);
dSUM1_v = (0.5034d12/pi)*(dSUM1_v.*TH3 + SUM1.*dTH3_v);
dSUM1_p = (0.5034d12/pi)*(dSUM1_p.*TH3 + SUM1.*dTH3_p);
SUM1    = (0.5034d12/pi)*(SUM1.*TH3);
% output from a RETURN here shows that the dSUM1 are correct to here, BEWARE:
% SUM1 = ...  statement MUST be after the derivatives are calculated

% ADD NON-OVERLAPPING SMMW O2 LINES}

for K=1:6
   DF     = A3(K)*THp2;
   dDF_t  = A3(K)*dTHp2_t;
   dDF_v  = A3(K)*dTHp2_v;
   dDF_p  = A3(K)*dTHp2_p;
   
   DF2    = DF.^2;
   dDF2_t = 2*DF.*dDF_t;
   dDF2_v = 2*DF.*dDF_v;
   dDF2_p = 2*DF.*dDF_p;

   f = (frequency/F0(K))*A1(K);  % FORTRAN convention ?
   e = exp(A2(K)*TH1);
   STR = f.*DF.*TH3.*e;
   dSTR_t = f.*(dDF_t.*TH3 + DF.*dTH3_t + DF.*TH3.*(A2(K)*dTH1_t)).*e;
   dSTR_v = f.*(dDF_v.*TH3 + DF.*dTH3_v).*e;
	dSTR_p = f.*(dDF_p.*TH3 + DF.*dTH3_p).*e;

   f1sq = (frequency - F0(K)).^2;
   f2sq = (frequency + F0(K)).^2;

   den = f1sq + DF2;
   den2 = den.^2;
   SF1   = 1./den;
   dSF1_t = -dDF2_t./den2;
   dSF1_v = -dDF2_v./den2;
   dSF1_p = -dDF2_p./den2;
   
   den = f2sq + DF2;
   den2 = den.^2;
   SF2   = 1./den;
   dSF2_t = -dDF2_t./den2;
   dSF2_v = -dDF2_v./den2;
   dSF2_p = -dDF2_p./den2;

	SUMAND    = Fp04.*  STR  .*(SF1 + SF2);
	dSUMAND_t = Fp04.*(dSTR_t.*(SF1 + SF2) + STR.*(dSF1_t + dSF2_t));
	dSUMAND_v = Fp04.*(dSTR_v.*(SF1 + SF2) + STR.*(dSF1_v + dSF2_v));
	dSUMAND_p = Fp04.*(dSTR_p.*(SF1 + SF2) + STR.*(dSF1_p + dSF2_p));

   SUM1    =  SUM1   +  SUMAND;
   dSUM1_t = dSUM1_t + dSUMAND_t;
   dSUM1_v = dSUM1_v + dSUMAND_v;
   dSUM1_p = dSUM1_p + dSUMAND_p;
   
end

g = SUM1;
dgdt = dSUM1_t;
dgdv = dSUM1_v;
dgdp = dSUM1_p;

return


