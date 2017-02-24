function [g, dgdt,dgdp, dgdv] = abs_H2O(temperature, pressure, vapor_density, frequency)
%  adopted for MATLAB from DOTLTR code VGI 8/17/2006
% (* version 10.1 2/17/91 AJG
% rewritten for Pascal by MK, March 1997
% adapted from version 10, pwr
% differentiated and tested by R. Hill, Aug. 2003
% assumes vvw line shape

%outputs:
% abh2o = absorption coeff. in the atmosphere due to water vapor (Nepers/km)
% derivative of abh2o with respect to temperature: dabh2o_t  (Nepers/km/K)
% derivative of abh2o with respect to water vapor density: dabh2o_v  (Nepers/km/(grams/cubic meter))
% derivative of abh2o with respect to pressure dabh2o_p    (Nepers/km/mb)

%inputs:
% temperature (K)
% pressure (mb)
% RHO = water vapor density (g/m3), named vapor_density in o2abs
% frequency (GHz)
%
% line parameters from H. J. Liebe, Radio Science V.20(5), pp.1069-1089 (1985)
% updated in FREQUENZ V.41, pp. 31-36 (1987)
% includes 30 MPM H2O lines and appropriate continuum term for this line base
% 
S1 = [0.0109d0,  0.00011d0, 0.00007d0, 0.23000d0,  0.00464d0, 0.1540d0, 0.0001d0,...
      1.1900d0,  0.00044d0, 0.00637d0, 0.09210d0,  0.01940d0, 1.0600d0, 0.0330d0,...
      0.1280d0,  0.02530d0, 0.00374d0, 0.00125d0, 51.00000d0, 0.5090d0, 0.0274d0,...
     25.0000d0,  0.00130d0, 0.01330d0, 0.00550d0,  0.00380d0, 0.0183d0, 0.8560d0,...
      0.9160d0, 13.80000d0];

B2 = [2.143d0, 8.730d0, 8.347d0, 0.653d0, 6.156d0, 1.515d0, 9.802d0,...
      1.018d0, 7.318d0, 5.015d0, 3.561d0, 5.015d0, 1.370d0, 3.561d0,...
      2.342d0, 2.814d0, 6.693d0, 6.693d0, 0.114d0, 2.150d0, 7.767d0,...
      0.336d0, 8.113d0, 7.989d0, 7.845d0, 8.360d0, 5.039d0, 1.369d0,...
      1.842d0, 0.178d0];

W3 = [2.784d-3, 2.760d-3, 2.700d-3, 3.164d-3, 2.140d-3, 2.970d-3, 2.650d-3,...
      3.036d-3, 1.900d-3, 1.370d-3, 1.640d-3, 1.440d-3, 2.380d-3, 1.820d-3,...    
      1.980d-3, 2.490d-3, 1.150d-3, 1.190d-3, 3.000d-3, 2.230d-3, 3.000d-3,...
      2.860d-3, 1.410d-3, 2.860d-3, 2.860d-3, 2.640d-3, 2.340d-3, 2.530d-3,...
      2.400d-3, 2.860d-3];

FL = [22.235080d0,  67.813960d0, 119.995940d0, 183.310117d0, 321.225644d0, 325.152919d0, 336.187000d0, ...
     380.197372d0, 390.134508d0, 437.346667d0, 439.150812d0, 443.018295d0, 448.001075d0, 470.888947d0, ...
     474.689127d0, 488.491133d0, 503.568532d0, 504.482692d0, 556.936002d0, 620.700807d0, 658.006500d0, ...
     752.033227d0, 841.073593d0, 859.865000d0, 899.407000d0, 902.555000d0, 906.205524d0, 916.171582d0, ...
     970.315022d0, 987.926764d0];

RHO = vapor_density;

if RHO <= 0 
    g    = 0;
    dgdt = 0;
    dgdv = 0;
    dgdp = 0;
    return
end

%RC = 217;
RC = 216.49;

PVAP = RHO.*temperature/RC;   % 216.49?
dPVAP_t = RHO/RC;
dPVAP_v = temperature/RC;
dPVAP_p = 0.;

PDA = pressure - PVAP;
dPDA_t = - dPVAP_t;
dPDA_v = - dPVAP_v;
dPDA_p = 1.;

TI = 300.0./temperature;
dTI_t = -TI./temperature;

TI3 = TI.^3.5;
dTI3_t = (3.5*TI3./TI).*dTI_t;

TI1 = 1.0 - TI;
dTI1_t = - dTI_t;

TIP3 = PVAP.*TI3;
dTIP3_t = dPVAP_t.*TI3 + PVAP.*dTI3_t;
dTIP3_v = dPVAP_v.*TI3;
dTIP3_p = dPVAP_p.*TI3;

TI08 = TI.^0.8;
WFAC = PDA.*TI08 + 4.8*PVAP.*TI;
dWFAC_t = dPDA_t.*TI08 + 4.8*dPVAP_t.*TI +  PDA.*((0.8*TI08./TI).*dTI_t) + 4.8*PVAP.*dTI_t;
dWFAC_v = dPDA_v.*TI08 + 4.8*dPVAP_v.*TI;
dWFAC_p = dPDA_p.*TI08 + 4.8*dPVAP_p.*TI;

TI30  = 1.13d-8*frequency.*TI.^3;
TI105 = 3.57d-7*frequency.*TI.^10.5;
den = PDA.*TI30 + PVAP.*TI105;
SUM1 = PVAP.*den;
dSUM1_t = dPVAP_t.*den + PVAP.*(dPDA_t.*TI30 + dPVAP_t.*TI105 ...
					            + PDA.*(3*TI30./TI).*dTI_t + PVAP.*(10.5.*TI105./TI).*dTI_t);
dSUM1_v = dPVAP_v.*den + PVAP.*(dPDA_v.*TI30 + dPVAP_v.*TI105);
dSUM1_p = dPVAP_p.*den + PVAP.*(dPDA_p.*TI30 + dPVAP_p.*TI105);

for I=1:30
   WIDTH    = W3(I).*WFAC; 
   dWIDTH_t = W3(I).*dWFAC_t;
   dWIDTH_v = W3(I).*dWFAC_v;
   dWIDTH_p = W3(I).*dWFAC_p;

   W2   = WIDTH.^2;
   dW2_t = 2*WIDTH.*dWIDTH_t;
   dW2_v = 2*WIDTH.*dWIDTH_v;
   dW2_p = 2*WIDTH.*dWIDTH_p;

   f1 = (frequency - FL(I));
   f2 = (frequency + FL(I));
   e = S1(I)*exp(B2(I).*TI1);
   S   = TIP3.*e;
	dS_t = e.*(dTIP3_t  + TIP3 .*(B2(I).*dTI1_t));
	dS_v = e.*dTIP3_v;
   dS_p = e.*dTIP3_p;
   
   freq = frequency/FL(I);
   den1 = (f1.^2 + W2);
   den2 = (f2.^2 + W2);
   wf1 = freq.*WIDTH./den1;
   wf2 = freq.*WIDTH./den2;
   
   SUM1 = SUM1 + S.*(wf1 + wf2);
   dSUM1_t = dSUM1_t + (wf1 + wf2).*(dS_t + S.*dWIDTH_t./WIDTH) - S.*(wf1./den1 + wf2./den2).*dW2_t;
   dSUM1_v = dSUM1_v + (wf1 + wf2).*(dS_v + S.*dWIDTH_v./WIDTH) - S.*(wf1./den1 + wf2./den2).*dW2_v;
   dSUM1_p = dSUM1_p + (wf1 + wf2).*(dS_p + S.*dWIDTH_p./WIDTH) - S.*(wf1./den1 + wf2./den2).*dW2_p;
end

g    = 0.0419.*frequency.*SUM1;
dgdt = 0.0419.*frequency.*dSUM1_t;
dgdv = 0.0419.*frequency.*dSUM1_v;
dgdp = 0.0419.*frequency.*dSUM1_p;

return
