function [c_O2, c_H2O] = loadO2H2O;
% load coefficients
load h2o.par;
c_H2O = h2o;
load o2.par;
c_O2 = o2;
