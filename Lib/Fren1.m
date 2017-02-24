function [rv, rh] = Fren1(eps,tet)
%
ct = cos(tet);
st = sin(tet);
s = sqrt(eps - st.^2);
rv = (eps*ct - s)./(eps*ct + s);
rh = (ct - s)./(ct + s);