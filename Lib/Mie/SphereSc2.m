function [S,Sc,Ab] = SphereSc2(m, x, ang)

nc = ceil(x+4.05*(x^(1/3))+2);
%nc = 40;
n=(1:nc)';
E(n,1) = (2*n+1)./(n.*(n+1));
[p,t] = ALegendr(ang,nc);
[a,b] = ScatCoef(m,x,nc);

nn = 2*n+1;
Ex = 2*pi*real(sum(nn.*(a + b)));
Sc = 2*pi*sum(nn.*(abs(a).^2 + abs(b).^2));
Ab = Ex - Sc;

a = a.*E;
b = b.*E;

S1 = a'*p + b'*t;
S2 = a'*t + b'*p;
S11 = ((S2.*conj(S2))+(S1.*conj(S1)))/2;
S12 = ((S2.*conj(S2))-(S1.*conj(S1)))/2;
S33 = ((S1.*conj(S2))+(S2.*conj(S1)))/2;
S34 = i*((S1.*conj(S2))-(S2.*conj(S1)))/2;
%S = [I(1)*S11 + I(2)*S12; I(1)*S12 + I(2)*S11; I(3)*S33 + I(4)*S34; -I(3)*S34 + I(4)*S33];
S = [S11; S12; S33; S34]; 