eps = 3.17+0.5*i;
ka = 0.05;
ang=(0:180)*pi/90;
[S,Sc,Ab]=SphereSc(sqrt(eps), ka, [1;1;0;0], ang);
polar(ang, (S(1,:)+S(2,:))/2,'g');
