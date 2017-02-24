ka = [0.05:0.02:2.0];
m = [2.0:0.1:8];

for i=1:length(ka)
    for j=1:length(m)
        [S,Sc,Ab]=SphereSc(m(j), ka(i), [1;1;0;0], 0);
        %Sx(i,j) = Sc(1)/(ka(i).^2);
        Sx(i,j) = Sc(1);
    end
end

