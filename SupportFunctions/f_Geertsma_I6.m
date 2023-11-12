%% Auxiliaty Function or I_6
function [I6] = f_Geertsma_I6(q,r,R)

m = 4*R*r/(q^2 + (r+R)^2);
[K, E] = ellipke(m); 

I6 = q*m^1.5/(8*pi*R*(1-m)*(r*R)^1.5)*( 3*E + ...
    m*(R^2 - r^2 - q^2)/(r*R)* ( (1-m/2)/(1-m)*E - ...
    0.25*K));

end