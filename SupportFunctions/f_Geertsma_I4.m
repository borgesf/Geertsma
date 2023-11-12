%% Auxiliaty Function or I_4
function [I4] =  f_Geertsma_I4(q,r,R)

m = 4*R*r/(q^2 + (r+R)^2);
[K, E] = ellipke(m); 

I4 = (m^1.5)*(R^2 - r^2 - q^2)*E/(8*pi*(r*R)^1.5*R*(1-m)) + ...
    sqrt(m)*K/(2*pi*sqrt(r*R)*R);

end