%% Auxiliary Function for I_2
function [I2] = f_Geertsma_I2(q,r,R)

m = 4*R*r/(q^2 + (r+R)^2);
[K, E] = ellipke(m); 

I2 = (q*sqrt(m)/(2*pi*(r*R)^1.5))*( ((1 - m/2)/(1-m))*E - ...
    K);

end