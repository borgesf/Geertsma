%% Auxiliaty Function or I_7
function [I7] = f_Geertsma_I7(q,r,R)

m = 4*R*r/(q^2 + (r+R)^2);
[K, E] = ellipke(m); 

I7 = (-1/q)*f_Geertsma_I2(q,r,R) + q^2*m^1.5/(8*pi*(1-m)*(r*R)^2.5)*(...
     ((1 - m + m^2)/(1-m))*E - (1 - 0.5*m)*K);

end