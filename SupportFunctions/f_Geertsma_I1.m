%% Auxiliary Function for I_1
function [I1] = f_Geertsma_I1(q,r,R)

m = 4*R*r/(q^2 + (r+R)^2);
[K, E] = ellipke(m); 

I1 = (2/(pi*sqrt(m*r*R)))*( (1 - m/2)*K - E);

end