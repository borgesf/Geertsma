%% Auxiliaty Function or I_3
function [I3] = f_Geertsma_I3(q,r,R)

m = 4*R*r/(q^2 + (r+R)^2);
Beta = asin(q/sqrt(q^2 + (R-r)^2));

I3 = -(q*sqrt(m)*ellipke(m))/( 2*pi*sqrt(r*R)*R ) + ...
     0.5*( sign(r-R) -  sign(R-r) )*...
     f_Heuman_Lambda(Beta, m)/(2*R) + 0.5*(1/R)*(1 + sign(R-r));

end