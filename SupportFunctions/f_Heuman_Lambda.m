%% Auxiliaty Function for Heuman-Lambda Function

function [L] = f_Heuman_Lambda(Beta,m)
[K, E] = ellipke(m); 
[F, EInc, ~] = elliptic12(Beta,m);

L = (2/pi)*(E*F + K*EInc - K*F);

end