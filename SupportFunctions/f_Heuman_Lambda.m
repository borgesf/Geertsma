%% Auxiliaty Function for Heuman-Lambda Function

% Reference 1: https://math.stackexchange.com/questions/628216/heuman-lambda-function-in-matlab
% Reference 2: https://math.stackexchange.com/questions/108575/what-does-mean/108659#108659

function [L] = f_Heuman_Lambda(Beta, m)

m_dash = (1-m);

[K,E] = ellipke(m);

[incF, incE, ~] = elliptic12(Beta,m_dash);

L = 2/pi * (E*incF + K*incE - K*incF );


end