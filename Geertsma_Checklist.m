%% This is a script to check the consistency of the auxiliary functions 
% Geertsma's model. It simply runs all the functions with some pre-set parameters
% and compare then to the values provided by Fjær et al., appendex D5.
% 
% Reference:
%	   Fjær, E., R. M. Holt, A. Raaen, R. Risnes, and P. Horsrud,
%        2008, Petroleum related rock mechanics: Elsevier, 53.
%
% Author: Filipe Borges (filipe.borges@ntnu.no)
% Date: 29/10/2018

    clear;
    clc;
    
%% Calling function

r_1 = 0.2;
q_1 = 0.4;
R_1 = 1.2;

r_2 = 0.4;
q_2 = 0.4;
R_2 = 1.0;

FUNCTIONS = [f_Geertsma_I1(q_1,r_1,R_1), f_Geertsma_I2(q_1,r_1,R_1), f_Geertsma_I3(q_1,r_1,R_1),...
            f_Geertsma_I4(q_1,r_1,R_1), f_Geertsma_I6(q_1,r_1,R_1), f_Geertsma_I7(q_1,r_1,R_1);
            f_Geertsma_I1(q_2,r_2,R_2), f_Geertsma_I2(q_2,r_2,R_2), f_Geertsma_I3(q_2,r_2,R_2),...
            f_Geertsma_I4(q_2,r_2,R_2), f_Geertsma_I6(q_2,r_2,R_2), f_Geertsma_I7(q_2,r_2,R_2)];
    
%% Table of values from book

CHKLIST = [0.0595689 0.0461035 0.565282 0.598431 0.477736 -0.0544387;...
            0.162188 0.197103 0.592418 0.818683 1.15435 -0.0445251];
        
%% Comparing 
 
error = abs(FUNCTIONS - CHKLIST)./(abs(CHKLIST)); 

if sum(sum(error))<10^-4
   disp('Bro, Functions are OK o/') 
end

%%

%%