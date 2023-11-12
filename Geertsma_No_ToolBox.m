function varargout = Geertsma_No_ToolBox(D,R,h,Delta_p,E,nu,K_mineral,Z,r,N_Layers)

% Function for calculation fo stress and displacement based on Geertsma's
% model. The equations implemented are the ones presented in the book
% Petroleum Related Rock Mechanics, Appendix D5.

% Syntax:
% [u_Z,u_r,Sigma_Z,Sigma_r,Sigma_t, Tau_rz] = Geertsma_No_ToolBox(D,R,h,Delta_p,...
%                                E,nu,K_mineral,Z,r,N_Layers)
%
% Inputs:
% - Reservoir Geometry and Pressure
%   D          - Reservoir Depth (meters) - center of disc-shaped reservoir
%   R          - Reservoir Radius (meters)
%   h          - Reservoir Thickness (meters)
%   Delta_Pp   - Change in Pore Pressure (Pascal)
%   N_Layers   - (Optional) Allows for calculation inside reservoir by dividing
%                the thickness 'h' in several layers. Default is N_Layers=1.
%
% - Medium Parameters
%   E          - Young's Modulus (Pa)
%   nu         - Poisson's Ratio (dimensionless)
%   K_mineral  - Mineral's Bulk Modulus (Pascal)
%
% - Coordinates for output
%   Z          - Array of depth values (meters). Strictly Positive.
%   r          - Array of radial distance values (meters). Strictly Positive.
%
% Outputs:
%   u_Z        - Vertical Displacement (meters)
%   u_r        - Radial Displacement (meters)
%   Sigma_Z    - Change in vertical stress (Pascal)
%   Sigma_r    - Change in radial stress (Pascal)
%   Sigma_t    - Change in tangential stress (Pascal)
%   Tau_rz     - Shear Stress (Pascal) 

% Author: Filipe Borges
%	First Version: 2018 
%	Last Update: 29/10/2018
%
%	References:
%	   Fjær, E., R. M. Holt, A. Raaen, R. Risnes, and P. Horsrud,
%        2008, Petroleum related rock mechanics: Elsevier, 53.
%
%     Geertsma, J., 1973, A basic theory of subsidence due to
%       reservoir compaction: the homogeneous case: Verhandelingen
%       Kon. Ned. Geol. Mijnbouwk. Gen, 28, 43-62.

%	Acknowledgments: Rodrigo Link, Filipe Leão

%	Development history:
%	   [29/10/2018]: Some updates in the modeling code for optimization.
% 
% Copyright (c) 2018, Filipe Borges
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the <organization> nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER  BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

arguments
    D (1,1) double {mustBePositive}
    R (1,1) double {mustBePositive}
    h (1,1) double {mustBePositive}
    Delta_p (1,1) double
    E (1,1) double {mustBePositive}
    nu (1,1) double {mustBeNonnegative, mustBeLessThan(nu,0.5)}
    K_mineral (1,1) double {mustBePositive}
    Z (:,1) double {mustBePositive}
    r (:,1) double {mustBePositive}
    N_Layers (1,1) double = 1
end

%% 1 - Checking consistency of input parameters

% Default values
if nargin < 10
    N_Layers = 1;
end

%% 2 - Deriving some ellastic parameters for modeling

G_fr = E/(2*(1+nu));                    % Frame Shear Modulus
K_fr = (2/3)*G_fr*(1+nu)/(1-2*nu);      % Frame Bulk Modulus 
Cm = (1+nu)*(1-2*nu)/(E*(1-nu));        % Uniaxial Compressibility
alpha = 1 - K_fr/K_mineral;             % Biot Coefficient

%% 3 - Modeling
% Here, some funtions are pre-called to accelerate the process. All 
% intermediary functions are defined in the end of the script.

u_Z = zeros(length(Z),length(r));          % Vertical Displacement
u_r = zeros(length(Z),length(r));          % Radial Displacement
Sigma_Z = zeros(length(Z),length(r));      % Change in vertical Stress
Sigma_r = zeros(length(Z),length(r));      % Change in radial Stress
Sigma_t = zeros(length(Z),length(r));      % Change in tangential Stress
Tau_rz = zeros(length(Z),length(r));       % Change in shear Stress

h_Layer = h/N_Layers;                      % In case of multiple layers

K_Disp = alpha*Cm*R*h_Layer*Delta_p/2;   % Constant factor for displacements
K_Stress = G_fr*Cm*R*h_Layer*Delta_p;    % Constant factor for stresses

for kk=1:N_Layers
    D_Layer = D - h/2 + (kk-1/2)*h_Layer;        
    
    for ii=1:length(Z)
        disp(['Layer ',num2str(kk),' of ',num2str(N_Layers),...
            ', step = ',num2str(ii),' of ',num2str(length(Z)),'.'])  
        
        for jj=1:length(r)       
            I_4_ZpD =  f_Geertsma_I4(Z(ii) + D_Layer,r(jj),R);               
            I_4_aZmD = f_Geertsma_I4(abs(Z(ii) - D_Layer),r(jj),R);
            I_6_ZpD =  f_Geertsma_I6(Z(ii) + D_Layer,r(jj),R);

                %% Strains       
                u_Z(ii,jj) = u_Z(ii,jj) + K_Disp*(...
                            sign(Z(ii) - D_Layer)*f_Geertsma_I3(abs(Z(ii) - D_Layer),r(jj),R) - ...
                            (3 - 4*nu)*f_Geertsma_I3(Z(ii) + D_Layer,r(jj),R) - 2*Z(ii)*I_4_ZpD);

                u_r(ii,jj) = u_r(ii,jj) + K_Disp*(...
                             f_Geertsma_I1(abs(Z(ii) - D_Layer),r(jj),R) + (3-...
                             4*nu)*f_Geertsma_I1(Z(ii) + D_Layer,r(jj),R) - ...
                             2*Z(ii)*f_Geertsma_I2(Z(ii) + D_Layer,r(jj),R));         


                %% Stresses 

                Sigma_Z(ii,jj) = Sigma_Z(ii,jj) + K_Stress*(...
                                -I_4_aZmD + I_4_ZpD + 2*Z(ii)*I_6_ZpD);

                Sigma_r(ii,jj) = Sigma_r(ii,jj) + K_Stress*(...
                            I_4_aZmD + 3*I_4_ZpD - 2*Z(ii)*I_6_ZpD - ...
                             (1/r(jj))*u_r(ii,jj)/K_Disp);

                Sigma_t(ii,jj) = Sigma_t(ii,jj) + K_Stress*(...
                                4*nu*f_Geertsma_I4(Z(ii) + D_Layer,r(jj),R) + ...
                                (1/r(jj))*u_r(ii,jj)/K_Disp);         
                            
                Tau_rz = -K_Stress*(-sign(Z(ii) - D_Layer)*f_Geertsma_I2(abs(Z(ii) - ...
                        D_Layer),r(jj),R) - f_Geertsma_I2(Z(ii) + D_Layer,r(jj),R) + ...
                        Z(ii)*f_Geertsma_I7(Z(ii) + D_Layer,r(jj),R));             
        end
    end
end

%% outputs

% Format outputs
varargout = cell(1, nargout);
outputs = {u_Z, u_r, Sigma_Z, Sigma_r, Sigma_t, Tau_rz};

for k = 1:nargout
    varargout{k} = outputs{k};
end

end

