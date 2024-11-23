function varargout = Geertsma_No_ToolBox(D, R, h, Delta_p, E, nu, K_mineral, Z, r, N_Layers)
% GEERTSMA_NO_TOOLBOX Computes stress and displacement based on Geertsma's model.
%
% This function calculates stress and displacement fields using Geertsma's
% model as described in "Petroleum Related Rock Mechanics" (Appendix D5).
%
% Syntax:
% [u_Z, u_r, Sigma_Z, Sigma_r, Sigma_t, Tau_rz] = Geertsma_No_ToolBox(...
%     D, R, h, Delta_p, E, nu, K_mineral, Z, r, N_Layers)
%
% Inputs:
%   D         - (1,1) double, reservoir depth (meters), center of the disc-shaped reservoir.
%   R         - (1,1) double, reservoir radius (meters).
%   h         - (1,1) double, reservoir thickness (meters).
%   Delta_p   - (1,1) double, change in pore pressure (Pascal).
%   N_Layers  - (1,1) double, number of layers for internal reservoir division (default: 1).
%   E         - (1,1) double, Young's modulus (Pa).
%   nu        - (1,1) double, Poisson's ratio (dimensionless, must be < 0.5).
%   K_mineral - (1,1) double, mineral bulk modulus (Pa).
%   Z         - (n,1) double array, depth values (meters), strictly positive.
%   r         - (m,1) double array, radial distance values (meters), strictly positive.
%
% Outputs:
%   u_Z       - (n,m) double array, vertical displacement (meters).
%   u_r       - (n,m) double array, radial displacement (meters).
%   Sigma_Z   - (n,m) double array, vertical stress change (Pascal).
%   Sigma_r   - (n,m) double array, radial stress change (Pascal).
%   Sigma_t   - (n,m) double array, tangential stress change (Pascal).
%   Tau_rz    - (n,m) double array, shear stress (Pascal).
%
% References:
%   Fjær, E., et al., 2008, "Petroleum Related Rock Mechanics": Elsevier.
%   Geertsma, J., 1973, "Basic theory of subsidence due to reservoir compaction".
%
% Acknowledgments:
%   Rodrigo Link, Filipe Leão
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

% Argument validation
arguments
    D(1, 1) double{mustBePositive}
    R(1, 1) double{mustBePositive}
    h(1, 1) double{mustBePositive}
    Delta_p(1, 1) double
    E(1, 1) double{mustBePositive}
    nu(1, 1) double{mustBeNonnegative, mustBeLessThan(nu, 0.5)}
    K_mineral(1, 1) double{mustBePositive}
    Z(:, 1) double{mustBePositive}
    r(:, 1) double{mustBePositive}
    N_Layers(1, 1) double = 1
end

%% 1. Check Consistency and Initialize
% Default layer handling
if nargin < 10
    N_Layers = 1;
end

% Calculate elastic parameters
G_fr = E / (2 * (1 + nu)); % Frame shear modulus
K_fr = (2 / 3) * G_fr * (1 + nu) / (1 - 2 * nu); % Frame bulk modulus
Cm = (1 + nu) * (1 - 2 * nu) / (E * (1 - nu)); % Uniaxial compressibility
alpha = 1 - K_fr / K_mineral; % Biot coefficient

% Reservoir compaction and displacement
Reservoir_compaction = 0.06; % Hardcoded compaction value (can be adjusted)
Displacement_resTop = -0.5 * Reservoir_compaction * ( ...
    3 - 4 * nu + 1 - 2 * D * (3 - 4 * nu) / sqrt(R^2+4*D^2) + ...
    2 * D * R^2 / (R^2 + 4 * D^2)^1.5);
Displacement_resBottom = -0.5 * Reservoir_compaction * ( ...
    3 - 4 * nu - 1 - 2 * D * (3 - 4 * nu) / sqrt(R^2+4*D^2) + ...
    2 * D * R^2 / (R^2 + 4 * D^2)^1.5);
Displacement_Surface = -2 * Reservoir_compaction * (1 - nu) * ...
    (1 - D / sqrt(D^2+R^2));

fprintf('Reservoir compaction: %.2f cm.\n', -100*Reservoir_compaction);
fprintf('Reservoir top displacement: %.2f cm.\n', 100*Displacement_resTop);
fprintf('Reservoir bottom displacement: %.2f cm.\n', 100*Displacement_resBottom);
fprintf('Surface center displacement: %.2f cm.\n', 100*Displacement_Surface);

%% 2. Initialize Output Arrays
u_Z = zeros(length(Z), length(r)); % Vertical displacement
u_r = zeros(length(Z), length(r)); % Radial displacement
Sigma_Z = zeros(length(Z), length(r)); % Vertical stress
Sigma_r = zeros(length(Z), length(r)); % Radial stress
Sigma_t = zeros(length(Z), length(r)); % Tangential stress
Tau_rz = zeros(length(Z), length(r)); % Shear stress

h_Layer = h / N_Layers; % Thickness of each layer
K_Disp = Reservoir_compaction * R / (2 * N_Layers); % Displacement factor
K_Stress = G_fr * Reservoir_compaction * R / (alpha * N_Layers); % Stress factor

%% 3. Model Stress and Displacement
for kk = 1:N_Layers
    D_Layer = D - h / 2 + (kk - 1 / 2) * h_Layer; % Layer depth
    fprintf('Processing Layer %d of %d...\n', kk, N_Layers);

    for ii = 1:length(Z)
        for jj = 1:length(r)
            % Intermediate calculations
            I_4_ZpD = f_Geertsma_I4(Z(ii)+D_Layer, r(jj), R);
            I_4_aZmD = f_Geertsma_I4(abs(Z(ii)-D_Layer), r(jj), R);
            I_6_ZpD = f_Geertsma_I6(Z(ii)+D_Layer, r(jj), R);

            % Displacements
            u_Z(ii, jj) = u_Z(ii, jj) + K_Disp * ( ...
                sign(Z(ii)-D_Layer) * f_Geertsma_I3(abs(Z(ii)-D_Layer), r(jj), R) - ...
                (3 - 4 * nu) * f_Geertsma_I3(Z(ii)+D_Layer, r(jj), R) - ...
                2 * Z(ii) * I_4_ZpD);

            u_r(ii, jj) = u_r(ii, jj) + K_Disp * ( ...
                f_Geertsma_I1(abs(Z(ii)-D_Layer), r(jj), R) + ...
                (3 - 4 * nu) * f_Geertsma_I1(Z(ii)+D_Layer, r(jj), R) - ...
                2 * Z(ii) * f_Geertsma_I2(Z(ii)+D_Layer, r(jj), R));

            % Stresses
            Sigma_Z(ii, jj) = Sigma_Z(ii, jj) + K_Stress * (-I_4_aZmD + I_4_ZpD + 2 * Z(ii) * I_6_ZpD);

            Sigma_r(ii, jj) = Sigma_r(ii, jj) + K_Stress * ( ...
                I_4_aZmD + 3 * I_4_ZpD - 2 * Z(ii) * I_6_ZpD - ...
                (1 / r(jj)) * u_r(ii, jj) / K_Disp);

            Sigma_t(ii, jj) = Sigma_t(ii, jj) + K_Stress * ( ...
                4 * nu * I_4_ZpD + (1 / r(jj)) * u_r(ii, jj) / K_Disp);

            Tau_rz(ii, jj) = -K_Stress * ( ...
                -sign(Z(ii)-D_Layer) * f_Geertsma_I2(abs(Z(ii)-D_Layer), r(jj), R) - ...
                f_Geertsma_I2(Z(ii)+D_Layer, r(jj), R) + ...
                Z(ii) * f_Geertsma_I7(Z(ii)+D_Layer, r(jj), R));
        end
    end
end

%% 4. Prepare Outputs
outputs = {u_Z, u_r, Sigma_Z, Sigma_r, Sigma_t, Tau_rz};
varargout = outputs(1:nargout);
end
