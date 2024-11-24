% This script showcases how the Geertsma function works for displacement
% and stress analysis in a reservoir, along with visualization techniques.
%
% Reference:
%	   Fjær, E., R. M. Holt, A. Raaen, R. Risnes, and P. Horsrud,
%        2008, Petroleum related rock mechanics: Elsevier, 53.
%
% Author: Filipe Borges (filipe.borges.7@gmail.com)

% Clear environment and prepare workspace
clear;
clc;

%% Reservoir Parameters
D = 500; % Reservoir Depth [m]
R = 500; % Reservoir Radius [m]
h = 10; % Reservoir Thickness [m]
deltaP = -10e6; % Change in Pore Pressure [Pa]
numLayers = 1; % Number of layers in the reservoir

%% Medium Properties
E = 2e9; % Young's Modulus [Pa]
nu = 0.25; % Poisson's Ratio [unitless]
kMineral = 37e9; % Bulk Modulus of Mineral [Pa]

%% Output Coordinates
zCoords = linspace(eps, 1000, 101); % Vertical positions [m]
rCoords = linspace(eps, 1500, 101); % Radial positions [m]

%% Function Evaluation
% Use either symbolic or non-symbolic version based on availability

[uz, ur, sigmaZ, sigmaR] = Geertsma_No_ToolBox(D, R, h, deltaP, E, nu, kMineral, zCoords, rCoords, numLayers);
% [uz, ur, sigmaZ, sigmaR] = Geertsma_Symbolic_Toolbox(D, R, h, deltaP, E, nu, kMineral, zCoords, rCoords, numLayers);

%% Symmetry and Grid Preparation
rSymmetric = [-fliplr(rCoords), rCoords];
uzSymmetric = [fliplr(uz), uz];
urSymmetric = [fliplr(ur), ur];
sigmaZSymmetric = -[fliplr(sigmaZ), sigmaZ];
sigmaRSymmetric = -[fliplr(sigmaR), sigmaR];

[rGrid, zGrid] = meshgrid(rSymmetric, zCoords);

%% Interpolation for Visualization
[rPlot, zPlot] = meshgrid(linspace(rSymmetric(1), rSymmetric(end), 500), linspace(zCoords(1), zCoords(end), 500));

uzInterp = interpolateField(uzSymmetric, rGrid, zGrid, rPlot, zPlot);
urInterp = interpolateField(urSymmetric, rGrid, zGrid, rPlot, zPlot);
sigmaZInterp = interpolateField(sigmaZSymmetric, rGrid, zGrid, rPlot, zPlot);
sigmaRInterp = interpolateField(sigmaRSymmetric, rGrid, zGrid, rPlot, zPlot);

%% Visualization
%visualizeResults(rPlot, zPlot, uzInterp, urInterp, sigmaZInterp, sigmaRInterp, R, D, h);
visualizeResultsQuiver(rPlot, zPlot, uzInterp, urInterp, R, D, h); % Call quiver visualization

%% Helper Functions

function result = interpolateField(field, rGrid, zGrid, rPlot, zPlot)
% INTERPOLATEFIELD Interpolates field data for visualization grids.
%
% This function interpolates a 2D field dataset from a given grid to a
% specified plotting grid. It handles missing or infinite values by
% assigning them as NaN and uses interpolation methods accordingly.
%
% Inputs:
%   field   - (n, m) double array of field data to interpolate
%   rGrid   - (n, m) double array, radial grid coordinates for field
%   zGrid   - (n, m) double array, vertical grid coordinates for field
%   rPlot   - (p, q) double array, radial coordinates for plotting grid
%   zPlot   - (p, q) double array, vertical coordinates for plotting grid
%
% Outputs:
%   result  - (p, q) double array, interpolated field data
%
% Example:
%   result = interpolateField(fieldData, rGrid, zGrid, rPlot, zPlot);

% Validate inputs using arguments block
arguments
    field(:, :) double{mustBeNumeric, mustBeReal}
    rGrid(:, :) double{mustBeNumeric, mustBeReal, mustBeEqualSize(field, rGrid)}
    zGrid(:, :) double{mustBeNumeric, mustBeReal, mustBeEqualSize(field, zGrid)}
    rPlot(:, :) double{mustBeNumeric, mustBeReal}
    zPlot(:, :) double{mustBeNumeric, mustBeReal, mustBeEqualSize(rPlot, zPlot)}
end

% Replace infinities with NaN for stability
field(isinf(field)) = NaN;

% Check for and handle NaN values
if any(isnan(field), 'all')
    validIndices = ~isnan(field);
    result = griddata(rGrid(validIndices), zGrid(validIndices), field(validIndices), rPlot, zPlot, 'natural');
else
    result = interp2(rGrid, zGrid, field, rPlot, zPlot, 'spline');
end
end

function mustBeEqualSize(A, B)
% MUSTBEEQUALSIZE Inline validation function to ensure matrices have the same size.
if ~isequal(size(A), size(B))
    error('Input matrices must have the same size.');
end
end


function visualizeResults(rPlot, zPlot, uz, ur, sigmaZ, sigmaR, R, D, h)
% VISUALIZERESULTS Visualizes displacement and stress fields in a reservoir.
%
% This function generates a 2x2 subplot showing:
%   1. Vertical displacement field
%   2. Radial displacement field
%   3. Vertical stress change
%   4. Radial stress change
%
% Inputs:
%   rPlot   - (n, m) double array of radial coordinates for plotting [m]
%   zPlot   - (n, m) double array of vertical coordinates for plotting [m]
%   uz      - (n, m) double array of vertical displacements [cm]
%   ur      - (n, m) double array of radial displacements [cm]
%   sigmaZ  - (n, m) double array of vertical stress changes [Pa]
%   sigmaR  - (n, m) double array of radial stress changes [Pa]
%   R       - (1, 1) double, reservoir radius [m]
%   D       - (1, 1) double, reservoir depth [m]
%   h       - (1, 1) double, reservoir thickness [m]
%
% Outputs:
%   None. Displays the generated plots.

% Create figure with specified position
figure('Position', [400, 100, 1400, 1000]);

% Plot vertical displacement
subplot(2, 2, 1);
plotField(rPlot, zPlot, -100*uz, 'Vertical Displacement (cm)', [], R, D, h);
title('Vertical Displacement', 'FontSize', 14, 'Interpreter', 'latex');

% Plot radial displacement
subplot(2, 2, 2);
plotField(rPlot, zPlot, -100*ur, 'Radial Displacement (cm)', [], R, D, h);
title('Radial Displacement', 'FontSize', 14, 'Interpreter', 'latex');

% Plot vertical stress change
subplot(2, 2, 3);
plotField(rPlot, zPlot, -sigmaZ/1e6, 'Vertical Stress Change ($\Delta\sigma_V$) [MPa]', [], R, D, h);
title('Vertical Stress Change', 'FontSize', 14, 'Interpreter', 'latex');

% Plot radial stress change
subplot(2, 2, 4);
plotField(rPlot, zPlot, -sigmaR/1e6, 'Radial Stress Change ($\Delta\sigma_r$) [MPa]', [], R, D, h);
title('Radial Stress Change', 'FontSize', 14, 'Interpreter', 'latex');

%exportgraphics(gcf, 'Geertsma_2.png', 'Resolution', 400);
end

function plotField(x, y, field, colorbarLabel, cLimits, R, D, h)
% PLOTFIELD Helper function to plot individual fields.
%
% This function creates a filled contour plot for the provided field data
% and adds a colorbar, axis labels, and a reservoir rectangle overlay.
%
% Inputs:
%   x            - (n, m) double array of horizontal coordinates [m]
%   y            - (n, m) double array of vertical coordinates [m]
%   field        - (n, m) double array of field values to plot
%   colorbarLabel- (1, 1) string, label for the colorbar
%   cLimits      - (1, 2) double array specifying color scale limits (optional)
%                  If empty or not provided, symmetric scaling around zero
%                  is applied based on the maximum absolute value.
%   R            - (1, 1) double, reservoir radius [m]
%   D            - (1, 1) double, reservoir depth [m]
%   h            - (1, 1) double, reservoir thickness [m]
%
% Outputs:
%   None. Displays the generated plot.

% Default symmetric color scale around zero if cLimits not provided
if isempty(cLimits)
    maxAbs = max(abs(field(:)), [], 'omitnan');
    cLimits = [-maxAbs, maxAbs];
end

% Create the contour plot
contourf(x, y, field, 25, 'LineColor', 'k');
hcb = colorbar; % Create the colorbar
ylabel(hcb, colorbarLabel, 'FontSize', 15, 'Interpreter', 'latex'); % Add label to the colorbar
cmocean('balance')
%colormap jet; % Not recomended; please use a proper, divergent
%colormap

% Apply color limits
clim(cLimits);

% Plot reservoir rectangle overlay
hold on;
set(gca, 'YDir', 'reverse', 'FontSize', 12, 'Box', 'on');
rectangle('Position', [-R, (D - h / 2), 2 * R, h], 'FaceColor', 'k'); % Reservoir rectangle

% Add axis labels
xlabel('Horizontal Coordinate (m)', 'Interpreter', 'latex');
ylabel('Depth (m)', 'Interpreter', 'latex'); % Subsea true vertical depth
ylim([0 max(y(:))])

% Release hold
hold off;
end

function visualizeResultsQuiver(rPlot, zPlot, uz, ur, R, D, h)
% VISUALIZERESULTSQUIVER Visualizes displacement and stress fields using quiver plots.
%
% This function generates 2x1 quiver plots for:
%   1. Displacement field (horizontal and vertical displacements)
%   2. Stress change field (horizontal and vertical stress components)
%
% Inputs:
%   rPlot   - (n, m) double array of radial coordinates for plotting [m]
%   zPlot   - (n, m) double array of vertical coordinates for plotting [m]
%   uz      - (n, m) double array of vertical displacements [cm]
%   ur      - (n, m) double array of radial displacements [cm]
%   R       - (1, 1) double, reservoir radius [m]
%   D       - (1, 1) double, reservoir depth [m]
%   h       - (1, 1) double, reservoir thickness [m]

    % Subsampling for quiver plot readability
    decimationFactor = 50; % Reduce density of arrows
    rSub = rPlot(1:decimationFactor:end, 1:decimationFactor:end);
    zSub = zPlot(1:decimationFactor:end, 1:decimationFactor:end);
    urSub = ur(1:decimationFactor:end, 1:decimationFactor:end);
    uzSub = uz(1:decimationFactor:end, 1:decimationFactor:end);
   

    % Adjust displacements and stresses to zero inside the reservoir
    insideReservoir = (zSub >= D - h / 2) & (zSub <= D + h / 2);
    uzSub(insideReservoir) = 0;
    urSub(insideReservoir) = 0;

    % Adjust ur for horizontal displacements
    urSub(rSub < 0) = -urSub(rSub < 0);

    % Create figure
    figure('Position', [400, 100, 700, 500]);

    % Plot displacement quiver plot
    quiver(rSub, zSub, -urSub, -uzSub, 'r', 'LineWidth', 1.5, 'AutoScaleFactor', 0.3);
    hold on;
    grid on;
    box on;
    rectangle('Position', [-R, D - h / 2, 2 * R, h], 'FaceColor', [0.5, 0.5, 0.5, 0.5]); % Reservoir rectangle
    set(gca, 'YDir', 'reverse');
    xlabel('Horizontal Coordinate (m)');
    ylabel('Depth (m)');
    title('Displacement Field', 'Interpreter', 'latex', 'FontSize', 14);
    xlim(max(rPlot(:))*[-1 1])
    ylim([0 max(zPlot(:))])
    exportgraphics(gcf, 'Geertsma_quiver.png', 'Resolution', 400); 
    hold off
 
end
