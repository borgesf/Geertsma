function u_Z = Geertsma_Exact(D, R, h, Delta_p, E, nu, K_mineral, Z, shouldPlot)
% GEERTSMA_EXACT Calculate the exact vertical compaction above the axis of 
% a compacting cylinder based on Geertsma's model in SI units.
%
% Inputs:
% D         - Depth to the center of the reservoir (meters, must be positive)
% R         - Radius of the cylinder representing the reservoir (meters, must be positive)
% h         - Thickness of the reservoir (meters, must be positive)
% Delta_p   - Change in pore pressure (Pascals)
% E         - Young's modulus of the reservoir rock (Pascals, must be positive)
% nu        - Poisson's ratio of the reservoir rock (dimensionless, 0 <= nu < 0.5)
% K_mineral - Bulk modulus of the minerals in the rock (Pascals, must be positive)
% Z         - Vector of depths at which to calculate displacement (meters, must be positive)
% shouldPlot - Flag to plot the results (logical, true/false, default true)
%
% Output:
% u_Z       - Calculated vertical displacement at depths Z (meters)
%
% The function assumes all parameters are in SI units.
% Author: Filipe Borges (filipe.borges.7@gmail.com)


arguments
    D (1,1) double {mustBePositive} = 2000
    R (1,1) double {mustBePositive} = 100
    h (1,1) double {mustBePositive} = 100
    Delta_p (1,1) double = 5e6
    E (1,1) double {mustBePositive} = 20e9
    nu (1,1) double {mustBeNonnegative, mustBeLessThan(nu,0.5)} = 0.3
    K_mineral (1,1) double {mustBePositive} = 20e9
    Z (:,1) double {mustBePositive} = eps:2500
    shouldPlot logical = true
end

% Derive some elastic parameters for modeling
G_fr = E/(2*(1+nu));                    % Frame Shear Modulus
K_fr = (2/3)*G_fr*(1+nu)/(1-2*nu);      % Frame Bulk Modulus 
Cm = (1+nu)*(1-2*nu)/(E*(1-nu));        % Uniaxial Compressibility
alpha = 1 - K_fr/K_mineral;             % Biot Coefficient

% Calculate the vertical displacement
u_Z = -0.5 * Cm * h * alpha * Delta_p *  ( 3 - 4*nu + sign(D-Z) - ...
       (D-Z)./sqrt(R^2 + (D-Z).^2)  - (D+Z).*(3-4*nu) ./ sqrt(R^2 + (D+Z).^2) + ...
       (2 * R^2 * Z) ./ ((R^2 + (D+Z).^2).^(3/2)));


% Conditional plotting based on shouldPlot flag
if shouldPlot
    figure
    hold on
    box on
    set(gcf, 'Position', [100, 200, 1200, 800]);
    grid on
    set(gca,'FontSize',14)

    % Define the reservoir boundaries
    reservoirTop = D - h/2;
    reservoirBottom = D + h/2;

    % Find indices for Z that are above and below the reservoir boundaries
    aboveReservoir = Z < reservoirTop;
    belowReservoir = Z > reservoirBottom;

    % Plot the displacement for these indices in two separate lines
    plot(100*u_Z(aboveReservoir), Z(aboveReservoir), 'r', 'LineWidth', 1.8,'DisplayName','Vertical displacement - Geertsma')
    plot(100*u_Z(belowReservoir), Z(belowReservoir), 'r', 'LineWidth', 1.8,'HandleVisibility','off')

    xlabel('Vertical displacement (cm)')
    ylabel('Depth (m)')

    XLIM = xlim;
    Y_Patch = [reservoirBottom, reservoirBottom, reservoirTop, reservoirTop]; 
    X_Patch = [XLIM(1), XLIM(2), XLIM(2), XLIM(1)];
    patch(X_Patch, Y_Patch, 'k', 'FaceAlpha', 0.25, 'EdgeColor', 'none','DisplayName','Reservoir Area');

    set(gca,'YDir','reverse')

    % Get the Y limits of the left y-axis
    leftAxisLimits = get(gca, 'YLim');

    yyaxis right
    ax = gca;
    ax.YColor = 'k';
    set(gca,'YDir','reverse')

    % Match the limits and ticks of the right y-axis with the left y-axis
    set(ax, 'YLim', leftAxisLimits)

    legend('Location','Northeast');
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);

    hold off
end
