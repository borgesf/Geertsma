%% Simple example of how Geertsma function works, and some suggestions 
% for displaying the results.
% 
% Author: Filipe Borges (filipe.borges.7@gmail.com)
% Date: 29/10/2018

    clear;
    clc;
    
% Reservoir geometry and Pressure change

    D = 400;                    % Reservoir Depth (meters)
    R = 300;                    % Reservoir Radius (meters)
    h = 100;                    % Reservoir Thickness (meters)
    Delta_p = -10*10^6;         % Change in Pore Pressure (Pascal)
    N_Layers = 10;               %

% Medium Parameters

    E = 2*10^9;                 % Young's Modulus (Pa)
    Nu = 0.25;                  % Nu - Poisson's Ratio (dimensionless)
    K_mineral = 37*10^9;        % K_mineral - Mineral's Bulk Modulus (Pascal) (Values for quartz)

% Coordinates for output
    Z = linspace(eps,1000,101);
    r = linspace(eps,1000,101);
        
% Running function (selec the correct one, whether you have or not the
% symbolic toolbox package)
    [uz,ur,sigmaz,sigmar] = Geertsma_No_ToolBox(D,R,h,Delta_p,E,Nu,K_mineral,Z,r,N_Layers);
    %[uz,ur,sigmaz,sigmar] = Geertsma_Symbolic_Toolbox(D,R,h,Delta_p,E,Nu,K_mineral,Z,r,N_Layers);
    

%% Interpolating results for display
% Tip: define "symmetric" arrays to allow for better visualization

r_symmetric = [-fliplr(r),r];

u_z_symmetric = [fliplr(uz),uz];
u_r_symmetric = [fliplr(ur),ur];
Sigma_z_symmetric = -[fliplr(sigmaz),sigmaz];   % Minus sign for having positive == compression
Sigma_r_symmetric = -[fliplr(sigmar),sigmar];   % Minus sign for having positive == compression

[RGRID,ZGRID] = meshgrid(r_symmetric,Z);        % Model Grid

rplot = linspace(r_symmetric(1),r_symmetric(end),500);
zplot = linspace(Z(1),Z(end),500);

[RPLOT,ZPLOT] = meshgrid(rplot,zplot);          % Display Grid

% Identifying possible NaN or infinity values, can cause numerical instability
u_z_symmetric(u_z_symmetric==inf | u_z_symmetric==-inf) = NaN;
u_r_symmetric(u_r_symmetric==inf | u_r_symmetric==-inf) = NaN;
Sigma_z_symmetric(Sigma_z_symmetric==inf | Sigma_z_symmetric==-inf) = NaN;
Sigma_r_symmetric(Sigma_z_symmetric==inf | Sigma_z_symmetric==-inf) = NaN;


if sum(sum(isnan(u_z_symmetric)))>0
   indx = ~(isnan(u_z_symmetric));  
   Interpolated_uz = griddata(RGRID(indx),ZGRID(indx),u_z_symmetric(indx),RPLOT,ZPLOT,'natural');
else
   Interpolated_uz = interp2(RGRID,ZGRID,u_z_symmetric,RPLOT,ZPLOT,'spline');
end

if sum(sum(isnan(u_r_symmetric)))>0
   indx = ~(isnan(u_r_symmetric));  
   Interpolated_ur = griddata(RGRID(indx),ZGRID(indx),u_r_symmetric(indx),RPLOT,ZPLOT,'natural');
else
   Interpolated_ur = interp2(RGRID,ZGRID,u_r_symmetric,RPLOT,ZPLOT,'spline');
end

if sum(sum(isnan(Sigma_z_symmetric)))>0
   indx = ~(isnan(Sigma_z_symmetric));  
   Interpolated_sz = griddata(RGRID(indx),ZGRID(indx),Sigma_z_symmetric(indx),RPLOT,ZPLOT,'natural');
else
   Interpolated_sz = interp2(RGRID,ZGRID,Sigma_z_symmetric,RPLOT,ZPLOT,'spline');
end

if sum(sum(isnan(Sigma_r_symmetric)))>0
   indx = ~(isnan(Sigma_r_symmetric));  
   Interpolated_sr = griddata(RGRID(indx),ZGRID(indx),Sigma_r_symmetric(indx),RPLOT,ZPLOT,'natural');
else
   Interpolated_sr = interp2(RGRID,ZGRID,Sigma_r_symmetric,RPLOT,ZPLOT,'spline');
end


%% Displaying Results

figure

subplot(221)
hold on
box on
set(gca,'FontSize',12)
title('Vertical Displacement')
contourf(RPLOT,ZPLOT,-100*Interpolated_uz,25)
%contourf(RPLOT,ZPLOT,-100*Interpolated_uz_grid,25)
colormap jet
hcb = colorbar;
ylabel(hcb,'Vertical Displacement (cm)','FontSize',15)
xlim([r_symmetric(1) r_symmetric(end)])
ylim([0 Z(end)])
set(gca,'YDir','reverse')
%caxis([-75 80])
xlabel('Horizontal Distance (m)')
ylabel('SSTVD (m)')
set(gca,'XTick',linspace(-1000,1000,5),'XTickLabel',abs(linspace(-1000,1000,5)));
rectangle('Position',[-R, (D - h/2),2*R,h],'FaceColor','k')
hold off

subplot(222)
hold on
box on
set(gca,'FontSize',12)
title('Radial Displacement')
contourf(RPLOT,ZPLOT,-100*Interpolated_ur,25)
colormap jet
hcb = colorbar;
ylabel(hcb,'Radial Displacement (cm)','FontSize',15)
xlim([r_symmetric(1) r_symmetric(end)])
ylim([0 Z(end)])
set(gca,'YDir','reverse')
%caxis([0 1])
xlabel('Horizontal Distance (m)')
ylabel('SSTVD (m)')
set(gca,'XTick',linspace(-1000,1000,5),'XTickLabel',abs(linspace(-1000,1000,5)));
rectangle('Position',[-R, (D - h/2),2*R,h],'FaceColor','k')
hold off


subplot(223)
hold on
box on
set(gca,'FontSize',12)
title('Vertical Stress Change')
contourf(RPLOT,ZPLOT,Interpolated_sz/10^6,25)
colormap jet
caxis([-2 2])
hcb = colorbar;
ylabel(hcb,'\Delta\sigma_V (MPa)','FontSize',15)
xlim([r_symmetric(1) r_symmetric(end)])
ylim([0 Z(end)])
set(gca,'YDir','reverse')
xlabel('Horizontal Distance (m)')
ylabel('SSTVD (m)')
set(gca,'XTick',linspace(-1000,1000,5),'XTickLabel',abs(linspace(-1000,1000,5)));
rectangle('Position',[-R, (D - h/2),2*R,h],'FaceColor','k')
hold off

subplot(224)
hold on
box on
set(gca,'FontSize',12)
%title('Radial Stress Change')
contourf(RPLOT,ZPLOT,Interpolated_sr/10^6,50)
colormap jet
%caxis([-2 2])
hcb = colorbar;
ylabel(hcb,'\Delta\sigma_r (MPa)','FontSize',15)
xlim([r_symmetric(1) r_symmetric(end)])
ylim([0 Z(end)])
set(gca,'YDir','reverse')
xlabel('Horizontal Distance (m)')
ylabel('SSTVD (m)')
set(gca,'XTick',linspace(-1000,1000,5),'XTickLabel',abs(linspace(-1000,1000,5)));
rectangle('Position',[-R, (D - h/2),2*R,h],'FaceColor','k')
hold off

%%