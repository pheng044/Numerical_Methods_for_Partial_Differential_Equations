% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 3 Mar 2025
% Sript to plot ketene flow and power consumption of PFR for varying
% cA0 and Ta
% -------------------------------------------------------------------
close all; clear variables; clc

% Colormaps from '200 Colormaps' on the MATLAB file exchange,
% SWITCH THIS TO A DEFAULT COLORMAP ('parula') IF YOU DO NOT HAVE THE 
% FILE INSTALLED
cmap=slanCM(100);

% Load optimization surface data
load('PFR_optimization')
% Mesh to plot over
[C,T] = meshgrid(c_vec,T_vec);

% -------------------------------------------------------------------
% Surface plot of ketene flow for varying cA0 and Ta
figure
surf(C,T,F',linestyle='none')

% Pretty plot parameters
xlabel('$c_{A,0}$ (mol/m$^3$)',interpreter='latex')
ylabel('$T_a$ (K)',interpreter='latex')
zlabel('$F_B$ (mol/s)',interpreter='latex')
colormap(cmap)
cb = colorbar;
ylabel(cb,'$F_B$ (mol/s)',interpreter='latex')
axis tight
fontname('Serif'); fontsize(16,'points')

% -------------------------------------------------------------------
% Surface projection plot of ketene flow for varying Ta along slices
% of constant cA0
figure
hold on
for i = 1:2:length(c_vec)
    plot(T_vec,F(i,:),linewidth=2,color= ...
                cmap(floor(256*i/length(c_vec)),:))
end

% Pretty plot parameters
xlabel('$T_a$ (K)',interpreter='latex')
ylabel('$F_B$ (mol/s)',interpreter='latex')
colormap(cmap)
cb = colorbar;
ylabel(cb,'$c_{A,0}$ (mol/m$^3$)',interpreter='latex')
clim([min(c_vec), max(c_vec)])
grid on
box on
fontname('Serif'); fontsize(16,'points')

% -------------------------------------------------------------------
% Surface projection plot of ketene flow for varying cA0 along slices
% of constant Ta
figure
hold on
for i = 1:2:length(T_vec)
    plot(c_vec,F(:,i),linewidth=2,color= ...
                cmap(floor(256*i/length(T_vec)),:))
end

% Pretty plot parameters
xlabel('$c_{A,0}$ (mol/m$^3$)',interpreter='latex')
ylabel('$F_B$ (mol/s)',interpreter='latex')
colormap(cmap)
cb = colorbar;
ylabel(cb,'$T_a$ (K)',interpreter='latex')
clim([min(T_vec), max(T_vec)])

axis tight
grid on
box on
fontname('Serif'); fontsize(16,'points')

% -------------------------------------------------------------------
% Surface plot of PFR power consumption for varying cA0 and Ta
figure
surf(C,T,W',linestyle='none')

% Pretty plot parameters
xlabel('$c_{A,0}$ (mol/m$^3$)',interpreter='latex')
ylabel('$T_a$ (K)',interpreter='latex')
zlabel('$\dot{W}$ (W)',interpreter='latex')
set(gca, 'ZScale', 'log')
colormap(cmap)
cb = colorbar;
ylabel(cb,'$\dot{W}$ (W)',interpreter='latex')
set(gca,'ColorScale','log')
axis tight
fontname('Serif'); fontsize(16,'points')
