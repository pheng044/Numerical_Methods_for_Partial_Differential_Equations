% -------------------------------------------------------------------
% Patrick Heng
% 24 Feb. 2025
% Test script to PLOT data for numerical electrolocation case study.
% MUST run 'HW_1_electrolocation.m' script beforehand to generate the
% data.
% -------------------------------------------------------------------

% Colormaps from '200 Colormaps' on the MATLAB file exchange
% SWITCH THIS TO A DEFAULT COLORMAP ('parula') IF YOU DO NOT HAVE THE 
% FILE INSTALLED
cmap=slanCM(1);

load('Electrolocation_study_data.mat')

DEG = THETA*180/pi;

% Contour plot
figure 

colormap(cmap)
contourf(HN,DEG,dEy')
xlim([1,10])
xticks(1:2:10)
ylim([0,70])

xlabel('$h_n$', interpreter='latex')
ylabel('$\theta\ (^\mathrm{o})$', interpreter='latex')

cb = colorbar;
ylabel(cb,'$\Phi$',interpreter='latex')

box on; grid on
fontsize(16,'points'); fontname('Serif');

% 3D surface plot
figure 

colormap(cmap)
surf(HN,DEG,dEy',linestyle='none')
xlim([1,10])
xticks(1:2:10)
ylim([0,70])

xlabel('$h_n$', interpreter='latex')
ylabel('$\theta\ (^\mathrm{o})$',interpreter='latex')
zlabel('$\Phi$', interpreter='latex')

box on; grid on
fontsize(16,'points'); fontname('Serif');



% Surface plot projected on lines of constant x
figure
colormap(cmap);
levels = 10;    % Number of slices to plot
hn = HN(1,:);

for i = 1:levels
    k = floor(i*size(hn,2)/levels);
    plot(hn,dEy(:,k),color=cmap(floor(256*(levels-i+1)/levels),:), ...
        linewidth=2)
    hold on
end

xlim([1,10])

xlabel('$h_n$',interpreter='latex')
ylabel('$\Phi$',interpreter='latex')

cb = colorbar;
clim([0,70])
ylabel(cb,'$\theta\ (^\mathrm{o})$',interpreter='latex')

box on; grid on
fontsize(16,'points'); fontname('Serif');

% Surface plot projected on lines of constant y
figure
levels = 10;    % Number of slices to plot
colormap(cmap);
deg = DEG(:,1);

for i = 1:levels
    k = floor(i*size(deg,1)/levels);
    plot(deg,dEy(k,:),color=cmap(floor(256*(levels-i+1)/levels),:),linewidth=2)
    hold on
end

xlim([0,70])

xlabel('$\theta\ (^\mathrm{o})$',interpreter='latex')
ylabel('$\Phi$',interpreter='latex')

cb = colorbar;
ylabel(cb,'$h_n$',interpreter='latex')
clim([1,10])

box on; grid on
fontsize(16,'points'); fontname('Serif');

