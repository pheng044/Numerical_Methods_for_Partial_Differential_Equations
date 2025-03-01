% -------------------------------------------------------------------
% Patrick Heng
% 24 Feb. 2025
% Plot diagram for electrolocation case study.
% -------------------------------------------------------------------

close all; clear all; clc;

% Normal distance and angle of top boundary
theta = atan2(2,8);
h_n = 2 + 2*tan(theta);

figure
hold on

% Plot domain
line([-4,-4],[0,2], linewidth=1,color='k')
line([4,4],[0,3], linewidth=1,color='k')
line([-4,4],[0,0], linewidth=1,color='k')
line([-4,4],[2,3], linewidth=1,color='k')
line([-4,4],[2,2], linewidth=1,color='k')
line([0,0],[0,h_n], linewidth=1,color='k')

% Source
line([-2,0],[0,0],linewidth=3,color='b')
line([0,2],[0,0],linewidth=3,color='r')

% Varical arrows
text(1.5,2.85,'$\theta$',interpreter='latex')
text(0.7,1.6,'$h_n$',interpreter='latex')
text(-2.9,1.6,'$h_c$',interpreter='latex')
text(2.7,1.6,'$h_b$',interpreter='latex')
annotation('arrow',[0.6 0.63],[0.9 0.85],colo='r')
annotation('arrow',[.55 .55],[0.15 .75])
annotation('arrow',[.17 .17],[0.15 .62])
annotation('arrow',[.86 .86],[0.15 .87])

% Bottom arrow
text(0.7,0.16,'$L$',interpreter='latex')
annotation('arrow',[0.2 0.9],[0.13 0.13])
annotation('arrow',[0.9 0.135],[0.13 0.13])


xticks([]); yticks([])
xlabel('$x$',interpreter='latex'); ylabel('$y$',interpreter='latex')
fontsize(16,'points'); fontname('Serif'); box on; grid on;
