close all, clear all; clc;

load('Refined_mesh_1000_1000.mat')

% Potential field solution in physical domain
figure
surf(X,Y,u,linestyle='none')

xlabel('$x$', interpreter='latex'); ylabel('$y$', interpreter='latex'); 
zlabel('$\phi$', interpreter='latex'); 

fontsize(16,'points'); fontname('Serif');

% Vector field of the electric potential
figure
quiver(X(1:1000:end),Y(1:1000:end),Ex(1:1000:end)./E(1:1000:end),Ey(1:1000:end)./E(1:1000:end))

% Surface plot of the magnitude of the electric field
figure
surf(X,Y,E,linestyle ='none')

% Contour plot of potential field

figure
hold on 
sc = 25;
quiver(X(1:sc:end,1:sc:end),Y(1:sc:end,1:sc:end),Ex(1:sc:end,1:sc:end)./E(1:sc:end,1:sc:end),Ey(1:sc:end,1:sc:end)./E(1:sc:end,1:sc:end))
levels = linspace(min(u,[],'all'),max(u,[],'all'),25);
contour(X,Y,u,levels,linewidth=2,color='k')
