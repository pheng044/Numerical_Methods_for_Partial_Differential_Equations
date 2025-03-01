% -------------------------------------------------------------------
% Patrick Heng
% 23 Feb. 2025
% Test script to run the example from our homework. Essentially just
% calling all of the previously developed functions.
% -------------------------------------------------------------------
close all; clear all; clc;

% Colormaps from '200 Colormaps' on the MATLAB file exchange
% SWITCH THIS TO A DEFAULT COLORMAP ('parula') IF YOU DO NOT HAVE THE 
% FILE INSTALLED
cmap=slanCM(1);

% -------------------------------------------------------------------
% ----- INPUTS -----
tic
% Forcing function as a function of the physical meshgrid variables,
% 0 for the Laplace equation
F = @(X,Y) 0*X;

% Define mesh shape
h_b = 3;
h_c = 2;

% Number of nodes in the computational domain, n -> x, m -> y
n = 512;
m = 512;
% -------------------------------------------------------------------
% ----- SOLVER -----

% Total number of nodes
nodes = n*m;

% Generate computational mesh in xi, eta
xi = linspace(0,1,n);
eta = linspace(0,1,m);
[Xi,Eta] = meshgrid(xi,eta);

% Finite differences
dxi = xi(2) - xi(1);
deta = eta(2) - eta(1);

% Get bilinear map coefficients, generate physical mesh (X,Y)
[X,Y,coeffs_x,coeffs_y] = bilinear_map(Xi,Eta,[-4,4,4,-4],[0,0,h_b,h_c]);

% Compute derivatives of the coordinate transformation
[x_xi,x_eta,y_xi,y_eta,x_xi_eta,...
            y_xi_eta,J,a,b,c,d,e,alpha,beta] = ...
            bilinear_map_derivatives(Xi,Eta,coeffs_x,coeffs_y);

% Generate 9 point stencil for the main diagonals of the A matrix 
% (discrete Laplacian)
diag = stencil_coefficients(J,a,b,c,d,e,dxi,deta,nodes);
diag_idx = [-n-1, -n, -n+1, -1, 0, 1, n-1, n , n+1];

% Compute the A matrix over the internal nodes. Initially generate the A 
% matrix with one extra column to get the correct indexing for the 
% stencil diagonals (see spdiags documentation)
A = spdiags(diag,diag_idx,nodes,nodes+1);

% Remove the extra column
A = A(:,1:nodes);

% Vectorize transformation derivatives
x_xi = reshape(x_xi',nodes,1);
x_eta = reshape(x_eta',nodes,1);
y_xi = reshape(y_xi',nodes,1);
y_eta = reshape(y_eta',nodes,1);
J = reshape(J',nodes,1);

% Impose boundary conditions on the A matrix, generate the corresponding
% forcing vector, f
f = reshape(F(X,Y)',nodes,1);

[A,f] = impose_BCs(A,f,x_xi,x_eta,y_xi,y_eta,J,dxi,deta,X);

% Linear solve for nodes
u = A\f;

% Generate the electric field by taking the negative gradient 
% of the potential
[Ex,Ey] = grad_field(u,x_xi,x_eta,y_xi,y_eta,J,dxi,deta,n,m);
Ex = -Ex; Ey = -Ey;

% Magnitude of the E field
E = sqrt(Ex.^2+Ey.^2);

% Reshape u to match meshgrid variables
u = flip(reshape(u,n,m)');

toc
% -------------------------------------------------------------------
% ----- PLOTTING -----
%
% Sparsity plot of A
figure
spy(A)

% Potential field solution in the physical domain
figure
surf(X,Y,u,linestyle='none')
colormap(cmap)
box on; grid on

xlabel('$x$', interpreter='latex'); ylabel('$y$', interpreter='latex'); 
zlabel('$\phi$', interpreter='latex'); 

fontsize(16,'points'); fontname('Serif');

% Contour plot of the potential field in the physical domain
figure
contourf(X,Y,u,20,FaceAlpha=0.95)
colormap(cmap)

cb = colorbar();
xlabel('$x$', interpreter='latex'); ylabel('$y$', interpreter='latex'); 
ylabel(cb,'$\phi$',interpreter='latex')
clim([-1,1])
box on

fontsize(16,'points'); fontname('Serif');

% Surface plot of the magnitude of the electric field
figure
colormap(cmap)
surf(X,Y,E,linestyle ='none')

cb = colorbar();
xlabel('$x$', interpreter='latex'); ylabel('$y$', interpreter='latex'); 
ylabel(cb,'$|E|$',interpreter='latex')
box on; grid on

fontsize(16,'points'); fontname('Serif');

% Contour plot of potential field
figure
hold on 
sc = 3;        % Arrow density for quiver plot, higher sc = less arrows
colormap(cmap)

% Plot normalized E vector field with set arrow density
quiver(X(1:sc:end,1:sc:end),Y(1:sc:end,1:sc:end), ...
        Ex(1:sc:end,1:sc:end)./E(1:sc:end,1:sc:end), ...
        Ey(1:sc:end,1:sc:end)./E(1:sc:end,1:sc:end),color='k', ...
        linewidth=0.3)

% Plot 30 equipotential lines 
levels = linspace(min(u,[],'all'),max(u,[],'all'),30);
contour(X,Y,u,levels,linewidth=0.8)
xlim([-4,4]); ylim([0,max([h_b,h_c])]); 
xlabel('$x$', interpreter='latex'); ylabel('$y$', interpreter='latex');

cb = colorbar;
ylabel(cb,'$\phi$',interpreter='latex')

box on; grid on
fontsize(16,'points'); fontname('Serif');

% Surface plot projected on lines of constant x
figure
colormap(cmap);
levels = 32;    % Number of slices to plot
x = X(1,:);

for i = 1:levels
    k = ceil(i*size(x,2)/levels);
    plot(x,u(k,:),color=flip(cmap(floor(256*i/levels),:),1),linewidth=2)
    hold on
end

xlabel('$x$',interpreter='latex'); ylabel('$\phi$',interpreter='latex')

cb = colorbar;
clim([0,max([h_b,h_c])])
ylabel(cb,'$y$',interpreter='latex')

box on; grid on
fontsize(16,'points'); fontname('Serif');

% Surface plot projected on lines of constant y
figure
levels = 32;    % Number of slices to plot
colormap(cmap);

for i = 1:levels
    k = ceil(i*size(Y(:,1),1)/levels);
    y = Y(:,k);
    plot(y,u(:,k),color=cmap(floor(256*i/levels),:),linewidth=2)
    hold on
end

xlabel('$y$',interpreter='latex'); ylabel('$\phi$',interpreter='latex')

cb = colorbar;
clim([-4,4])
ylabel(cb,'$x$',interpreter='latex')

box on; grid on
fontsize(16,'points'); fontname('Serif');
%}


% -------------------------------------------------------------------
% IF YOU DESIRE TO SAVE THE SOLUTION: UNCOMMENT, RENAME
%save('Refined_mesh_1024_1024.mat','X','Y','Xi','Eta','u','Ex','Ey','E')
