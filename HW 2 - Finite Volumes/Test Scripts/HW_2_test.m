% -------------------------------------------------------------------
% 'HW_2_test.m'
% Patrick Heng
% 29 Mar. 2025
% Test script to run the finite volume solver for the convection-
% diffusion equations.
% -------------------------------------------------------------------

close all; clear variables; clc;

% Colormaps from '200 Colormaps' on the MATLAB file exchange
% SWITCH THIS TO A DEFAULT COLORMAP ('parula') IF YOU DO NOT HAVE THE 
% FILE INSTALLED
cmap=slanCM(2);

% -------------------------------------------------------------------
% ----- INPUTS -----
upwind = true;
diffusion = false;
step_IC = false;
periodic = true;

% Number of nodes (cell centers) in the computational domain, 
% n -> x, m -> y
n = 120;
m = 60;

% Simulation time
t_end = 12;

% Courant number, used to determine time step size
Co = 1;

% Initial concentration distribution
IC = @(X,Y) 0.25*exp(-50*((X-0.5).^2+(Y-0.5).^2));

% Diffusion constant
D = 0.005;

% Sink strength, 0 = off
q = 0;

% Flow velocity components in absence of sources/sinks
u0 = 2;
v0 = 0;

% Velocity components as functions of the grid points
U = @(X,Y) u0 - q*(X-2)./((X-2).^2+(Y-2).^2);
V = @(X,Y) v0 - q*(Y-2)./((X-2).^2+(Y-2).^2);

% Length of domain
Lx = 2;
Ly = 1;

% -------------------------------------------------------------------
% ----- SPATIAL DISCRETIZATION -----
% Total number of cell centers
nodes = n*m;

% Spatial discretization steps
dx = Lx/n;
dy = Ly/m;

% Create grid of cell centers
x = linspace(dx/2,Lx-dx/2,n);
y = linspace(dy/2,Ly-dy/2,m);
[X,Y] = meshgrid(x,y);

% Calculate step size based on Courant number
%dt = Co*min(dx,dy)/max(abs(u0),abs(v0));
dt = Co/(u0/dx+v0/dy);

% Determine flux function
if upwind == true
    F = upwind_flux(X,Y,U,V,dx,dy,n,m,nodes,periodic);
else
    [Fx,Fy] = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,m,nodes,periodic);
end

% Turn on/off diffusion
if diffusion == true
    G = diffusive_flux(D,dx,dy,n,m,nodes);
    if upwind == true
        A = speye(nodes) - dt*F + dt*G;
    else
        Ax = speye(nodes) - dt*Fx;
        Ay = speye(nodes) - dt*Fy + dt*G;
        A = Ay*Ax;
    end
else
    if upwind == true
        A = speye(nodes) - dt*F;
    else
        Ax = speye(nodes) - dt*Fx;
        Ay = speye(nodes) - dt*Fy;
        A = Ay*Ax;
    end
end

% Evaluate initial condition at cell centers
phi = IC(X,Y);

% If true, any values below 1% of the peak amplitude are set to 0 and
% andy values above are set to the peak amplitude
if step_IC == true
    phi(phi<0.0025) = 0;
    phi(phi>=0.0025) = 0.25;
end

phi0 = phi;

% Determine the amount of time steps to get to tf
time_steps = ceil(t_end/dt);

% Store total pollutant mass, peak value of concentration
phi_tot = zeros(time_steps,1);
phi_max = zeros(time_steps,1);

% -------------------------------------------------------------------
% ----- TIME INTEGRATION -----
% Loop through time in an explicit marching scheme
for k = 1:time_steps
    
    % Store current pollutant mass, peak amplitude
    phi_tot(k) = sum(phi,'all')*dx*dy;
    phi_max(k) = max(phi,[],'all');

    % For any of plots, remove '{'  from the '%{' below to use

    %
    colormap(cmap);
    surf(X,Y,phi,linestyle='none')
    axis equal 
    view(0,90)
    cb = colorbar;
    ylabel(cb,'$\phi$',interpreter='latex')
    xlabel('$x$',interpreter='latex')
    ylabel('$y$',interpreter='latex')
    fontsize(16,'points'); fontname('Serif')
    box on
    grid on
    pause(0.01)
    %}

    %{
    colormap(cmap);
    surf(X,Y,phi,linestyle='none')
    xlabel('$x$',interpreter='latex')
    ylabel('$y$',interpreter='latex')
    zlabel('$\phi$',interpreter='latex')
    axis equal 
    colorbar
    cb = colorbar;
    ylabel(cb,'$\phi$',interpreter='latex')
    fontsize(16,'points'); fontname('Serif')
    box on
    grid on
    pause(0.01)
    %}

    %}

    % Toggle for centerline plot animation
    %{
    plot(x,phi(round(m/2),:),linewidth=2,marker='none',color='b')
    ylim([-0.05,0.3])
    xticks(0:6); yticks(0:0.1:0.3)
    xlabel('$x$',interpreter='latex')
    ylabel('$\phi$',interpreter='latex')
    fontsize(16,'points'); fontname('Serif')
    pause(0.01)
    %}

    % Vectorize phi
    phi = reshape(phi',nodes,1);
    % Update solution
    phi = A*phi;
    % Convert phi to meshgrid for plotting
    phi = reshape(phi,n,m)';
    
end

% -------------------------------------------------------------------
% ----- POST-PROCESSING -----

% Final contour plot of solution
%
figure
colormap(cmap);
surf(X,Y,phi,linestyle='none')
axis equal 
view(0,90)
cb = colorbar;
ylabel(cb,'$\phi$',interpreter='latex')
xlabel('$x$',interpreter='latex')
ylabel('$y$',interpreter='latex')
fontsize(16,'points'); fontname('Serif')
box on
grid on
%}

% Final surface plot of solution
%
figure
colormap(cmap);
surf(X,Y,phi,linestyle='none')
xlabel('$x$',interpreter='latex')
ylabel('$y$',interpreter='latex')
zlabel('$\phi$',interpreter='latex')
axis equal 
view(-20,33)
colorbar
cb = colorbar;
ylabel(cb,'$\phi$',interpreter='latex')
fontsize(16,'points'); fontname('Serif')
box on
grid on
%}


% Plot centerline solution at final time against initial condition to
% compare numerical diffusion/dispersion
%
figure
hold on
plot(x,phi0(round(m/2),:),linewidth=2,marker='none',color='r')
plot(x,phi(round(m/2),:),linewidth=2,marker='none',color='b')
ylim([-0.05,0.3])
xticks(0:6); yticks(0:0.1:0.3)
xlabel('$x$',interpreter='latex')
ylabel('$\phi$',interpreter='latex')
fontsize(16,'points'); fontname('Serif')
legend('$\phi(0)$','$\phi(t)$',interpreter='latex',location='northeast')
box on
grid on
%}

% Conservation plot of total mass of pollutant at each iteration compared
% to starting mass of pollutant
%
figure

r = abs(phi_tot-phi_tot(1));
semilogy(1:time_steps,r,linewidth=2,marker='none',color='b')
hold on
semilogy([1,time_steps],[0,0]+eps,linewidth=2,marker='none',color='r')
xlabel('Iteration',interpreter='latex')
ylabel('$r_{\phi}$',interpreter='latex')
fontsize(16,'points'); fontname('Serif')
legend('$r_{\phi}$','$\varepsilon_{mach}$',interpreter='latex', ...
        location='northwest')
box on
grid on
ylim padded
%}

% Plot peak value of phi as a measure of the diffusion of the flux scheme
%
figure
plot(1:time_steps,phi_max,linewidth=2,marker='none',color='b')
xlabel('Iteration',interpreter='latex')
ylabel('$\phi_{max}$',interpreter='latex')
fontsize(16,'points'); fontname('Serif')
box on
grid on
%}