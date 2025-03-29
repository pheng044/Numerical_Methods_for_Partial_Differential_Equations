% -------------------------------------------------------------------
% 'HW_2_pumps_test.m'
% Patrick Heng
% 27 Mar. 2025
% Test script to run the finite volume solver for the convection-
% diffusion equations with potential flow sink.
% -------------------------------------------------------------------

close all; clear all; clc;

% Colormaps from '200 Colormaps' on the MATLAB file exchange
% SWITCH THIS TO A DEFAULT COLORMAP ('parula') IF YOU DO NOT HAVE THE 
% FILE INSTALLED
cmap=slanCM(2);

% -------------------------------------------------------------------
% ----- INPUTS -----
upwind = true;
pumps = false;
diffusion = false;
step_IC = false;
periodic = false;


% Number of nodes (cell centers) in the computational domain, 
% n -> x, m -> y
n = 240;
m = 40;

t_end = 20;


% Courant number, used to determine time step size
Co = 0.7;

% Initial concentration distribution
IC = @(X,Y) 0.25*exp(-50*((X-0.5).^2+(Y-0.5).^2));

% Diffusion constant
D = 0.005;

% Sink strength, 0 = off
Q = 0:0.1:1;

% Flow velocity components in absence of sources/sinks
u0 = 2;
v0 = 0;

% Length of domain
Lx = 6;
Ly = 1;

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
%


figure
hold on
step = floor(256/length(Q));
col = cmap(1:step:256,:);

i = 1;
for q = Q
% Velocity components as functions of the grid points
U = @(X,Y) u0 - q*(X-2)./((X-2).^2+(Y-2).^2);
V = @(X,Y) v0 - q*(Y-2)./((X-2).^2+(Y-2).^2);


% Determine flux function
if upwind == true
    F = upwind_flux(X,Y,U,V,dx,dy,n,m,nodes,periodic);
else
    [Fx,Fy] = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,m,nodes,periodic);
end

% Turn on diffusion
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


% Determine the amount of time steps to get to tf
time_steps = ceil(t_end/dt);

% Store total pollutant mass, peak value of concentration
phi_tot = zeros(time_steps,1);

% Evaluate initial condition at cell centers
phi = IC(X,Y);
phi = reshape(phi',nodes,1);

% ----- TIME INTEGRATION -----
% Loop through time in an explicit marching scheme
for k = 1:time_steps
    
    % Store current pollutant mass, peak amplitude
    phi_tot(k) = sum(phi,'all')*dx*dy;
    % update solution
    phi = A*phi;
    
end

% Conservation plot of total mass of pollutant at each iteration
%

plot(0:dt:t_end,phi_tot,color=col(i,:),linewidth=2)
i = i + 1;
%plot(0:dt:t_end,phi_tot,linewidth=2,marker='none')
end

xlabel('$t$ (s)',interpreter='latex')
ylabel('$\phi_{tot}$',interpreter='latex')
ylim([0,1.01*max(phi_tot)])

cb = colorbar;
colormap(cmap);
ylabel(cb,'$q$',Interpreter='latex')
clim([0,1])

fontsize(16,'points'); fontname('Serif')
box on
grid on
%}
