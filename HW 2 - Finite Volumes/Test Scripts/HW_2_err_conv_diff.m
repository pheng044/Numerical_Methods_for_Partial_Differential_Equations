% -------------------------------------------------------------------
% Patrick Heng
% 28 Mar. 2025
% Test script to run the error data generation of finite volume solver 
% for the convection-diffusion equations.
% -------------------------------------------------------------------

close all; clear all; clc;

% -------------------------------------------------------------------
% ----- INPUTS -----
upwind = false;
diffusion = true;
step_IC = false;
periodic = true;

% Simulation length
t_end = 0.5;

% Courant number, used to determine time step size
Co = 0.1;

% Initial concentration distribution
IC = @(X,Y) 0.25*exp(-50*((X-0.5).^2+(Y-0.5).^2));

% Flow velocity components in absence of sources/sinks
u0 = 2;
v0 = 0;

% Velocity conmponents as functions of the grid points
q = 0;
U = @(X,Y) u0 - q*(X-2)./((X-2).^2+(Y-2).^2);
V = @(X,Y) v0 - q*(Y-2)./((X-2).^2+(Y-2).^2);


for kk = 2.^(2:9)

    % Number of nodes (cell centers) in the computational domain, 
    % n -> x, m -> y
    n = 2*kk;
    m = kk;
    
    % Length of domain
    Lx = 2;
    Ly = 1;
    
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
        A = speye(nodes) - dt*F;
    else
        [Fx,Fy] = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,m,nodes,periodic);
        Ax = speye(nodes) - dt*Fx;
        Ay = speye(nodes) - dt*Fy;
        A = Ay*Ax;
    end
    
    phi = IC(X,Y);
    time_steps = ceil(t_end/dt);
    
    % Loop through time in an explicit marching scheme
    for k = 1:time_steps
        phi = reshape(phi',nodes,1);
        phi = A*phi;
        phi = reshape(phi,n,m)';
    end
    % Save the physical mesh and solution
    save(['phi_LWD_' num2str(n) '_' num2str(m) '.mat'],'X','Y','phi','dx','dy')

end