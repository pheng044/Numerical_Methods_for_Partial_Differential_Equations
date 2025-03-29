% -------------------------------------------------------------------
% 'HW_2_A_Eigenvalues.m'
% Patrick Heng
% 08 Mar. 2025
% Test script to study the stability of explicit time schemes for
% upwind and Lax-Wendroff flux schemes.
% -------------------------------------------------------------------

clear all; close all; clc;

% -------------------------------------------------------------------
% ----- INPUTS -----
upwind = false;
diffusion = false;
step_IC = false;


% Number of nodes in the computational domain, n -> x, m -> y
n = 60;
m = 30;

% Simulation time
t_end = 1;

% Courant numbers, used to determine time step size
tests = 0.7:0.1:1.1;

% Diffusion constant
D = 0.005;

% Domain lengths
Lx = 2;
Ly = 1;

% Sink strength, 0 = off
q = 0;

% Flow velocity components in absence of sources/sinks
u0 = 2;
v0 = 0;

% Velocity components as functions of the grid points
U = @(X,Y) u0 - q*(X-2)./((X-2).^2+(Y-2).^2);
V = @(X,Y) v0 - q*(Y-2)./((X-2).^2+(Y-2).^2);

% Total number of cell centers
nodes = n*m;

% Spatial discretization steps
dx = Lx/n;
dy = Ly/m;

% Create grid of cell centers
x = linspace(dx/2,Lx-dx/2,n);
y = linspace(dy/2,Ly-dy/2,m);
[X,Y] = meshgrid(x,y);

% Set up counter
i = 0;
% Store Courant number sfor legend
names = cell(numel(tests),1);

% Loop through Courant numbers
for Co = tests
    % Calculate time step
    dt = Co*dx/2;
    % Increment counter
    i = i + 1;
    % At Co to legend names
    names{i} = [num2str(Co)];
    
    % Covective flux update matrices
    if upwind == true
        F = upwind_flux(X,Y,U,V,dx,dy,n,m,nodes,'periodic');
    else
        [Fx,Fy] = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,m,nodes, ...
                                        'periodic');
    end
    
    % Calculate diffusive flux update matrix
    if diffusion == true
        G = diffusive_flux(D,dx,dy,n,m,nodes);
    else
        G = 0;
    end
    
    % Total update matrix
    if upwind == true
        A = speye(nodes) - dt*F + dt*G;
    else
        Ax = speye(nodes) - dt*Fx;
        Ay = speye(nodes) - dt*Fy + dt*G;
        A = Ay*Ax;
    end
    
    % Get eigenvalues of A matrix
    EIGS = eig(full(A));
    Re = real(EIGS);
    Im = imag(EIGS);
    plot(Re,Im,linestyle='none',marker='o',markersize=4)
    hold on

end

% Plot instable region in red
theta = linspace(0,2*pi,101);
theta2 = flip(theta);
shade_x = [cos(theta),4*cos(theta2)]; 
shade_y = [sin(theta),4*sin(theta2)];
fill(shade_x,shade_y,'r',FaceAlpha=0.15,linestyle='none');

% Plot principle axes
xline(0); yline(0)

% Axis limits a = size of limits
a = 1.5;
xlim([-a,a]); ylim([-a,a])
xticks(-a:a); yticks(-a:a)

% Pretty plot parameters
xlabel('Re($\lambda$)',interpreter='latex')
ylabel('Im($\lambda$)',interpreter='latex')
leg = legend(names,interpreter='latex',location='northeastoutside');
title(leg,'Co')
axis equal
box on
fontname('Serif'); fontsize(12,'points')
