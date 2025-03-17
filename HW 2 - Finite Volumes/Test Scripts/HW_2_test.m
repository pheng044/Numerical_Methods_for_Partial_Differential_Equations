% -------------------------------------------------------------------
% Patrick Heng
% 08 Mar. 2025
% Test script to run the finite volume solver for the convection-
% diffusion equations.
% -------------------------------------------------------------------

close all; clear all; clc;

% -------------------------------------------------------------------
% ----- INPUTS -----
upwind = false;
pumps = false;
diffusion = false;
step_IC = false;


% Number of nodes in the computational domain, n -> x, m -> y
n = 240;
m = 40;

t_end = 12;

Co = 0.9;

IC = @(X,Y) 0.25*exp(-50*((X-0.5).^2+(Y-0.5).^2));
D = 0.005;
q = 0;

U = @(X,Y) 2 - q*(X-2)./((X-2).^2+(Y-2).^2);
V = @(X,Y) -q*(Y-2)./((X-2).^2+(Y-2).^2);


Lx = 6;
Ly = 1;

nodes = n*m;

x = linspace(0,Lx,n);
y = linspace(0,Ly,m);

[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
dt = Co*dx/2;

if upwind == true
    F = upwind_flux(X,Y,U,V,dx,dy,n,nodes);
else
    [Fx,Fy] = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,nodes);
end

if diffusion == true
    G = diffusive_flux(D,dx,dy,n,m,nodes);
    if upwind == true
        A = speye(nodes) - dt*F + dt*G;
    else
        Ax = speye(nodes) - dt*Fx + dt*G;
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

A = periodic_BC(A,U(X(:,end),Y(:,end)),V(X(end,:),Y(end,:)),n,m);

if pumps == true
    NN = @(i,j) n*(j-1) + i;
    A(NN(1:n,1),:) = 0;
    A(NN(1:n,m),:) = 0;
end

phi = IC(X,Y);

if step_IC == true
    phi(phi<0.0025) = 0;
    phi(phi>=0.0025) = 0.25;
end


time_steps = ceil(t_end/dt);
phi_tot = zeros(time_steps,1);


% Loop through time in an explicit marching scheme
for k = 1:time_steps
    
    phi_tot(k) = sum(phi(1:end-1,1:end-1),'all')*dx*dy;
    %{
    surf(X,Y,phi,linestyle='none')
    axis equal
    view(0,90)
    colorbar
    %}
    plot(x,phi(round(m/2),:),linewidth=2,marker='none',color='b')
    ylim([-0.05,0.3])
    xticks(0:6); yticks(0:0.1:0.3)
    xlabel('$x$',interpreter='latex')
    ylabel('$\phi$',interpreter='latex')
    fontsize(16,'points'); fontname('Serif')
    pause(0.1)

    phi = reshape(phi',nodes,1);
    
  
    phi = A*phi;
    phi = reshape(phi,n,m)';
    
end


box on
%{
figure
plot(1:time_steps,abs(phi_tot(1)-phi_tot)/phi_tot(1),linewidth=2,marker='none',color='b')
xlabel('Interation',interpreter='latex')
ylabel('$r_{\phi}$',interpreter='latex')
fontsize(16,'points'); fontname('Serif')
box on
%}
%pause(0.1)