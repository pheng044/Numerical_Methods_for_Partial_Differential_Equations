% -------------------------------------------------------------------
% Patrick Heng
% 08 Mar. 2025
% Test script to study the stability of explicit time schemes for
% upwind and Lax-Wendroff flux schemes.
% -------------------------------------------------------------------

close all; clc;

% -------------------------------------------------------------------
% ----- INPUTS -----
upwind = true;
pumps = false;
diffusion = false;
step_IC = false;


% Number of nodes in the computational domain, n -> x, m -> y
n = 60;
m = 20;

t_end = 12;


D = 0.005;

q = 0;
U = @(X,Y) 2 - q*(X-2)./((X-2).^2+(Y-2).^2);
V = @(X,Y) -q*(Y-2)./((X-2).^2+(Y-2).^2);


Lx = 3;
Ly = 1;

nodes = n*m;

x = linspace(0,Lx,n);
y = linspace(0,Ly,m);

[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);

tests = 0.7:0.1:1;

i = 0;
names = cell(numel(tests),1);


for Co = tests
    dt = Co*dx/2;
    i = i + 1;
    names{i} = [num2str(Co)]; 
    if upwind == true
        F = upwind_flux(X,Y,U,V,dx,dy,n,nodes);
    else
        F = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,nodes);
    end
    
    if diffusion == true
        G = diffusive_flux(D,dx,dy,n,m,nodes);
        A = speye(nodes) - dt*F + dt*G;
    else
        A = speye(nodes) - dt*F;
    end
    
    A = periodic_BC(A,U(X(:,end),Y(:,end)),V(X(end,:),Y(end,:)),n,m);
    
    if pumps == true
        NN = @(i,j) n*(j-1) + i;
        A(NN(1:n,1),:) = 0;
        A(NN(1:n,m),:) = 0;
    end
    
    EIGS = eig(full(A));
    Re = real(EIGS);
    Im = imag(EIGS);
    plot(Re,Im,linestyle='none',marker='o',markersize=4)
    hold on

end

theta = linspace(0,2*pi,101);
theta2 = flip(theta);
shade_x = [cos(theta),4*cos(theta2)]; shade_y = [sin(theta),4*sin(theta2)];

fill(shade_x,shade_y,'r',FaceAlpha=0.15,linestyle='none');



xline(0); yline(0)
xlabel('Re($\lambda$)',interpreter='latex')
ylabel('Im($\lambda$)',interpreter='latex')

leg = legend(names,interpreter='latex',location='northeastoutside');
title(leg,'Co')
xlim([-1.3,1.3]); ylim([-1.3,1.3])
xticks(-1:1); yticks(-1:1)

axis equal
box on

fontname('Serif'); fontsize(12,'points')


