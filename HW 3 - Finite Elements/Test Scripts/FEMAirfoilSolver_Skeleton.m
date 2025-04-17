% -------------------------------------------------------------------
% Patrick Heng --> A matrix assembly and solution
% 3 Apr. 2025
% FEA/FEM script for potential flow problem over a bluff body, 
% adapted from code provided by Prof. D.J. Willis.
% -------------------------------------------------------------------

% ===================================================================
% Start by defining your domain: Meshing parameters (done for you)
% ===================================================================
% The following parameters are used for defining the mesh of the 
% object
%
% Shape: This selects an ellispoid or a NACA-4 Digit airfoil
%     -- Ellipsoid: Set the shape to a single value (the aspect 
%                   ratio of  the ellipsoid shape. eg: Shape = [1];
%    
%     -- NACA-4-Digit airfoil: Set the shape parameter a 2-valued 
%                              entry,
%     eg: Shape = [1 0012], where the second value is the 4-digit 
%         reference for the NACA airfoil.
%
% DomainSize: This is a basic parameter that defines how far the 
% domain should extend (approx.) from your shape. A value of 4 should
% be sufficient, but you can play around with this.
%
% ref: This is the refinement of the mesh. The higher the value, the 
% more elements and vertices you will have.
% 
% powerRef: This is the mesh refinment control. The larger the 
% number, the more refined the mesh is near the object. A value of 
% 1 produces a uniform refiment. A value of 1.75-2.0 should work well 
% for most of your problems.
%
% ===================================================================

close all; clear variables; clc;

% Colormaps from '200 Colormaps' on the MATLAB file exchange
% SWITCH THIS TO A DEFAULT COLORMAP ('parula') IF YOU DO NOT HAVE THE 
% FILE INSTALLED
cmap=slanCM(100);


% Viewing window for velocity plot - set 0 for full view
window = 0;

% Flag to display max velocity text
show_max_velocity = true;


% -------------------------------------------------------------------
% ----- INPUTS -----

% Forcing function for Poisson equation - accounts for 
% compressibility in the case of potential flow. Set to @ 0*X for 
% incompressible flow (Laplace equation)
RHS = @(X,Y) 0*X;

% Mesh generation parameters
Shape = [1 0065];
DomainSize = 4;
ref = 5;
powerRef = 2;


% -------------------------------------------------------------------
% ----- MESH GENERATION -----
[TRI, Nodes, DirichletNodes] = getDiscreteGeometry(Shape, ... 
                                    DomainSize, ref, powerRef);

% Plot triangulation
%
figure
trimesh(TRI, Nodes(:,1),Nodes(:,2),0*Nodes(:,1),edgecolor='k')
hold on
trimesh(TRI, Nodes(:,1),-Nodes(:,2),0*Nodes(:,1),edgecolor='k')
view([0,0,1])
axis equal


xlabel('$x$',interpreter='latex')
ylabel('$y$',interpreter='latex')
fontname('Serif'); fontsize(16,'points')
box on
%}

% -------------------------------------------------------------------
% ----- A MATRIX ASSEMBLY -----
tic
% Determine the number of nodes and elements in the domain
num_nodes = length(Nodes);
num_elem = length(TRI);

% Initialize the global A matrix and RHS
global_A_mtx = spalloc(num_nodes, num_nodes, 12*num_nodes);
f = RHS(Nodes(:,1),Nodes(:,2));   % Evaluate forcing function
F = zeros(num_nodes, 1);                % Preallocate RHS
% Matrix to solve for elemental basis coefficients
basis_mtx = ones(3,3);  
I = eye(3);                             % Identity matrix

% Stamp global A matrix with elemental stiffness matrices and 
% elemental load vectors by looping through all elements 
for i = 1:num_elem
    % Triangle node numbers
    node_ID = TRI(i,:);
    % Coordinates of triangle node numbers
    basis_mtx(:,1) = Nodes(node_ID,1);
    basis_mtx(:,2) = Nodes(node_ID,2);

    % Element basis function coefficients:
    % Linear element basis: phi(x,y) = ax + by + c
    basis_coeffs = basis_mtx\I; 
    % Get a and b coefficients for stiffness matrix calculation
    a = basis_coeffs(1,:);
    b = basis_coeffs(2,:);
    
    % Calculate element size (area)
    side1 = [Nodes(node_ID(1),:) - Nodes(node_ID(2),:),0];
    side2 = [Nodes(node_ID(2),:) - Nodes(node_ID(3),:),0];
    tri_area = 0.5*norm(cross(side1,side2));

    % Inner product of basis function gradients (stiffness matrix)
    elem_mtx = tri_area*(a'*a + b'*b);
    
    % Place elemental stiffness matrix into global A matrix
    global_A_mtx(node_ID,node_ID) = global_A_mtx(node_ID,node_ID) ...
                                        + elem_mtx;
    % Place elemental load vector into global load vector
    F(node_ID) = F(node_ID) - tri_area*f(node_ID)/3;
end

% -------------------------------------------------------------------
% ----- BC's -----

% For Dirichlet nodes, replace the rows corresponding to these nodes
% with the same rows of the identity matrix. For load vector, replace
% the Dirichlet nodes with desired boundary value
I = speye(num_nodes);
global_A_mtx(DirichletNodes,:) = I(DirichletNodes,:);
F(DirichletNodes) = Nodes(DirichletNodes,1);

% -------------------------------------------------------------------
% ----- SOLVE ----

phi = global_A_mtx\F;
toc

% -------------------------------------------------------------------
% ----- POST PROCESSING -----

% Scalar potential value at foward stagnation point (0,0)
phi_stag = phi(Nodes(Nodes(:,1)==0,2)==0)


% ===================================================================
% Gradient calculation
% ===================================================================
for i=1:num_elem
    element_number = i;

    N1 = TRI(i,1);
    N2 = TRI(i,2);
    N3 = TRI(i,3);
    
    X1 = Nodes(N1,1);
    X2 = Nodes(N2,1);
    X3 = Nodes(N3,1);
    
    Y1 = Nodes(N1,2);
    Y2 = Nodes(N2,2);
    Y3 = Nodes(N3,2);

    S1 = phi(N1);
    S2 = phi(N2);
    S3 = phi(N3);

    C = [1 X1 Y1; 1 X2 Y2; 1 X3 Y3]\[eye(3)];
    
    Gradient_IE(i,:) = [(S1*C(2,1)+S2*C(2,2)+S3*C(2,3)), ...
        (S1*C(3,1)+S2*C(3,2)+S3*C(3,3))];
    
    Centroid(i,:) = [(X1+X2+X3)/3, (Y1+Y2+Y3)/3];
end

% ===================================================================
% Figure 1
% ===================================================================
figure
trisurf(TRI, Nodes(:,1), Nodes(:,2), phi-500)
hold on
trisurf(TRI, Nodes(:,1),-Nodes(:,2),phi-500)

quiver(Centroid(:,1), Centroid(:,2), ...
            Gradient_IE(:,1), Gradient_IE(:,2),.75,'k')
quiver(Centroid(:,1), -Centroid(:,2), ...
            Gradient_IE(:,1), -Gradient_IE(:,2),.75,'k')

shading interp
%title(['Scalar Velocity Potential Distribution ' ...
%            '(Nodal) and the Velocity Vector Field'])
view([0 0 1])
axis equal
colormap(cmap)
cb = colorbar;
xlim([min(Nodes(:,1)),max(Nodes(:,1))])
ylim([-max(Nodes(:,2)),max(Nodes(:,2))])

box on
grid on
xlabel('$x$',interpreter='latex')
ylabel('$y$',interpreter='latex')
zlabel('$\phi$',interpreter='latex')
ylabel(cb,'$\phi$',interpreter='latex')
fontname('Serif'); fontsize(16,'points')

% ===================================================================
% Figure 2
% ===================================================================
figure
Vel = ((Gradient_IE(:,1).^2 + Gradient_IE(:,2).^2).^(.5))';
a = trisurf(TRI, Nodes(:,1), Nodes(:,2), Nodes(:,2)*0, Vel);
hold on

set(a,'edgealpha',0)
trisurf(TRI, Nodes(:,1), -Nodes(:,2), -Nodes(:,2)*0, Vel, ...
                linestyle='none')
%title('Velocity Distribution (Centroidal)')
view([0 0 1])
axis equal

colormap(cmap)
cb = colorbar;

if window == 0
    xlim([min(Nodes(:,1)),max(Nodes(:,1))])
    ylim([-max(Nodes(:,2)),max(Nodes(:,2))])
else
    xlim([-window,window])
    ylim([0,window])
end

if show_max_velocity == true
    u_max = max(Vel,[],'all');
    annotation('textbox',[0.19,0.65,0.1,0.1],...
        'string',['$u_{max}=' num2str(u_max,3) '$'],...
        interpreter='latex',backgroundcolor='w', ...
        verticalalignment='middle', ...
        HorizontalAlignment='center')
end

box on
grid on
xlabel('$x$',interpreter='latex')
ylabel('$y$',interpreter='latex')
zlabel('$u$',interpreter='latex')
ylabel(cb,'$u$',interpreter='latex')
fontname('Serif'); fontsize(16,'points')
