% -------------------------------------------------------------------
% Patrick Heng
% 23 Feb. 2025
% Test script to run error convergence from our homework. Essentially 
% just calling all of the previously developed functions in a BIG loop.
% -------------------------------------------------------------------
close all; clear all; clc;

% -------------------------------------------------------------------
% ----- INPUTS -----
tic
% Forcing function as a function of the physical meshgrid variables,
% 0 for the Laplace equation
F = @(X,Y) 0*X;

% Define mesh shape
h_b = 3;
h_c = 2;

% -------------------------------------------------------------------
% ----- SOLVER -----
% BIG loop to save mesh data to study error convergence.
% Grab a coffee, this may take a while to run unless you have a super
% computer
for n = 2.^(2:10)

    % Number of nodes in the computational domain, n -> x, m -> y
    m = n;
    
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
    [X,Y,coeffs_x,coeffs_y] = bilinear_map(Xi,Eta,[-4,4,4,-4], ...
                                [0,0,h_b,h_c]);
    
    % Compute derivatives of the coordinate transformation
    [x_xi,x_eta,y_xi,y_eta,x_xi_eta,...
                y_xi_eta,J,a,b,c,d,e,alpha,beta] = ...
                bilinear_map_derivatives(Eta,Xi,coeffs_x,coeffs_y);
    
    % Generate 9 point stencil for the main diagonals of the A matrix 
    % (discrete Laplacian)
    diag = stencil_coefficients(J,a,b,c,d,e,dxi,deta,nodes);
    diag_idx = [-n-1, -n, -n+1, -1, 0, 1, n-1,n , n+1];
    
    % Compute the A matrix over the internal nodes. Initially generate the 
    % A matrix with one extra column to get the correct indexing for the 
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
    
    % Impose boundary conditions on the A matrix, generate the 
    % corresponding forcing vector, f
    f = reshape(F(X,Y),nodes,1);
    
    [A,f] = impose_BCs(A,f,x_xi,x_eta,y_xi,y_eta,J,dxi,deta,X);
    
    % Linear solve for nodes
    u = A\f;
    
    % Reshape u to match meshgrid variables
    u = flip(reshape(u,n,m)');
    
    % Save the physical mesh and solution
    save(['u_mesh_' num2str(n) '_' num2str(m) '.mat'],'X','Y','u')

end
toc
