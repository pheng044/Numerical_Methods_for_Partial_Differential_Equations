% -------------------------------------------------------------------
% Patrick Heng
% 23 Feb. 2025
% Test script to run the electrolocation study. Essentially just
% calling all of the previously developed functions in a BIG loop.
% -------------------------------------------------------------------
close all; clear all; clc;

% Generate the normal distances and angles to test
H_n = linspace(1,10,75);
theta = linspace(0,pi/2-0.01,75);

[HN,THETA] = meshgrid(H_n,theta);

% Preallocate min/max values of E field 
Ey_max = NaN(size(H_n,2),size(theta,2));
Ey_min = NaN(size(H_n,2),size(theta,2));

% -------------------------------------------------------------------
% ----- INPUTS -----
tic
% Forcing function as a function of the physical meshgrid variables,
% 0 for the Laplace equation
F = @(X,Y) 0*X;
% Number of nodes in the computational domain, n -> x, m -> y
n = 64;
m = 64;

for i = 1:size(H_n,2)
    for j = 1:size(theta,2)
        % Define mesh shape
        h_b = H_n(i) + 4*tan(theta(j));
        h_c = H_n(i) - 4*tan(theta(j));
        if h_c < 0
            Ey_max(i,j) = NaN;
            break
        end
        
        
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
        
        Ey_max(i,j) = max(Ey(1,:),[],'all');
        Ey_min(i,j) = min(Ey(1,:),[],'all');
    end
end
toc

dEy = Ey_max-Ey_min;

surf(HN,THETA*180/pi,dEy',linestyle='none')

% -------------------------------------------------------------------
% IF YOU DESIRE TO SAVE THE SOLUTION: UNCOMMENT, RENAME
%save('Electrolocation_study_data.mat','HN','THETA','dEy', 
%       ... 'Ey_max','Ey_min')
