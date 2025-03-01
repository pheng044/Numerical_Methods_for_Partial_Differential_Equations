% -------------------------------------------------------------------
% Patrick Heng
% 8 Feb. 2025
% Function to impose boundary conditions for our specific HW problem.
% Neumann, No flux: Top, Left, Right, Bottom (partial)
% Dirichlet: Bottom over x = [-2,0] and [0,2] in the physcial domain
% Return modified A matrix, f vector 
% -------------------------------------------------------------------

function [A,f] = impose_BCs(A,f,x_xi,x_eta,y_xi,y_eta,J,dxi,deta,X)
    
    % Get domain size
    n = numel(X(1,:));
    m = numel(X(:,1));
    nodes = n*m;

    % ----- NEUMANN CONDITIONS -----

    % --- TOP BOUNDARY ---

    % Compute transformed normal derivative coefficients 
    % (notation defined in report)
    idx = NN(1:n,1,n);
    [N_eta_xi, N_eta_eta] = normals_eta(x_xi(idx), x_eta(idx), ...
        y_xi(idx), y_eta(idx), J(idx));

    % Central xi, foward eta derivative stencil and indexing
    diag = [-N_eta_xi/(2*dxi), -3*N_eta_eta/(2*deta), ... 
        N_eta_xi/(2*dxi), 4*N_eta_eta/(2*deta), -N_eta_eta/(2*deta)];
    diag_idx = [-1, 0, 1, n, 2*n];
    
    % Clear rows of A corresponding to the boundary
    A(idx,:) = 0;

    % Apply stencil to A matrix, avoid domain corners
    for i = 2:n-1
        A(NN(i,1,n),NN(i,1,n)+diag_idx) = diag(i,:);
    end


    % --- BOTTOM BOUNDARY ---
    
    % Compute transformed normal derivative coefficients 
    idx = NN(1:n,m,n);

    [N_eta_xi, N_eta_eta] = normals_eta(x_xi(idx), ... 
        x_eta(idx), y_xi(idx), y_eta(idx), J(idx));
    
    % Central xi, backward eta derivative stencil and indexing
    diag = [N_eta_eta/(2*deta), -4*N_eta_eta/(2*deta), ... 
        -N_eta_xi/(2*dxi), 3*N_eta_eta/(2*deta), N_eta_xi/(2*dxi)];
    diag_idx = [-2*n, -n, -1, 0, 1];
    
    
    % Clear rows of A corresponding to the boundary
    A(idx,:) = 0;

    % Apply stencil to A matrix, avoid domain corners
    for i = 2:n-1
        A(NN(i,m,n),NN(i,m,n)+diag_idx) = diag(i,:);
    end
    
    % --- LEFT BOUNDARY ---

    % Compute transformed normal derivative coefficients 
    idx = NN(1,1:m,n);
    [N_xi_xi, N_xi_eta] = normals_xi(x_xi(idx), ...
        x_eta(idx),y_xi(idx),y_eta(idx), J(idx));
    
    % Forward xi, central eta derivative stencil and indexing
    diag = [-N_xi_eta/(2*deta), -3*N_xi_xi/(2*dxi), ... 
        4*N_xi_xi/(2*dxi), -N_xi_xi/(2*dxi), N_xi_eta/(2*deta)];
    diag_idx = [-n, 0, 1, 2, n];
    
    % Clear rows of A corresponding to the boundary
    A(NN(1,2:m-1,n),:) = 0;
    
    % Apply stencil to A matrix, avoid domain corners
    for j = 2:m-1
        A(NN(1,j,n),NN(1,j,n)+diag_idx) = diag(j,:);
    end
    
    % --- RIGHT BOUNDARY ---

    % Compute transformed normal derivative coefficients 
    idx = NN(n,1:m,n);
    [N_xi_xi, N_xi_eta] = normals_xi(x_xi(idx), ... 
        x_eta(idx), y_xi(idx), y_eta(idx), J(idx));
    
    % Backward xi, central eta derivative stencil and indexing
    diag = [-N_xi_eta/(2*deta), N_xi_xi/(2*dxi), -4*N_xi_xi/(2*dxi), ... 
        3*N_xi_xi/(2*dxi), N_xi_eta/(2*deta)];
    diag_idx = [-n, -2, -1, 0, n];
    
    % Clear rows of A corresponding to the boundary
    A(NN(n,2:m-1,n),:) = 0;
    
    % Apply stencil to A matrix, avoid domain corners
    for j = 2:m-1
        A(NN(n,j,n),NN(n,j,n)+diag_idx) = diag(j,:);
    end
    
    
    % ----- CORNER POINTS -----
    % Mix of foward and backward differencing to maintain second order
    % accuracy. F = forwards, B = backwards

    % --- TOP LEFT ---

    % F xi, F eta
    diag = [-3*N_xi_xi(1)/(2*dxi)+3*N_xi_eta(1)/(2*deta), ...
        4*N_xi_xi(1)/(2*dxi), -N_xi_xi(1)/(2*dxi), ...
        -4*N_xi_eta(1)/(2*deta), N_xi_eta(1)/(2*deta)];
    diag_idx = [0,1,2,n,2*n];
    A(NN(1,1,n),:) = 0;
    A(NN(1,1,n), 1 + diag_idx) = diag;
    
    
    % --- BOTTOM LEFT ---

    % F xi, B eta
    diag = 0.5*[(N_xi_eta(end)+N_eta_eta(1))/(2*deta), ...
            -4*(N_xi_eta(end)+N_eta_eta(1))/(2*deta), ...
            -3*(N_xi_xi(end)+N_eta_xi(1))/(2*dxi) + ...
             3*(N_xi_eta(end)+N_eta_eta(1))/(2*deta), ...
             4*(N_xi_xi(end)+N_eta_xi(1))/(2*dxi), ...
            -(N_xi_xi(end)+N_eta_xi(1))/(2*dxi)];
    diag_idx = [-2*n,-n,0,1,2];
    A(NN(1,m,n),:) = 0;
    A(NN(1,m,n), NN(1,m,n) + diag_idx) = diag;
    
    
    % --- TOP RIGHT ---

    % B xi, F eta
    diag = [N_xi_xi(1)/(2*dxi), -4*N_xi_xi(1)/(2*dxi), ...
        3*N_xi_xi(1)/(2*dxi)-3*N_xi_eta(1)/(2*deta), ...
        4*N_xi_eta(1)/(2*deta), -N_xi_eta(1)/(2*deta)];
    diag_idx = [-2,-1,0,n,2*n];
    A(NN(n,1,n),:) = 0;
    A(NN(n,1,n), NN(n,1,n) + diag_idx) = diag;
    
    
    % --- BOTTOM RIGHT ---

    % B xi, B eta
    diag = 0.5*[(N_xi_eta(end)+N_eta_xi(end))/(2*deta), ...
            -4*(N_xi_eta(end)+N_eta_xi(end))/(2*deta), ...  
             (N_xi_xi(end)+N_eta_xi(end))/(2*dxi), ...
            -4*(N_xi_xi(end)+N_eta_xi(end))/(2*dxi), ...
             3*(N_xi_xi(end)+N_eta_xi(end))/(2*dxi) + ...
             3*(N_xi_eta(end)+N_eta_xi(end))/(2*deta)];
    diag_idx = [-2*n,-n,-2,-1,0];
    A(NN(n,m,n),:) = 0;
    A(NN(n,m,n), NN(n,m,n) + diag_idx) = diag;
    
    
    % ----- DIRICHLET CONDITIONS ----
    % Find the indices of the physical domain where Dirichlet boundaries 
    % should be imposed, use logical indexing
    
    indices = 1:n;

    % Left DBC
    dirich_idx_left = indices(-2 <= X(1,:) & X(1,:) <= 0);
    % Right DBC
    dirich_idx_right = indices(0 <= X(1,:) & X(1,:) <= 2);
    
    % Replace the rows of the Dirichlet nodes with the corresponding rows 
    % of the identity matrix
    
    I = speye(nodes);
    A(NN(dirich_idx_left,m,n),:) = I(NN(dirich_idx_left,m,n),:);
    A(NN(dirich_idx_right,m,n),:) = I(NN(dirich_idx_right,m,n),:);

    % ----- FORCING FUCNTION -----
    % Neumann conditions for RHS
    f(NN(1:n,1,n)) = 0;     % Top
    f(NN(1:n,m,n)) = 0;     % Bottom
    f(NN(1,1:m,n)) = 0;     % Left
    f(NN(n,1:m,n)) = 0;     % Right
    
    % Impose Dirichlet conditions for the source
    f(NN(dirich_idx_left,m,n)) = -1;
    f(NN(dirich_idx_right,m,n)) = 1;


    % Node numbering function, i -> xi, j -> eta, n = total xi nodes
    function NN = NN(i,j,n)
        NN = n*(j-1) + i;
    end

end