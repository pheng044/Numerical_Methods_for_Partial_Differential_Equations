% -------------------------------------------------------------------
% Patrick Heng
% 9 Feb. 2025
% Function to evaluate the gradient of a potential field, u, over a 
% uniform computational grid given the mapping derivatives to physical
% domain (notation defined in report).
% -------------------------------------------------------------------

function [Ex,Ey] = grad_field(u,x_xi,x_eta,y_xi,y_eta,J,dxi,deta,n,m)
    
    % Get total number of nodes, n -> xi, m -> eta
    nodes = n*m;
    
    % Make sure u is a column vector
    %u = flip(u,2);

    % C = Central difference, F = Forward difference, 
    % B = Backward difference
    
    % ----- X COMPONENT OF GRAD FIELD -----

    % Generate the main stencil for Ex with central differencing
    diag = 0.5.*(J.^-1).*[y_xi/deta, -y_eta/dxi, y_eta/dxi, -y_xi/deta];
    diag_idx = [-n,-1,1,n];

    % Place the stencil coefficients on the appropriate diagonals
    A = spdiags(diag,diag_idx,nodes,nodes+1);
    A = A(:,1:nodes);

    % Generate the boundary stencils for Ex with mixed central, forward,
    % and backwards differencing

    % --- BOTTOM ---
    idx = NN(1:n,m,n);
    
    % xi - C, eta - B
    diag = 0.5.*(J(idx).^-1).*[-y_xi(idx)/deta,4*y_xi(idx)/deta, ...
        -y_eta(idx)/dxi,-3*y_xi(idx)/deta,y_eta(idx)/dxi];
    diag_idx = [-2*n,-n,-1,0,1]; 
    
    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:n-1
        A(NN(i,m,n),NN(i,m,n)+diag_idx) = diag(i,:);
    end

    % --- TOP ---
    idx = NN(1:n,1,n);

    % xi - C, eta - F
    diag = 0.5*(J(idx).^-1).*[-y_eta(idx)/dxi,3*y_xi(idx)/deta, ...
        y_eta(idx)/dxi,-4*y_xi(idx)/deta,y_xi(idx)/deta];
    diag_idx = [-1,0,1,n,2*n]; 
    
    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:n-1
        A(NN(i,1,n),NN(i,1,n)+diag_idx) = diag(i,:);
    end

    % --- LEFT ---
    idx = NN(1,1:m,n);

    % xi - F, eta - C
    diag = 0.5*(J(idx).^-1).*[y_xi(idx)/deta,-3*y_eta(idx)/dxi, ...
        4*y_eta(idx)/dxi,-y_eta(idx)/dxi,-y_xi(idx)/deta];
    diag_idx = [-n,0,1,2,n]; 

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:m-1
        A(NN(1,i,n),NN(1,i,n)+diag_idx) = diag(i,:);
    end

    % --- RIGHT ---
    idx = NN(n,1:m,n);

    % xi - B, eta - C
    diag = 0.5*(J(idx).^-1).*[y_xi(idx)/deta,y_eta(idx)/dxi, ...
        -4*y_eta(idx)/dxi,3*y_eta(idx)/dxi,-y_xi(idx)/deta];

    diag_idx = [-n,-2,-1,0,n]; 
    
    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:m-1
        A(NN(n,i,n),NN(n,i,n)+diag_idx) = diag(i,:);
    end
    
    % --- BOTTOM RIGHT ---
    idx = NN(n,m,n);
    % xi - B, eta - B 
    diag = 0.5*(J(idx).^-1).*[-y_xi(idx)/deta,4*y_xi(idx)/deta, ...
            y_eta(idx)/dxi,-4*y_eta(idx)/dxi, ...
            3*y_eta(idx)/dxi-3*y_xi(idx)/deta];
    diag_idx = [-2*n,-n,-2,-1,0];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag; 

    % --- TOP RIGHT ---
    idx = NN(n,1,n);
    % xi - B, eta - F
    diag = 0.5*(J(idx).^-1).*[y_eta(idx)/dxi,-4*y_eta(idx)/dxi, ...
            3*y_eta(idx)/dxi+3*y_xi(idx)/deta,-4*y_xi(idx)/deta, ...
            y_xi(idx)/deta];
    diag_idx = [-2,-1,0,n,2*n];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag;

    % --- BOTTOM LEFT ---
    idx = NN(1,m,n);
    % xi - F, eta - B
    diag = 0.5*(J(idx).^-1).*[-y_xi(idx)/deta,4*y_xi(idx)/deta, ...
            -3*y_eta(idx)/dxi-3*y_xi(idx)/deta,4*y_eta(idx)/dxi, ...
            -y_eta(idx)/dxi];
    diag_idx = [-2*n,-n,0,1,2];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag;

    % --- TOP LEFT ---
    idx = NN(1,1,n);
    % xi - F, eta F
    diag = 0.5*(J(idx).^-1).*[-3*y_eta(idx)/dxi+3*y_xi(idx)/deta, ... 
            4*y_eta(idx)/dxi,-y_eta(idx)/dxi,-4*y_xi(idx)/deta, ...
            y_xi(idx)/deta];
    diag_idx = [0,1,2,n,2*n];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag;

    
    % Generate Ex by applying the finite difference matrix
    Ex = A*u;

    % ----- Y COMPONENT OF GRAD FIELD -----

    % Generate the stencil for Ey with central differencing
    diag = 0.5.*(J.^-1).*[-x_xi/deta, x_eta/dxi, -x_eta/dxi, x_xi/deta];
    diag_idx = [-n,-1,1,n];

    % Place the stencil coefficients on the appropriate diagonals
    A = spdiags(diag,diag_idx,nodes,nodes+1);
    A = A(:,1:nodes);

    % Generate the boundary stencils for Ey with mixed central, forward,
    % and backwards differencing

    % --- BOTTOM ---
    idx = NN(1:n,m,n);
    
    % xi - C, eta - B
    diag = -0.5.*(J(idx).^-1).*[-x_xi(idx)/deta,4*x_xi(idx)/deta, ...
            -x_eta(idx)/dxi,-3*x_xi(idx)/deta,x_eta(idx)/dxi];
    diag_idx = [-2*n,-n,-1,0,1]; 

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:n-1
        A(NN(i,m,n),NN(i,m,n)+diag_idx) = diag(i,:);
    end

    % --- TOP ---
    idx = NN(1:n,1,n);
    
    % xi - C, eta - F
    diag = -0.5*(J(idx).^-1).*[-x_eta(idx)/dxi,3*x_xi(idx)/deta, ...
            x_eta(idx)/dxi,-4*x_xi(idx)/deta,x_xi(idx)/deta];
    diag_idx = [-1,0,1,n,2*n]; 

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:n-1
        A(NN(i,1,n),NN(i,1,n)+diag_idx) = diag(i,:);
    end

    % --- LEFT ---
    idx = NN(1,1:m,n);

    % xi - F, eta - C
    diag = -0.5*(J(idx).^-1).*[x_xi(idx)/deta,-3*x_eta(idx)/dxi, ...
            4*x_eta(idx)/dxi,-x_eta(idx)/dxi,-x_xi(idx)/deta];
    diag_idx = [-n,0,1,2,n]; 

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:m-1
        A(NN(1,i,n),NN(1,i,n)+diag_idx) = diag(i,:);
    end

    % --- RIGHT ---
    idx = NN(n,1:m,n);
    
    % xi - B, eta - C
    diag = -0.5*(J(idx).^-1).*[x_xi(idx)/deta,x_eta(idx)/dxi, ... 
            -4*x_eta(idx)/dxi, 3*x_eta(idx)/dxi,-x_xi(idx)/deta];
    diag_idx = [-n,-2,-1,0,n]; 

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary

    % Fill the rows corresponding to the boundary
    for i = 2:m-1
        A(NN(n,i,n),NN(n,i,n)+diag_idx) = diag(i,:);
    end

    % --- BOTTOM RIGHT ---
    idx = NN(n,m,n);
    % xi - B, eta - B
    diag = 0.5*(J(idx).^-1).*[-x_xi(idx)/deta,4*x_xi(idx)/deta, ...
            x_eta(idx)/dxi,-4*x_eta(idx)/dxi, ...
            3*x_eta(idx)/dxi-3*x_xi(idx)/deta];
    diag_idx = [-2*n,-n,-2,-1,0];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag;

    % --- TOP RIGHT ---
    idx = NN(n,1,n);
    % xi - B, eta - F
    diag = 0.5*(J(idx).^-1).*[x_eta(idx)/dxi,-4*x_eta(idx)/dxi, ...
            3*x_eta(idx)/dxi+3*x_xi(idx)/deta, -4*x_xi(idx)/deta, ...
            x_xi(idx)/deta];
    diag_idx = [-2,-1,0,n,2*n];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag;

    % --- BOTTOM LEFT ---
    idx = NN(1,m,n);
    % xi - F, eta - B
    diag = 0.5*(J(idx).^-1).*[-x_xi(idx)/deta,4*x_xi(idx)/deta, ...
            -3*x_eta(idx)/dxi-3*x_xi(idx)/deta,4*x_eta(idx)/dxi, ...
            -x_eta(idx)/dxi];
    diag_idx = [-2*n,-n,0,1,2];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag;

    % --- TOP LEFT ---
    idx = NN(1,1,n);
    % xi - F, eta - F
    diag = 0.5*(J(idx).^-1).*[-3*x_eta(idx)/dxi+3*x_xi(idx)/deta, ...
            4*x_eta(idx)/dxi,-x_eta(idx)/dxi, ...
            -4*x_xi(idx)/deta,x_xi(idx)/deta];
    diag_idx = [0,1,2,n,2*n];

    A(idx,:) = 0;   % Clear rows of A corresponding to the boundary
    A(idx,idx+diag_idx) = diag;
    
    % Generate Ey by applying the finite difference matrix
    Ey = -A*u;

    % Reshape gradient components to match meshgrid variables
    Ex = flip(reshape(Ex,n,m)');
    Ey = flip(reshape(Ey,n,m)');
    
    % Node numbering function, i -> xi, j -> eta, n = total xi nodes
    function NN = NN(i,j,n)
        NN = n*(j-1) + i;
    end

end