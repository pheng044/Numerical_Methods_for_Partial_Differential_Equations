% -------------------------------------------------------------------
% Patrick Heng
% 22 Mar. 2025
% Diffusive flux scheme for the convection term in the convection-
% diffusion equations over a uniformly spaced rectangular grid.
% -------------------------------------------------------------------

function G = diffusive_flux(D,dx,dy,n,m,nodes)
    % Node numbering function
    NN = @(i,j) n*(j-1) + i;
    
    % Fill diffusive flux stencil
    diag = [1/dy^2,1/dx^2,-2/dx^2-2/dy^2,1/dx^2,1/dy^2];
    diag = repmat(diag,nodes,1);
    diag(1:n:nodes,2) = 0;
    diag(n:n:nodes-1,4) = 0;

    % Place stencil coefficients on diagonals
    G = spdiags(diag,[-n,-1,0,1,n],nodes,nodes+1);
    G = G(:,1:nodes);
    
    % Boundary Corrections for periodic boundary conditions
    top = [NN(1:n,1)',NN(1:n,m)',ones(n,1)/dy^2];
    btm = [NN(1:n,m)',NN(1:n,1)',ones(n,1)/dy^2];
    left = [NN(1,1:m)',NN(n,1:m)',ones(m,1)/dx^2];
    right = [NN(n,1:m)',NN(1,1:m)',ones(m,1)/dx^2];
    
    % Form sparse diagonal matrices for corrections
    top = sparse(top(:,1),top(:,2),top(:,3),nodes,nodes);
    btm = sparse(btm(:,1),btm(:,2),btm(:,3),nodes,nodes);
    left = sparse(left(:,1),left(:,2),left(:,3),nodes,nodes);
    right = sparse(right(:,1),right(:,2),right(:,3),nodes,nodes);
    
    % Apply corrections 
    G = D*(G + top + btm + left + right);
    
end