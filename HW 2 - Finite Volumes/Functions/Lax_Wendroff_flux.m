% -------------------------------------------------------------------
% Patrick Heng
% 08 Mar. 2025
% Lax-Wendroff flux scheme for the convection term in the convection-
% diffusion equations over a uniformly spaced rectangular grid.
% -------------------------------------------------------------------

function [Fx,Fy] = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,m,nodes,periodic)
    % Node numbering function
    NN = @(i,j) n*(j-1) + i;

    % Evaluate face velocities for fluxes
    uw = U(X-dx/2,Y); ue = U(X+dx/2,Y);
    vs = V(X,Y-dy/2); vn = V(X,Y+dy/2);

    % Vectorize velocities
    uw = reshape(uw',nodes,1); ue = reshape(ue',nodes,1); 
    vs = reshape(vs',nodes,1); vn = reshape(vn',nodes,1);
    
    % Ensure the left boundary velocity is the same as the right 
    % boundary velocity
    uw(NN(1,1:m)) = ue(NN(n,1:m));

    % Stencil coefficients
    W = -0.5*(uw+(uw.^2)*dt/dx)/dx;
    E = 0.5*(ue-(ue.^2)*dt/dx)/dx;
    S = -0.5*(vs+(vs.^2)*dt/dy)/dy;
    N = 0.5*(vn-(vn.^2)*dt/dy)/dy;
    P = 0.5*(ue+(ue.^2)*dt/dx-uw+(uw.^2)*dt/dx)/dx;

    % Left, right periodic condition
    diag = zeros(nodes,2);
    diag(NN(n,1:m),1) = E(NN(n,1:m));
    diag(NN(1,1:m),2) = W(NN(1,1:m));

    % Sparse corrector matrix
    Fx_BC = spdiags(diag,[-n+1,n-1],nodes,nodes+1);
    Fx_BC = Fx_BC(:,1:nodes);

    % Avoid double counting
    E(NN(n,1:m)) = 0;
    W(NN(1,1:m)) = 0;
    diag = [W,P,E];

    % Form Fx matrix with corrections
    Fx = spdiags(diag,[-1,0,1],nodes,nodes+1);
    Fx = Fx(:,1:nodes);
    Fx = Fx + Fx_BC;
    
    % Py coefficients
    P = 0.5*(vn+(vn.^2)*dt/dy-vs+(vs.^2)*dt/dy)/dy;
    
    if periodic == false
        % Bottom correction
        P(NN(1:n,1)) = P(NN(1:n,1)) - S(NN(1:n,1));
        S(NN(1:n,1)) = 0;
        % Top correction
        P(NN(1:n,m)) = P(NN(1:n,m)) + N(NN(1:n,m));
        N(NN(1:n,m)) = 0;
    end
    
    % Sparse correction matrix
    diag = [S,P,N];
    Fy = spdiags(diag,[-n,0,n],nodes,nodes+1);
    Fy = Fy(:,1:nodes);

    if periodic == true
        % Top, bottom periodic condition
        diag = [N,S];
        Fy_BC = spdiags(diag,[-nodes+n,nodes-n],nodes,nodes+1);
        Fy_BC = Fy_BC(:,1:nodes);
        Fy = Fy + Fy_BC;
    end
    
end
