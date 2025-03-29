% -------------------------------------------------------------------
% Patrick Heng
% 08 Mar. 2025
% Upwind flux scheme for the convection term in the convection-
% diffusion equations over a uniformly spaced rectangular grid.
% -------------------------------------------------------------------

function F = upwind_flux(X,Y,U,V,dx,dy,n,m,nodes,periodic)
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
    W = -0.5*(uw+abs(uw))/dx;
    E = 0.5*(ue-abs(ue))/dx;
    S = -0.5*(vs+abs(vs))/dy;
    N = 0.5*(vn-abs(vn))/dy;
    P = 0.5*(ue+abs(ue)-uw+abs(uw))/dx + 0.5*(vn+abs(vn)-vs+abs(vs))/dy;
    
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

    if periodic == false
        % Bottom correction
        P(NN(1:n,1)) = P(NN(1:n,1)) - S(NN(1:n,1));
        S(NN(1:n,1)) = 0;
        % Top correction
        P(NN(1:n,m)) = P(NN(1:n,m)) + N(NN(1:n,m));
        N(NN(1:n,m)) = 0;
        % Form flux matrix
        diag = [S,W,P,E,N];
        F = spdiags(diag,[-n,-1,0,1,n],nodes,nodes+1);
        F = F(:,1:nodes) + Fx_BC;
    else
        % Form flux matrix
        diag = [S,W,P,E,N];
        F = spdiags(diag,[-n,-1,0,1,n],nodes,nodes+1);
        F = F(:,1:nodes);
        % Sparse corrector matrix
        diag = [N,S];
        F_BC = spdiags(diag,[-nodes+n,nodes-n],nodes,nodes+1);
        F_BC = F_BC(:,1:nodes);
        F = F + F_BC + Fx_BC;
    end

end
