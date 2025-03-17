% -------------------------------------------------------------------
% Patrick Heng
% 08 Mar. 2025
% Upwind flux scheme for the convection term in the convection-
% diffusion equations over a uniformly spaced rectangular grid.
% -------------------------------------------------------------------

function F = upwind_flux(X,Y,U,V,dx,dy,n,nodes)
    
    %NN = @(i,j) n*(j-1) + i;
    
    uw = U(X-dx/2,Y); ue = U(X+dx/2,Y);
    vs = V(X,Y-dy/2); vn = V(X,Y+dy/2);

    uw = reshape(uw',nodes,1); ue = reshape(ue',nodes,1); 
    vs = reshape(vs',nodes,1); vn = reshape(vn',nodes,1);
    
    W = -0.5*(uw+abs(uw))/dx;
    E = 0.5*(ue-abs(ue))/dx;
    S = -0.5*(vs+abs(vs))/dy;
    N = 0.5*(vn-abs(vn))/dy;
    P = 0.5*(ue+abs(ue)-uw+abs(uw))/dx + 0.5*(vn+abs(vn)-vs+abs(vs))/dy;

    diag = [S,W,P,E,N];

    F = spdiags(diag,[-n,-1,0,1,n],nodes,nodes+1);
    F = F(:,1:nodes);

    
end