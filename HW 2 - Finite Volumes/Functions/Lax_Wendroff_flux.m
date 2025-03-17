% -------------------------------------------------------------------
% Patrick Heng
% 08 Mar. 2025
% Lax-Wendroff flux scheme for the convection term in the convection-
% diffusion equations over a uniformly spaced rectangular grid.
% -------------------------------------------------------------------

function [Fx,Fy] = Lax_Wendroff_flux(X,Y,U,V,dx,dy,dt,n,nodes)
    
    %NN = @(i,j) n*(j-1) + i;
    
    uw = U(X-dx/2,Y); ue = U(X+dx/2,Y);
    vs = V(X,Y-dy/2); vn = V(X,Y+dy/2);

    uw = reshape(uw',nodes,1); ue = reshape(ue',nodes,1); 
    vs = reshape(vs',nodes,1); vn = reshape(vn',nodes,1);
    
    W = -0.5*(uw+uw.^2*dt/dx)/dx;
    E = 0.5*(ue-ue.^2*dt/dx)/dx;
    S = -0.5*(vs+vs.^2*dt/dy)/dy;
    N = 0.5*(vn-vn.^2*dt/dy)/dy;
    P = 0.5*(ue+ue.^2*dt/dx-uw+uw.^2*dt/dx)/dx;

    diag = [W,P,E];

    Fx = spdiags(diag,[-1,0,1],nodes,nodes+1);
    Fx = Fx(:,1:nodes);

    P = 0.5*(vn+vn.^2*dt/dy-vs+vs.^2*dt/dy)/dy;
    
    diag = [S,P,N];

    Fy = spdiags(diag,[-n,0,n],nodes,nodes+1);
    Fy = Fy(:,1:nodes);

    
end