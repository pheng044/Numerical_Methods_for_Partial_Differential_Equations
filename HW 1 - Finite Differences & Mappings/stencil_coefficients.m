% -------------------------------------------------------------------
% Patrick Heng
% 7 Feb 2025
% 2D stencil coefficients for mapped Laplacian and derivatives for 
% finite differences over a uniform computational domain.
% -------------------------------------------------------------------

function coeffs = stencil_coefficients(J,a,b,c,d,e,dxi,deta,nodes)

    % Compute stencil coefficients
    SW = b/(2*dxi*deta)./(J.^2);
    S = -(c/deta^2 - d/(2*deta))./(J.^2);
    SE = -SW;
    W = -(a/dxi^2 - e/(2*dxi))./(J.^2);
    P = 2*(a/dxi^2+c/deta^2)./(J.^2);
    E = -(a/dxi^2 + e/(2*dxi))./(J.^2);
    NW = -SW;
    N = -(c/deta^2 + d/(2*deta))./(J.^2);
    NE = SW;
    
    % Vectorize stencil coefficients, vertiacally stack rows 
    % of contant eta
    SW = reshape(SW',nodes,1);
    S = reshape(S',nodes,1);
    SE = reshape(SE',nodes,1);
    W = reshape(W',nodes,1);
    P = reshape(P',nodes,1);
    E = reshape(E',nodes,1);
    NW = reshape(NW',nodes,1);
    N = reshape(N',nodes,1);
    NE = reshape(NE',nodes,1);
    
    % Return coefficients in convenient form for using spdiags
    coeffs = [SW,S,SE,W,P,E,NW,N,NE];
    
end
