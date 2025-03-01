% -------------------------------------------------------------------
% Patrick Heng
% 7 Feb. 2025
% Function to generate bilinear map coefficients from a unit square
% to a general affine grid.
% Inputs: Xi, Eta: Meshgrid variables defined over the unit square
%         x, y: Vector of x and y coordinates of verticies of transformed
%               coordinates. Specified in the order (d,a,b,c)          
%               c ----- b
%               |       |
%               |       |
%               d ----- a
% Outputs: X, Y: Meshgrid variables defined over the physical domain
%          coeffs_x, coeffs_y vectors of transformation coefficients
%          specified such that
%               x(xi,eta) = coeffs_x'*[xi.*eta; xi; eta; 1]
%          and similarly for y(xi,eta)
% -------------------------------------------------------------------

function [X,Y,coeffs_x,coeffs_y] = bilinear_map(Xi,Eta,x,y)
    
    % Make sure x, y are column vectors
    x = reshape(x,4,1);
    y = reshape(y,4,1);
    
    % Get size of mesh
    sz = size(Xi);
    N = numel(Xi);

    % Vectorize the mesh, calculate the points needed for mapping
    xi = reshape(Xi,N,1);
    eta =  reshape(Eta,N,1);
    xi_eta = reshape(Xi.*Eta,N,1);
    
    % Place vectorized mesh into matrix
    bilin =  [xi_eta, xi, eta, ones(N,1)]';
    
    % Transfer matrix to find mapping coefficients
    trans_mtx = [0,0,0,1; 0,1,0,1; 1,1,1,1; 0,0,1,1];
    
    % Solve for mapping coefficients
    coeffs_x = trans_mtx\x;
    coeffs_y = trans_mtx\y;

    % Calculate new coordinates
    coords(1,:) = coeffs_x'*bilin;
    coords(2,:) = coeffs_y'*bilin;

    % Return X and Y back as in the same style as meshgrid
    X = reshape(coords(1,:),sz(1),sz(2));
    Y = reshape(coords(2,:),sz(1),sz(2));

end